"""
CrossGeneEx GUI — Cross-species L2/3 IT gene expression explorer
Run: python3 app_gui.py
"""

import sys
import io
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("QtAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.lines as mlines

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QSplitter,
    QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QListWidget, QListWidgetItem,
    QTreeWidget, QTreeWidgetItem, QTabWidget, QScrollArea,
    QTableWidget, QTableWidgetItem, QHeaderView,
    QFileDialog, QMessageBox, QProgressDialog, QFrame,
    QAbstractItemView, QSizePolicy,
)
from PyQt6.QtCore import Qt, QSize
from PyQt6.QtGui import QFont, QColor

from goatools.obo_parser import GODag


WORKDIR   = Path(__file__).parent / "data"
SP_ORDER  = ["Human", "RhesusM", "Mouse"]
SP_LABELS = {"Human": "Human", "RhesusM": "NHP", "Mouse": "Mouse"}
SP_COLORS = {"Human": "#C1121F", "RhesusM": "#0077B6", "Mouse": "#1B4332"}


# ── Data ───────────────────────────────────────────────────────────────────────
class DataStore:
    def __init__(self):
        print("Loading expression matrix …", flush=True)
        self.expr = pd.read_csv(
            WORKDIR / "AllGenes_L23IT_10xOnly_species_comparison_strict_log2cpm.csv.gz",
            index_col=0,
        )
        meta = (
            pd.read_csv(
                WORKDIR / "AllGenes_L23IT_10xOnly_species_comparison_strict_sample_metadata.csv"
            ).set_index("pseudobulk_sample")
        )
        self.meta    = meta.loc[self.expr.columns]
        self.kw_all  = pd.read_csv(WORKDIR / "AllGenes_L23IT_10xOnly_kruskal_results.csv").set_index("gene")
        self.kw_sig  = pd.read_csv(WORKDIR / "AllGenes_L23IT_10xOnly_kruskal_significant.csv").set_index("gene")
        self.syngo   = pd.read_csv(WORKDIR / "AllGenes_L23IT_syngo_annotations.csv")

        self.kw_all_d  = self.kw_all.to_dict(orient="index")
        self.kw_sig_d  = self.kw_sig.to_dict(orient="index")
        self.sp_samples = {sp: self.meta.index[self.meta["species"] == sp].tolist() for sp in SP_ORDER}
        self.all_genes  = sorted(self.expr.index.tolist())

        self.syngo_terms = sorted(self.syngo["syngo_term"].unique())
        self.term_genes  = self.syngo.groupby("syngo_term")["human_gene"].apply(list).to_dict()
        self.children, self.term_order = self._build_hierarchy()
        print("Ready.", flush=True)

    def _build_hierarchy(self):
        print("Building SynGO hierarchy …", flush=True)
        term_set = set(self.syngo_terms)
        godag = GODag(
            str(WORKDIR / "go-basic.obo"),
            optional_attrs={"relationship"},
            load_obsolete=False,
            prt=io.StringIO(),
        )
        name2go = {t.name.lower(): t for t in godag.values()}

        parent_map: dict = {}
        for name in self.syngo_terms:
            go = name2go.get(name.lower())
            if go:
                is_a    = [p.name for p in go.parents if p.name in term_set]
                part_of = [t.name for t in go.relationship.get("part_of", []) if t.name in term_set]
                parent_map[name] = list(set(is_a + part_of))
            else:
                parent_map[name] = []

        children: dict = defaultdict(list)
        for term, parents in parent_map.items():
            for p in parents:
                children[p].append(term)
        children = {k: sorted(v) for k, v in children.items()}

        roots      = sorted(t for t in self.syngo_terms if not parent_map[t])
        term_order: list = []
        visited:    set  = set()

        def dfs(node, depth):
            if node in visited: return
            visited.add(node)
            term_order.append((node, depth))
            for child in children.get(node, []):
                dfs(child, depth + 1)

        for r in roots: dfs(r, 0)
        for t in self.syngo_terms:
            if t not in visited: term_order.append((t, 0))

        return dict(children), term_order

    def get_desc_genes(self, term: str) -> list:
        genes = set(self.term_genes.get(term, []))
        stack = list(self.children.get(term, []))
        while stack:
            t = stack.pop()
            genes.update(self.term_genes.get(t, []))
            stack.extend(self.children.get(t, []))
        return sorted(genes)


# ── Plot & table logic ─────────────────────────────────────────────────────────
def pstars(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    return "*"

def fmt_p(p: float) -> str:
    if np.isnan(p): return "N/A"
    if p < 0.0001:  return "< 0.0001"
    return f"{p:.4f}"


def make_figure(gene: str, ds: DataStore) -> Figure:
    vals  = {sp: ds.expr.loc[gene, ds.sp_samples[sp]].values.astype(float) for sp in SP_ORDER}
    means = {sp: float(v.mean()) for sp, v in vals.items()}
    kw_p  = float(ds.kw_all_d.get(gene, {}).get("kw_pval", np.nan))

    fig = Figure(figsize=(5, 5.5), dpi=100)
    ax  = fig.add_subplot(111)
    x   = np.arange(3)

    ax.bar(x, [means[sp] for sp in SP_ORDER], width=0.55,
           color=[SP_COLORS[sp] for sp in SP_ORDER], alpha=0.80, zorder=2)

    rng = np.random.default_rng(42)
    for i, sp in enumerate(SP_ORDER):
        jitter = rng.uniform(-0.14, 0.14, len(vals[sp]))
        ax.scatter(i + jitter, vals[sp], color="white",
                   edgecolors=SP_COLORS[sp], s=90, lw=2.0, zorder=4)

    all_vals = np.concatenate(list(vals.values()))
    y_max    = float(all_vals.max())
    gap      = max(y_max * 0.12, 0.10)

    if (not np.isnan(kw_p)) and kw_p < 0.05 and gene in ds.kw_sig_d:
        d = ds.kw_sig_d[gene]
        raw = [
            (0, 1, d.get("p_Human_NHP",   np.nan)),
            (0, 2, d.get("p_Human_Mouse",  np.nan)),
            (1, 2, d.get("p_NHP_Mouse",    np.nan)),
        ]
        sig_pairs = [(i, j, float(p)) for i, j, p in raw
                     if p and not np.isnan(float(p)) and float(p) < 0.05]
        if sig_pairs:
            y0 = y_max + gap * 0.5
            for level, (xi, xj, p) in enumerate(sig_pairs):
                y = y0 + level * gap * 1.4; tip = gap * 0.30
                ax.plot([x[xi], x[xi], x[xj], x[xj]], [y, y+tip, y+tip, y],
                        lw=1.8, color="black", zorder=5)
                ax.text((x[xi]+x[xj])/2, y+tip*1.1, pstars(p),
                        ha="center", va="bottom", fontsize=24)
            ax.set_ylim(0, y0 + len(sig_pairs)*gap*1.4 + gap*0.8)
        else:
            ax.set_ylim(0, y_max * 1.22)
    else:
        ax.set_ylim(0, y_max * 1.22)

    ax.set_xticks(x)
    ax.set_xticklabels([SP_LABELS[sp] for sp in SP_ORDER], fontsize=24)
    ax.set_ylabel("log₂(CPM + 1)", fontsize=20)
    ax.set_title(gene, fontsize=26, fontweight="bold")
    ax.tick_params(axis="y", labelsize=18)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(-0.5, 2.5)

    kw_handle = mlines.Line2D([], [], color="none", linestyle="none",
                               label=f"KW p = {fmt_p(kw_p)}")
    ax.legend(handles=[kw_handle], loc="upper right", fontsize=18,
              frameon=True, framealpha=0.8, edgecolor="#cccccc",
              handlelength=0, handletextpad=0)
    fig.tight_layout()
    return fig


def build_table_df(genes: list, ds: DataStore) -> pd.DataFrame:
    sp_key = {v: k for k, v in SP_LABELS.items()}
    rows   = []
    for gene in genes:
        if gene not in ds.expr.index: continue
        vals  = {sp: ds.expr.loc[gene, ds.sp_samples[sp]].values.astype(float) for sp in SP_ORDER}
        means = {sp: float(v.mean()) for sp, v in vals.items()}
        kw_p  = float(ds.kw_all_d.get(gene, {}).get("kw_pval", np.nan))
        posthoc = ""
        if (not np.isnan(kw_p)) and kw_p < 0.05 and gene in ds.kw_sig_d:
            d = ds.kw_sig_d[gene]
            parts = []
            for col, la, lb in [("p_Human_NHP","Human","NHP"),
                                 ("p_Human_Mouse","Human","Mouse"),
                                 ("p_NHP_Mouse","NHP","Mouse")]:
                p = d.get(col, np.nan)
                if p and not np.isnan(float(p)) and float(p) < 0.05:
                    ka, kb = sp_key[la], sp_key[lb]
                    hi, lo = (la, lb) if means[ka] >= means[kb] else (lb, la)
                    parts.append(f"{hi}>{lo} (p={float(p):.3f})")
            posthoc = "; ".join(parts)
        rows.append({
            "Gene":                  gene,
            "Human mean":            round(means["Human"],   4),
            "NHP mean":              round(means["RhesusM"], 4),
            "Mouse mean":            round(means["Mouse"],   4),
            "KW p-value":            round(kw_p, 6) if not np.isnan(kw_p) else "",
            "Post-hoc (sig. pairs)": posthoc,
        })
    return pd.DataFrame(rows)


# ── Main window ────────────────────────────────────────────────────────────────
class MainWindow(QMainWindow):
    def __init__(self, ds: DataStore):
        super().__init__()
        self.ds       = ds
        self.figs:    dict = {}
        self.table_df      = None

        self.setWindowTitle("CrossGeneEx — Cross-species L2/3 IT Gene Expression Explorer")
        self.resize(1700, 980)

        central = QWidget()
        self.setCentralWidget(central)
        root_lay = QHBoxLayout(central)
        root_lay.setContentsMargins(6, 6, 6, 6)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        root_lay.addWidget(splitter)
        splitter.addWidget(self._build_left_panel())
        splitter.addWidget(self._build_right_panel())
        splitter.setSizes([460, 1240])
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

    # ── Left panel ──────────────────────────────────────────────────────────────
    def _build_left_panel(self) -> QWidget:
        panel = QWidget()
        panel.setMinimumWidth(440)
        panel.setMaximumWidth(520)
        lay = QVBoxLayout(panel)
        lay.setSpacing(8)

        self.tabs = QTabWidget()
        lay.addWidget(self.tabs)

        # Tab 1 — Gene name search
        gene_tab = QWidget()
        gt_lay   = QVBoxLayout(gene_tab)
        gt_lay.addWidget(QLabel("Search gene symbol:"))

        self.gene_search = QLineEdit()
        self.gene_search.setPlaceholderText("e.g. DLG4, SHANK3 …")
        self.gene_search.textChanged.connect(self._filter_gene_list)
        gt_lay.addWidget(self.gene_search)

        self.gene_count_lbl = QLabel(f"{len(self.ds.all_genes):,} genes shown")
        gt_lay.addWidget(self.gene_count_lbl)

        self.gene_list = QListWidget()
        self._populate_gene_list(self.ds.all_genes)
        gt_lay.addWidget(self.gene_list)

        gt_btn_row = QHBoxLayout()
        for label, state in [("Check All", True), ("Uncheck All", False)]:
            btn = QPushButton(label)
            btn.clicked.connect(lambda _, s=state: self._set_all_checked(self.gene_list, s))
            gt_btn_row.addWidget(btn)
        gt_lay.addLayout(gt_btn_row)
        self.tabs.addTab(gene_tab, "Search by Gene Name")

        # Tab 2 — SynGO browser
        syngo_tab = QWidget()
        st_lay    = QVBoxLayout(syngo_tab)

        syngo_search = QLineEdit()
        syngo_search.setPlaceholderText("Filter terms …")
        syngo_search.textChanged.connect(self._filter_syngo_tree)
        st_lay.addWidget(syngo_search)

        self.syngo_tree = QTreeWidget()
        self.syngo_tree.setHeaderHidden(True)
        self.syngo_tree.itemSelectionChanged.connect(self._on_syngo_term_selected)
        self._populate_syngo_tree()

        self.syngo_gene_lbl  = QLabel("Select a term above:")
        self.syngo_gene_list = QListWidget()
        self.syngo_gene_list.setMaximumHeight(230)

        sg_btn_row = QHBoxLayout()
        for label, state in [("Check All", True), ("Uncheck All", False)]:
            btn = QPushButton(label)
            btn.clicked.connect(lambda _, s=state: self._set_all_checked(self.syngo_gene_list, s))
            sg_btn_row.addWidget(btn)

        # Vertical splitter for tree vs gene list
        v_split = QSplitter(Qt.Orientation.Vertical)
        v_split.addWidget(self.syngo_tree)
        gene_box = QWidget()
        gene_box_lay = QVBoxLayout(gene_box)
        gene_box_lay.setContentsMargins(0, 0, 0, 0)
        gene_box_lay.addWidget(self.syngo_gene_lbl)
        gene_box_lay.addWidget(self.syngo_gene_list)
        gene_box_lay.addLayout(sg_btn_row)
        v_split.addWidget(gene_box)
        v_split.setSizes([400, 280])
        st_lay.addWidget(v_split)
        self.tabs.addTab(syngo_tab, "Browse by SynGO Term")

        # Generate button
        sep = QFrame(); sep.setFrameShape(QFrame.Shape.HLine)
        lay.addWidget(sep)
        self.gen_btn = QPushButton("▶   Generate Plots")
        self.gen_btn.setFixedHeight(48)
        self.gen_btn.setStyleSheet(
            "QPushButton { background:#1a6fc4; color:white; font-weight:bold; border-radius:6px; }"
            "QPushButton:hover { background:#1457a0; }"
        )
        self.gen_btn.clicked.connect(self._on_generate)
        lay.addWidget(self.gen_btn)

        return panel

    # ── Right panel ─────────────────────────────────────────────────────────────
    def _build_right_panel(self) -> QScrollArea:
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)

        ph = QLabel("Select genes in the left panel, then click  ▶ Generate Plots")
        ph.setAlignment(Qt.AlignmentFlag.AlignCenter)
        ph.setStyleSheet("color: #888888;")
        self.scroll.setWidget(ph)
        return self.scroll

    # ── Gene list helpers ────────────────────────────────────────────────────────
    def _populate_gene_list(self, genes: list) -> None:
        self.gene_list.clear()
        for gene in genes:
            item = QListWidgetItem(gene)
            item.setCheckState(Qt.CheckState.Unchecked)
            self.gene_list.addItem(item)

    def _filter_gene_list(self, text: str) -> None:
        text = text.strip().lower()
        filtered = [g for g in self.ds.all_genes if text in g.lower()] if text else self.ds.all_genes
        self._populate_gene_list(filtered)
        self.gene_count_lbl.setText(f"{len(filtered):,} genes shown")

    def _set_all_checked(self, lw: QListWidget, state: bool) -> None:
        s = Qt.CheckState.Checked if state else Qt.CheckState.Unchecked
        for i in range(lw.count()):
            lw.item(i).setCheckState(s)

    def _get_selected_genes(self) -> list:
        lw = self.gene_list if self.tabs.currentIndex() == 0 else self.syngo_gene_list
        return [
            lw.item(i).text()
            for i in range(lw.count())
            if lw.item(i).checkState() == Qt.CheckState.Checked
        ]

    # ── SynGO tree helpers ───────────────────────────────────────────────────────
    def _populate_syngo_tree(self, filter_text: str = "") -> None:
        self.syngo_tree.clear()
        self._tree_items: dict = {}
        visited: set = set()

        def add_item(term: str, parent=None):
            if term in visited: return
            visited.add(term)
            n     = len(self.ds.get_desc_genes(term))
            label = f"{term}  [{n}]"
            item  = (QTreeWidgetItem(self.syngo_tree, [label]) if parent is None
                     else QTreeWidgetItem(parent, [label]))
            item.setData(0, Qt.ItemDataRole.UserRole, term)
            self._tree_items[term] = item
            for child in self.ds.children.get(term, []):
                add_item(child, item)

        ft = filter_text.strip().lower()
        for term, depth in self.ds.term_order:
            if depth == 0:
                if ft and ft not in term.lower():
                    # include root only if it or a descendant matches
                    desc = self.ds.get_desc_genes(term)
                    if not any(ft in t.lower() for t in
                               [term] + list(self.ds.children.get(term, []))):
                        if ft not in term.lower():
                            continue
                add_item(term)

        self.syngo_tree.expandToDepth(0)

    def _filter_syngo_tree(self, text: str) -> None:
        self._populate_syngo_tree(filter_text=text)

    def _on_syngo_term_selected(self) -> None:
        items = self.syngo_tree.selectedItems()
        if not items: return
        term  = items[0].data(0, Qt.ItemDataRole.UserRole)
        if not term: return
        genes = self.ds.get_desc_genes(term)
        self.syngo_gene_lbl.setText(f"{len(genes)} gene(s) under '{term}':")
        self.syngo_gene_list.clear()
        for gene in genes:
            item = QListWidgetItem(gene)
            item.setCheckState(Qt.CheckState.Checked)
            self.syngo_gene_list.addItem(item)

    # ── Generate ─────────────────────────────────────────────────────────────────
    def _on_generate(self) -> None:
        raw   = self._get_selected_genes()
        genes = [g for g in raw if g in self.ds.expr.index]
        if not genes:
            QMessageBox.warning(self, "No genes selected",
                                "Check at least one gene before clicking Generate Plots.")
            return

        self.figs.clear()
        self.table_df = None

        prog = QProgressDialog("Generating plots…", None, 0, len(genes), self)
        prog.setWindowModality(Qt.WindowModality.WindowModal)
        prog.setMinimumDuration(0)
        prog.setWindowTitle("CrossGeneEx")
        prog.setValue(0)
        QApplication.processEvents()

        figs_ordered = []
        for i, gene in enumerate(genes):
            prog.setValue(i)
            prog.setLabelText(f"Plotting {gene} ({i+1} / {len(genes)}) …")
            QApplication.processEvents()
            fig = make_figure(gene, self.ds)
            self.figs[gene] = fig
            figs_ordered.append((gene, fig))

        prog.setValue(len(genes))
        self.table_df = build_table_df(genes, self.ds)
        self._render_right(figs_ordered)
        self.setWindowTitle(f"CrossGeneEx — {len(genes)} gene(s) plotted")

    # ── Render right panel ───────────────────────────────────────────────────────
    def _render_right(self, figs_ordered: list) -> None:
        content = QWidget()
        v_lay   = QVBoxLayout(content)
        v_lay.setContentsMargins(12, 12, 12, 12)
        v_lay.setSpacing(14)

        # Plot grid — 3 columns
        NCOLS     = 3
        plot_grid = QGridLayout()
        plot_grid.setSpacing(18)

        for idx, (gene, fig) in enumerate(figs_ordered):
            row, col = divmod(idx, NCOLS)
            cell     = QWidget()
            c_lay    = QVBoxLayout(cell)
            c_lay.setSpacing(4)
            c_lay.setContentsMargins(0, 0, 0, 0)

            canvas = FigureCanvas(fig)
            canvas.setFixedSize(QSize(500, 550))
            c_lay.addWidget(canvas)

            dl_btn = QPushButton(f"⬇  Download  {gene}.png")
            dl_btn.setFixedHeight(38)
            dl_btn.clicked.connect(lambda _checked, g=gene, f=fig: self._download_plot(g, f))
            c_lay.addWidget(dl_btn)
            plot_grid.addWidget(cell, row, col)

        v_lay.addLayout(plot_grid)

        # Divider
        sep = QFrame(); sep.setFrameShape(QFrame.Shape.HLine); sep.setFrameShadow(QFrame.Shadow.Sunken)
        v_lay.addWidget(sep)

        # Table header
        tbl_hdr = QLabel("Summary Table")
        tbl_hdr.setStyleSheet("font-weight: bold;")
        v_lay.addWidget(tbl_hdr)

        v_lay.addWidget(self._make_qtable(self.table_df))

        dl_csv = QPushButton("⬇  Download Table (CSV)")
        dl_csv.setFixedHeight(38)
        dl_csv.clicked.connect(self._download_table)
        v_lay.addWidget(dl_csv)
        v_lay.addStretch()

        self.scroll.setWidget(content)

    # ── Table widget ─────────────────────────────────────────────────────────────
    def _make_qtable(self, df: pd.DataFrame) -> QTableWidget:
        tbl = QTableWidget(len(df), len(df.columns))
        tbl.setHorizontalHeaderLabels(df.columns.tolist())
        tbl.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        tbl.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        tbl.setAlternatingRowColors(True)
        tbl.horizontalHeader().setStretchLastSection(True)
        tbl.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        tbl.verticalHeader().setVisible(False)
        for r in range(len(df)):
            for c, val in enumerate(df.iloc[r]):
                tbl.setItem(r, c, QTableWidgetItem("" if val is None else str(val)))
        tbl.resizeColumnsToContents()
        tbl.setMinimumHeight(min(500, 36 + 34 * len(df)))
        return tbl

    # ── Downloads ─────────────────────────────────────────────────────────────────
    def _download_plot(self, gene: str, fig: Figure) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, f"Save {gene} plot", f"{gene}_expression.png", "PNG Image (*.png)"
        )
        if path:
            fig.savefig(path, dpi=160, bbox_inches="tight", facecolor="white")

    def _download_table(self) -> None:
        if self.table_df is None: return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save summary table", "CrossGeneEx_table.csv", "CSV file (*.csv)"
        )
        if path:
            self.table_df.to_csv(path, index=False)


# ── Entry point ────────────────────────────────────────────────────────────────
def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")

    font = QFont("Helvetica Neue", 14)
    app.setFont(font)

    app.setStyleSheet("""
        /* Base: always black text on light background */
        QMainWindow, QWidget          { background: #f5f5f5; color: #111111; }

        QTabWidget::pane              { border: 1px solid #cccccc; background: #ffffff; }
        QTabBar::tab                  { padding: 8px 18px; font-size: 13pt;
                                        color: #111111; background: #e0e0e0; }
        QTabBar::tab:selected         { background: #1a6fc4; color: white;
                                        border-radius: 4px 4px 0 0; }
        QTabBar::tab:hover:!selected  { background: #cccccc; color: #111111; }

        QLabel                        { font-size: 13pt; color: #111111; }

        QLineEdit                     { font-size: 13pt; padding: 5px;
                                        border: 1px solid #bbbbbb; border-radius: 4px;
                                        background: white; color: #111111; }

        QListWidget                   { font-size: 13pt; border: 1px solid #cccccc;
                                        background: white; color: #111111;
                                        alternate-background-color: #f0f4ff; }
        QListWidget::item             { color: #111111; }
        QListWidget::item:hover       { background: #dde8f8; color: #111111; }
        QListWidget::item:selected    { background: #1a6fc4; color: white; }

        QTreeWidget                   { font-size: 12pt; border: 1px solid #cccccc;
                                        background: white; color: #111111;
                                        alternate-background-color: #f0f4ff; }
        QTreeWidget::item             { color: #111111; }
        QTreeWidget::item:hover       { background: #dde8f8; color: #111111; }
        QTreeWidget::item:selected    { background: #1a6fc4; color: white; }

        QPushButton                   { font-size: 13pt; color: #111111;
                                        border: 1px solid #aaaaaa; border-radius: 5px;
                                        padding: 5px 12px; background: #e8e8e8; }
        QPushButton:hover             { background: #d0d8e8; color: #111111; }

        QTableWidget                  { font-size: 12pt; color: #111111;
                                        gridline-color: #dddddd;
                                        alternate-background-color: #f5f8ff; }
        QTableWidget::item            { color: #111111; }
        QTableWidget::item:selected   { background: #1a6fc4; color: white; }

        QHeaderView::section          { font-size: 13pt; font-weight: bold; color: #111111;
                                        background: #e0e6f0; padding: 6px; border: none; }

        QScrollArea                   { border: none; }

        QProgressDialog               { color: #111111; }
        QProgressDialog QLabel        { font-size: 13pt; color: #111111; }
        QMessageBox QLabel            { font-size: 13pt; color: #111111; }
        QFileDialog                   { color: #111111; }

        QSplitter::handle             { background: #cccccc; }
    """)

    ds  = DataStore()
    win = MainWindow(ds)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
