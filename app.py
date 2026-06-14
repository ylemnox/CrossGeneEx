"""
CrossGeneEx — Interactive cross-species L2/3 IT gene expression explorer
Run:  streamlit run app.py
"""

import io
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
import streamlit as st
from goatools.obo_parser import GODag

# ── Paths ──────────────────────────────────────────────────────────────────────
WORKDIR = Path(__file__).parent / "data"

SP_ORDER  = ["Human", "RhesusM", "Mouse"]
SP_LABELS = {"Human": "Human", "RhesusM": "NHP", "Mouse": "Mouse"}
SP_COLORS = {"Human": "#C1121F", "RhesusM": "#0077B6", "Mouse": "#1B4332"}

# ── Data loaders ───────────────────────────────────────────────────────────────
@st.cache_data(show_spinner="Loading expression matrix …")
def load_data():
    expr = pd.read_csv(
        WORKDIR / "AllGenes_L23IT_10xOnly_species_comparison_strict_log2cpm.csv.gz",
        index_col=0,
    )
    meta = (
        pd.read_csv(
            WORKDIR / "AllGenes_L23IT_10xOnly_species_comparison_strict_sample_metadata.csv"
        )
        .set_index("pseudobulk_sample")
    )
    meta   = meta.loc[expr.columns]
    kw_all = pd.read_csv(WORKDIR / "AllGenes_L23IT_10xOnly_kruskal_results.csv").set_index("gene")
    kw_sig = pd.read_csv(WORKDIR / "AllGenes_L23IT_10xOnly_kruskal_significant.csv").set_index("gene")
    syngo  = pd.read_csv(WORKDIR / "AllGenes_L23IT_syngo_annotations.csv")
    return expr, meta, kw_all, kw_sig, syngo


@st.cache_data(show_spinner="Building SynGO term hierarchy …")
def build_hierarchy(syngo_terms_tuple):
    """
    Match SynGO term names against go-basic.obo, build parent→children tree.
    Returns: children dict, term_order list of (term, depth) in DFS pre-order.
    """
    terms    = list(syngo_terms_tuple)
    term_set = set(terms)

    godag = GODag(
        str(WORKDIR / "go-basic.obo"),
        optional_attrs={"relationship"},
        load_obsolete=False,
        prt=io.StringIO(),           # suppress GO-loading output
    )
    name2go = {t.name.lower(): t for t in godag.values()}

    # Direct parents of each SynGO term that are also in our SynGO set.
    # Use both IS_A and PART_OF — SynGO's anatomical structure uses PART_OF
    # (e.g. presynaptic active zone PART_OF presynapse, not IS_A presynapse).
    parent_map: dict[str, list[str]] = {}
    for name in terms:
        go_entry = name2go.get(name.lower())
        if go_entry:
            is_a    = [p.name for p in go_entry.parents if p.name in term_set]
            part_of = [t.name for t in go_entry.relationship.get("part_of", [])
                       if t.name in term_set]
            parent_map[name] = list(set(is_a + part_of))
        else:
            parent_map[name] = []

    # Children map (sorted alphabetically at each level)
    children: dict[str, list[str]] = defaultdict(list)
    for term, parents in parent_map.items():
        for p in parents:
            children[p].append(term)
    children = {k: sorted(v) for k, v in children.items()}

    # DFS pre-order for UI display
    roots      = sorted(t for t in terms if not parent_map[t])
    term_order: list[tuple[str, int]] = []
    visited: set[str] = set()

    def dfs(node: str, depth: int) -> None:
        if node in visited:
            return
        visited.add(node)
        term_order.append((node, depth))
        for child in children.get(node, []):
            dfs(child, depth + 1)

    for root in roots:
        dfs(root, 0)
    for t in terms:           # catch any nodes not reached via roots
        if t not in visited:
            term_order.append((t, 0))

    return children, term_order


def get_desc_genes(
    term: str,
    children: dict[str, list[str]],
    term_genes: dict[str, list[str]],
) -> set[str]:
    """All genes annotated to this term or any of its descendants."""
    genes: set[str] = set(term_genes.get(term, []))
    stack = list(children.get(term, []))
    while stack:
        t = stack.pop()
        genes.update(term_genes.get(t, []))
        stack.extend(children.get(t, []))
    return genes


# ── Stat helpers ───────────────────────────────────────────────────────────────
def pstars(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    return "*"


def fmt_p(p: float) -> str:
    if np.isnan(p):  return "N/A"
    if p < 0.0001:   return "< 0.0001"
    return f"{p:.4f}"


# ── Bar plot per gene ──────────────────────────────────────────────────────────
def make_plot(
    gene: str,
    expr: pd.DataFrame,
    sp_samples: dict[str, list[str]],
    kw_all_d: dict,
    kw_sig_d: dict,
) -> plt.Figure:
    vals  = {sp: expr.loc[gene, sp_samples[sp]].values.astype(float) for sp in SP_ORDER}
    means = {sp: float(v.mean()) for sp, v in vals.items()}
    kw_p  = float(kw_all_d.get(gene, {}).get("kw_pval", np.nan))

    fig, ax = plt.subplots(figsize=(6, 6.5))
    x = np.arange(3)

    # Bars
    ax.bar(
        x,
        [means[sp] for sp in SP_ORDER],
        width=0.55,
        color=[SP_COLORS[sp] for sp in SP_ORDER],
        alpha=0.80,
        zorder=2,
    )

    # Individual donor dots with jitter
    rng = np.random.default_rng(42)
    for i, sp in enumerate(SP_ORDER):
        n      = len(vals[sp])
        jitter = rng.uniform(-0.14, 0.14, n)
        ax.scatter(
            i + jitter, vals[sp],
            color="white", edgecolors=SP_COLORS[sp],
            s=90, lw=2.0, zorder=4,
        )

    # Dunn post-hoc brackets
    all_vals = np.concatenate(list(vals.values()))
    y_max    = float(all_vals.max())
    gap      = max(y_max * 0.12, 0.10)

    if (not np.isnan(kw_p)) and kw_p < 0.05 and gene in kw_sig_d:
        d = kw_sig_d[gene]
        raw_pairs = [
            (0, 1, d.get("p_Human_NHP",   np.nan)),
            (0, 2, d.get("p_Human_Mouse",  np.nan)),
            (1, 2, d.get("p_NHP_Mouse",    np.nan)),
        ]
        sig_pairs = [
            (i, j, float(p)) for i, j, p in raw_pairs
            if p is not None and not np.isnan(float(p)) and float(p) < 0.05
        ]

        if sig_pairs:
            y0 = y_max + gap * 0.5
            for level, (xi, xj, p) in enumerate(sig_pairs):
                y   = y0 + level * gap * 1.4
                tip = gap * 0.30
                ax.plot(
                    [x[xi], x[xi], x[xj], x[xj]],
                    [y, y + tip, y + tip, y],
                    lw=1.8, color="black", zorder=5,
                )
                ax.text(
                    (x[xi] + x[xj]) / 2,
                    y + tip * 1.1,
                    pstars(p),
                    ha="center", va="bottom", fontsize=24, zorder=5,
                )
            top = y0 + len(sig_pairs) * gap * 1.4 + gap * 0.8
            ax.set_ylim(0, top)
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

    # KW p-value as legend entry (no line/marker, just text)
    kw_handle = mlines.Line2D(
        [], [], color="none", linestyle="none",
        label=f"KW p = {fmt_p(kw_p)}",
    )
    ax.legend(
        handles=[kw_handle],
        loc="upper right",
        fontsize=18,
        frameon=True,
        framealpha=0.8,
        edgecolor="#cccccc",
        handlelength=0,
        handletextpad=0,
    )

    plt.tight_layout()
    return fig


def fig_to_png(fig: plt.Figure) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    return buf.read()


# ── Summary table ──────────────────────────────────────────────────────────────
def build_table(
    genes: list[str],
    expr: pd.DataFrame,
    sp_samples: dict[str, list[str]],
    kw_all_d: dict,
    kw_sig_d: dict,
) -> pd.DataFrame:
    sp_key = {v: k for k, v in SP_LABELS.items()}   # "NHP" → "RhesusM"
    rows   = []

    for gene in genes:
        if gene not in expr.index:
            continue
        vals  = {sp: expr.loc[gene, sp_samples[sp]].values.astype(float) for sp in SP_ORDER}
        means = {sp: float(v.mean()) for sp, v in vals.items()}
        kw_p  = float(kw_all_d.get(gene, {}).get("kw_pval", np.nan))

        posthoc = ""
        if (not np.isnan(kw_p)) and kw_p < 0.05 and gene in kw_sig_d:
            d = kw_sig_d[gene]
            pairs = [
                ("p_Human_NHP",   "Human", "NHP"),
                ("p_Human_Mouse", "Human", "Mouse"),
                ("p_NHP_Mouse",   "NHP",   "Mouse"),
            ]
            parts = []
            for col, la, lb in pairs:
                p = d.get(col, np.nan)
                if p is None or np.isnan(float(p)) or float(p) >= 0.05:
                    continue
                ka, kb     = sp_key[la], sp_key[lb]
                hi, lo     = (la, lb) if means[ka] >= means[kb] else (lb, la)
                parts.append(f"{hi}>{lo} (p={float(p):.3f})")
            posthoc = "; ".join(parts)

        rows.append({
            "Gene":                  gene,
            "Human mean":            round(means["Human"],   4),
            "NHP mean":              round(means["RhesusM"], 4),
            "Mouse mean":            round(means["Mouse"],   4),
            "KW p-value":            round(kw_p, 6) if not np.isnan(kw_p) else None,
            "Post-hoc (sig. pairs)": posthoc,
        })

    return pd.DataFrame(rows)


# ══════════════════════════════════════════════════════════════════════════════
# Streamlit UI
# ══════════════════════════════════════════════════════════════════════════════
st.set_page_config(page_title="CrossGeneEx", page_icon="🧬", layout="wide")

# ── Global UI font scaling (2.5×) ──────────────────────────────────────────────
st.markdown("""
<style>
/* Increase base font for all Streamlit text elements */
html, body, [class*="css"] {
    font-size: 22px !important;
}

/* Widget labels (selectbox, multiselect, radio) */
label,
.stRadio label,
[data-testid="stWidgetLabel"] p,
[data-testid="stWidgetLabel"] span {
    font-size: 22px !important;
    line-height: 1.5 !important;
}

/* Sidebar header */
[data-testid="stSidebarContent"] h1,
[data-testid="stSidebarContent"] h2,
[data-testid="stSidebarContent"] h3 {
    font-size: 28px !important;
}

/* General paragraph / markdown text */
[data-testid="stMarkdownContainer"] p,
[data-testid="stMarkdownContainer"] li,
[data-testid="stCaptionContainer"] p,
.stCaption {
    font-size: 20px !important;
    line-height: 1.5 !important;
}

/* Info / warning banners */
[data-testid="stAlert"] p {
    font-size: 20px !important;
}

/* Dropdown options and multiselect pills */
[data-baseweb="select"] *,
[data-baseweb="tag"] *,
[data-baseweb="menu"] * {
    font-size: 20px !important;
}

/* Buttons */
.stButton > button,
.stDownloadButton > button {
    font-size: 22px !important;
    padding: 0.55rem 1.2rem !important;
}

/* Divider label and subheader */
h2, h3 {
    font-size: 26px !important;
}

/* Dataframe table text */
[data-testid="stDataFrame"] * {
    font-size: 18px !important;
}
</style>
""", unsafe_allow_html=True)

st.title("🧬 CrossGeneEx")
st.caption(
    "Cross-species L2/3 IT gene expression explorer · Human (n=4) · NHP (n=3) · Mouse (n=10)"
)

# Load all data
expr, meta, kw_all, kw_sig, syngo_df = load_data()

sp_samples = {sp: meta.index[meta["species"] == sp].tolist() for sp in SP_ORDER}
all_genes  = sorted(expr.index.tolist())

# Convert KW dataframes to dicts for fast lookup
kw_all_d: dict = kw_all.to_dict(orient="index")
kw_sig_d: dict = kw_sig.to_dict(orient="index")

# SynGO structures
syngo_terms = tuple(sorted(syngo_df["syngo_term"].unique()))
term_genes  = syngo_df.groupby("syngo_term")["human_gene"].apply(list).to_dict()

children, term_order = build_hierarchy(syngo_terms)

# Precompute descendant gene counts for each SynGO term
term_counts = {t: len(get_desc_genes(t, children, term_genes)) for t, _ in term_order}

# Session state: persist plots across rerenders
if "genes_to_plot" not in st.session_state:
    st.session_state.genes_to_plot = []

# ── Sidebar ────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Gene Selection")

    mode = st.radio(
        "Selection mode",
        ["Search by gene name", "Browse by SynGO term"],
        label_visibility="collapsed",
    )

    if mode == "Search by gene name":
        selected: list[str] = st.multiselect(
            f"Select genes  ({len(all_genes):,} available)",
            options=all_genes,
            help="Start typing a gene symbol to filter.",
        )

    else:  # Browse by SynGO term
        # Build dropdown labels with visual hierarchy indentation
        labels: list[str] = []
        names:  list[str] = []
        for term_name, depth in term_order:
            prefix = "  " * depth + ("└─ " if depth > 0 else "")
            n      = term_counts.get(term_name, 0)
            labels.append(f"{prefix}{term_name}  [{n}]")
            names.append(term_name)

        chosen_label = st.selectbox(
            "SynGO term (type to search)",
            options=labels,
            help="[N] = total genes including sub-terms",
        )
        chosen_term = names[labels.index(chosen_label)]

        pool = sorted(get_desc_genes(chosen_term, children, term_genes))
        st.caption(f"**{len(pool)}** gene(s) under *{chosen_term}* (incl. sub-terms)")

        default = pool[:20] if len(pool) > 20 else pool
        if len(pool) > 20:
            st.info(f"Showing first 20 of {len(pool)} by default — adjust below.")

        selected = st.multiselect(
            "Genes to plot",
            options=pool,
            default=default,
        )

    st.divider()

    if st.button("▶  Generate Plots", type="primary", use_container_width=True):
        valid   = [g for g in selected if g in expr.index]
        missing = [g for g in selected if g not in expr.index]
        if missing:
            st.warning(f"Not found in dataset: {', '.join(missing)}")
        st.session_state.genes_to_plot = valid

    if st.session_state.genes_to_plot:
        if st.button("✕  Clear", use_container_width=True):
            st.session_state.genes_to_plot = []
            st.rerun()

# ── Main area ──────────────────────────────────────────────────────────────────
genes = st.session_state.genes_to_plot

if not genes:
    st.info("Select genes in the sidebar and click **▶ Generate Plots**.")
    st.stop()

st.subheader(f"Expression plots — {len(genes)} gene(s)")

ncols = min(3, len(genes))
cols  = st.columns(ncols)

for idx, gene in enumerate(genes):
    with cols[idx % ncols]:
        fig = make_plot(gene, expr, sp_samples, kw_all_d, kw_sig_d)
        st.pyplot(fig, use_container_width=True)
        png = fig_to_png(fig)
        plt.close(fig)
        st.download_button(
            label="⬇ Download PNG",
            data=png,
            file_name=f"{gene}_expression.png",
            mime="image/png",
            key=f"dl_png_{gene}_{idx}",
            use_container_width=True,
        )

st.divider()
st.subheader("Summary table")

table = build_table(genes, expr, sp_samples, kw_all_d, kw_sig_d)
st.dataframe(table, use_container_width=True, hide_index=True)

st.download_button(
    label="⬇ Download table (CSV)",
    data=table.to_csv(index=False).encode(),
    file_name="CrossGeneEx_table.csv",
    mime="text/csv",
)
