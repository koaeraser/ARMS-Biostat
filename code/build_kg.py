"""
Session 2: Knowledge Graph Construction
========================================
Builds a NetworkX MultiDiGraph from trials_data.csv and drug_target_map.csv.
Computes pairwise inter-trial similarity matrices (drug Jaccard, target Jaccard,
population similarity, composite). Generates an interactive pyvis visualization.

Outputs:
    data/kg_graph.gpickle
    data/similarity_matrices.pkl
    figures/kg_visualization.html
"""

import os
import pickle
import warnings

import numpy as np
import pandas as pd
import networkx as nx
from pyvis.network import Network

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
FIG_DIR = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(FIG_DIR, exist_ok=True)

TRIALS_CSV = os.path.join(DATA_DIR, "trials_data.csv")
DRUG_TARGET_CSV = os.path.join(DATA_DIR, "drug_target_map.csv")
GRAPH_OUT = os.path.join(DATA_DIR, "kg_graph.gpickle")
SIM_OUT = os.path.join(DATA_DIR, "similarity_matrices.pkl")
VIS_OUT = os.path.join(FIG_DIR, "kg_visualization.html")


# ===================================================================
# 1. Load data
# ===================================================================
def load_data():
    """Load trials and drug-target data."""
    trials = pd.read_csv(TRIALS_CSV)
    drug_target = pd.read_csv(DRUG_TARGET_CSV)
    return trials, drug_target


def parse_targets(raw: str) -> str:
    """Normalise target string to a canonical short name.

    drug_target_map has entries like 'PSMB5 (20S proteasome)' or
    'CRBN (cereblon)'.  We keep only the short symbol before the
    parenthetical for node identity, so two PIs map to the same target.
    """
    return raw.strip().split("(")[0].strip().split("/")[0].strip()


# ===================================================================
# 2. Build the knowledge graph
# ===================================================================
def build_graph(trials: pd.DataFrame, drug_target: pd.DataFrame) -> nx.MultiDiGraph:
    """Construct a MultiDiGraph with trial, drug, target, endpoint, and
    population nodes connected by typed edges.

    Node attributes
    ---------------
    All nodes: node_type ∈ {trial, drug, target, endpoint, population}

    Trial nodes carry the full feature vector X(v_t) from the plan:
        p_hat, n, phase, asct_required, median_age, frac_high_risk,
        frac_iss3, mrd_threshold, triplet_or_quad, trial_name, arm_label

    Drug nodes: drug_class, mechanism
    Target nodes: (name only)
    Endpoint: (name only)
    Population: (name only)

    Edge types (relation)
    ---------------------
    uses_drug   : trial → drug
    targets     : drug  → target
    measures    : trial → endpoint
    enrolls     : trial → population
    """
    G = nx.MultiDiGraph()

    # --- Drug & target nodes + edges ---
    target_set = set()
    for _, row in drug_target.iterrows():
        drug = row["drug_name"]
        G.add_node(drug, node_type="drug",
                   drug_class=row["drug_class"],
                   mechanism=row["mechanism"])
        tgt = parse_targets(row["targets"])
        target_set.add(tgt)
        G.add_node(tgt, node_type="target")
        G.add_edge(drug, tgt, relation="targets")

    # --- Endpoint node (single: MRD-negativity) ---
    G.add_node("MRD-negativity", node_type="endpoint")

    # --- Population nodes ---
    for pop in ("TE", "TIE", "mixed"):
        G.add_node(pop, node_type="population")

    # --- Trial nodes + edges ---
    for _, row in trials.iterrows():
        tid = row["trial_id"]
        G.add_node(tid, node_type="trial",
                   trial_name=row["trial_name"],
                   arm_label=row["arm_label"],
                   p_hat=row["p_hat"],
                   n=int(row["n"]),
                   y=int(row["y"]),
                   phase=int(row["phase"]),
                   asct_required=int(row["asct_required"]),
                   median_age=float(row["median_age"]),
                   frac_high_risk=float(row["frac_high_risk"]),
                   frac_iss3=float(row["frac_iss3"]),
                   mrd_threshold=row["mrd_threshold"],
                   triplet_or_quad=row["triplet_or_quad"],
                   population=row["population"],
                   nct_number=row.get("nct_number", ""))

        # uses_drug edges
        drugs = [d.strip() for d in row["drugs"].split(",")]
        for drug in drugs:
            G.add_edge(tid, drug, relation="uses_drug")

        # measures edge
        G.add_edge(tid, "MRD-negativity", relation="measures")

        # enrolls edge
        G.add_edge(tid, row["population"], relation="enrolls")

    return G


# ===================================================================
# 3. Compute similarity matrices
# ===================================================================
def get_trial_ids(G: nx.MultiDiGraph) -> list:
    """Return sorted list of trial node IDs."""
    return sorted([n for n, d in G.nodes(data=True) if d["node_type"] == "trial"])


def drug_neighbours(G: nx.MultiDiGraph, trial_id: str) -> set:
    """D(t) = {drug nodes reachable via uses_drug from trial t}."""
    return {v for _, v, d in G.out_edges(trial_id, data=True)
            if d["relation"] == "uses_drug"}


def target_neighbours(G: nx.MultiDiGraph, trial_id: str) -> set:
    """T(t) = 2-hop target neighbourhood via uses_drug → targets."""
    targets = set()
    for drug in drug_neighbours(G, trial_id):
        for _, tgt, d in G.out_edges(drug, data=True):
            if d["relation"] == "targets":
                targets.add(tgt)
    return targets


def jaccard(set_a: set, set_b: set) -> float:
    """Jaccard similarity.  Returns 0 if both sets empty."""
    if not set_a and not set_b:
        return 0.0
    return len(set_a & set_b) / len(set_a | set_b)


def compute_drug_jaccard(G: nx.MultiDiGraph, trial_ids: list) -> np.ndarray:
    """40×40 drug-level Jaccard matrix."""
    n = len(trial_ids)
    drug_sets = {t: drug_neighbours(G, t) for t in trial_ids}
    J = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            val = jaccard(drug_sets[trial_ids[i]], drug_sets[trial_ids[j]])
            J[i, j] = J[j, i] = val
    return J


def compute_target_jaccard(G: nx.MultiDiGraph, trial_ids: list) -> np.ndarray:
    """40×40 target-level Jaccard matrix (2-hop)."""
    n = len(trial_ids)
    tgt_sets = {t: target_neighbours(G, t) for t in trial_ids}
    J = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            val = jaccard(tgt_sets[trial_ids[i]], tgt_sets[trial_ids[j]])
            J[i, j] = J[j, i] = val
    return J


def _pop_feature_vector(G: nx.MultiDiGraph, trial_id: str) -> np.ndarray:
    """Encode population characteristics as a numeric vector.

    Components (all scaled to [0, 1]):
        0: TE indicator (1 if TE, 0.5 if mixed, 0 if TIE)
        1: frac_high_risk  (already [0, 1])
        2: median_age / 100 (reasonable normaliser for age 50-80)
        3: frac_iss3        (already [0, 1])
        4: mrd_threshold indicator (1 if 1e-5, 0 if 1e-4)
        5: asct_required (0 or 1)
    """
    attrs = G.nodes[trial_id]
    pop_map = {"TE": 1.0, "mixed": 0.5, "TIE": 0.0}
    te = pop_map.get(attrs["population"], 0.5)
    mrd = 1.0 if attrs["mrd_threshold"] == "1e-5" else 0.0
    return np.array([
        te,
        attrs["frac_high_risk"],
        attrs["median_age"] / 100.0,
        attrs["frac_iss3"],
        mrd,
        float(attrs["asct_required"]),
    ])


def compute_population_similarity(G: nx.MultiDiGraph, trial_ids: list) -> np.ndarray:
    """Population similarity: 1 − ||x_pop(i) − x_pop(j)||_1 / C.

    C is the maximum possible L1 distance across the 5-d feature vector
    (each component ranges [0, 1], so C = 5).
    """
    n = len(trial_ids)
    vecs = {t: _pop_feature_vector(G, t) for t in trial_ids}
    C = 6.0  # max L1 for 6 components each in [0, 1]
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            l1 = np.sum(np.abs(vecs[trial_ids[i]] - vecs[trial_ids[j]]))
            val = 1.0 - l1 / C
            S[i, j] = S[j, i] = val
    return S


def compute_composite_similarity(
    J_drug: np.ndarray,
    J_target: np.ndarray,
    S_pop: np.ndarray,
    alpha: float = 1 / 3,
    beta: float = 1 / 3,
    gamma: float = 1 / 3,
) -> np.ndarray:
    """Composite similarity: s = α·J_drug + β·J_target + γ·pop_sim."""
    assert abs(alpha + beta + gamma - 1.0) < 1e-9, "Weights must sum to 1"
    return alpha * J_drug + beta * J_target + gamma * S_pop


# ===================================================================
# 4. Interactive visualisation (pyvis)
# ===================================================================

# Colour palette per node type
_COLORS = {
    "trial": "#4e79a7",
    "drug": "#f28e2b",
    "target": "#e15759",
    "endpoint": "#76b7b2",
    "population": "#59a14f",
}

_SHAPES = {
    "trial": "dot",
    "drug": "diamond",
    "target": "triangle",
    "endpoint": "star",
    "population": "square",
}

_EDGE_COLORS = {
    "uses_drug": "#888888",
    "targets": "#e15759",
    "measures": "#76b7b2",
    "enrolls": "#59a14f",
}


def _trial_label(attrs: dict) -> str:
    """Short label for a trial node."""
    name = attrs.get("trial_name", "")
    arm = attrs.get("arm_label", "")
    p = attrs.get("p_hat", 0)
    return f"{name}\n{arm}\np̂={p:.2f}"


def visualise(G: nx.MultiDiGraph, out_path: str):
    """Create an interactive HTML visualisation using pyvis."""
    net = Network(height="900px", width="100%", directed=True,
                  bgcolor="#ffffff", font_color="#333333")
    net.barnes_hut(gravity=-8000, central_gravity=0.3,
                   spring_length=200, spring_strength=0.01)

    # Add nodes
    for node, attrs in G.nodes(data=True):
        ntype = attrs.get("node_type", "trial")
        color = _COLORS.get(ntype, "#999999")
        shape = _SHAPES.get(ntype, "dot")

        if ntype == "trial":
            label = _trial_label(attrs)
            size = 15 + attrs.get("n", 50) / 30  # scale by sample size
            title = (f"<b>{attrs.get('trial_name','')} — {attrs.get('arm_label','')}</b><br>"
                     f"n={attrs.get('n','')}, y={attrs.get('y','')}, p̂={attrs.get('p_hat',0):.3f}<br>"
                     f"Phase {attrs.get('phase','')}, ASCT={attrs.get('asct_required','')}<br>"
                     f"Pop: {attrs.get('population','')}, MRD: {attrs.get('mrd_threshold','')}<br>"
                     f"HR: {attrs.get('frac_high_risk',0):.0%}, ISS3: {attrs.get('frac_iss3',0):.0%}")
        elif ntype == "drug":
            label = node
            size = 20
            title = f"<b>{node}</b><br>Class: {attrs.get('drug_class','')}<br>{attrs.get('mechanism','')}"
        elif ntype == "target":
            label = node
            size = 18
            title = f"<b>Target: {node}</b>"
        else:
            label = node
            size = 16
            title = f"<b>{node}</b> ({ntype})"

        net.add_node(node, label=label, color=color, shape=shape,
                     size=size, title=title)

    # Add edges
    for u, v, d in G.edges(data=True):
        rel = d.get("relation", "")
        color = _EDGE_COLORS.get(rel, "#cccccc")
        net.add_edge(u, v, title=rel, color=color, arrows="to", width=1.5)

    net.save_graph(out_path)
    print(f"  Visualisation saved → {out_path}")


# ===================================================================
# 5. Summary statistics
# ===================================================================
def print_graph_summary(G: nx.MultiDiGraph):
    """Print basic graph statistics."""
    type_counts = {}
    for _, d in G.nodes(data=True):
        t = d.get("node_type", "unknown")
        type_counts[t] = type_counts.get(t, 0) + 1
    edge_counts = {}
    for _, _, d in G.edges(data=True):
        r = d.get("relation", "unknown")
        edge_counts[r] = edge_counts.get(r, 0) + 1

    print("\n=== Knowledge Graph Summary ===")
    print(f"  Total nodes : {G.number_of_nodes()}")
    for t, c in sorted(type_counts.items()):
        print(f"    {t:12s}: {c}")
    print(f"  Total edges : {G.number_of_edges()}")
    for r, c in sorted(edge_counts.items()):
        print(f"    {r:12s}: {c}")


def print_similarity_summary(trial_ids, J_drug, J_target, S_pop, S_comp):
    """Print summary statistics of similarity matrices."""
    n = len(trial_ids)
    # Extract upper triangle (excluding diagonal)
    mask = np.triu_indices(n, k=1)

    print("\n=== Similarity Matrix Summary (off-diagonal) ===")
    for name, M in [("Drug Jaccard", J_drug), ("Target Jaccard", J_target),
                    ("Population", S_pop), ("Composite", S_comp)]:
        vals = M[mask]
        print(f"  {name:16s}: mean={vals.mean():.3f}, "
              f"median={np.median(vals):.3f}, "
              f"min={vals.min():.3f}, max={vals.max():.3f}")

    # Most similar pair (composite)
    np.fill_diagonal(S_comp, -1)
    idx = np.unravel_index(S_comp.argmax(), S_comp.shape)
    np.fill_diagonal(S_comp, 1.0)
    print(f"\n  Most similar pair (composite): "
          f"{trial_ids[idx[0]]} ↔ {trial_ids[idx[1]]} = {S_comp[idx]:.3f}")

    # Least similar pair (composite)
    np.fill_diagonal(S_comp, 2)
    idx = np.unravel_index(S_comp.argmin(), S_comp.shape)
    np.fill_diagonal(S_comp, 1.0)
    print(f"  Least similar pair (composite): "
          f"{trial_ids[idx[0]]} ↔ {trial_ids[idx[1]]} = {S_comp[idx]:.3f}")


# ===================================================================
# Main
# ===================================================================
def main():
    print("Session 2: Knowledge Graph Construction")
    print("=" * 50)

    # 1. Load
    trials, drug_target = load_data()
    print(f"Loaded {len(trials)} trial arms, {len(drug_target)} drugs")

    # 2. Build graph
    G = build_graph(trials, drug_target)
    print_graph_summary(G)

    # 3. Similarity matrices
    trial_ids = get_trial_ids(G)
    J_drug = compute_drug_jaccard(G, trial_ids)
    J_target = compute_target_jaccard(G, trial_ids)
    S_pop = compute_population_similarity(G, trial_ids)
    S_comp = compute_composite_similarity(J_drug, J_target, S_pop)

    print_similarity_summary(trial_ids, J_drug, J_target, S_pop, S_comp)

    # 4. Save graph
    with open(GRAPH_OUT, "wb") as f:
        pickle.dump(G, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"\n  Graph saved → {GRAPH_OUT}")

    # 5. Save similarity matrices
    sim_data = {
        "trial_ids": trial_ids,
        "J_drug": J_drug,
        "J_target": J_target,
        "S_pop": S_pop,
        "S_composite": S_comp,
        "weights": {"alpha": 1 / 3, "beta": 1 / 3, "gamma": 1 / 3},
    }
    with open(SIM_OUT, "wb") as f:
        pickle.dump(sim_data, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"  Similarity matrices saved → {SIM_OUT}")

    # 6. Visualisation
    visualise(G, VIS_OUT)

    print("\nSession 2 complete.")


if __name__ == "__main__":
    main()
