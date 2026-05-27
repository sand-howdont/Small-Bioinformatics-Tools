#!/usr/bin/env python3
import argparse
from collections import defaultdict, Counter, deque

import matplotlib.pyplot as plt
import matplotlib as mpl


BASE_Y = {
    "A": 3,
    "T": 2,
    "C": 1,
    "G": 0,
}


def parse_gfa(gfa_file):
    """
    Parse GFA file.

    S line: node sequence
    L line: graph edge
    P line: path through graph

    Returns:
        nodes: dict, node_id -> base
        edges: set of (source, target)
        paths: dict, path_name -> list of node_id
    """
    nodes = {}
    edges = set()
    paths = {}

    with open(gfa_file, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            record_type = fields[0]

            if record_type == "S":
                node_id = fields[1]
                seq = fields[2].upper()
                base = seq[0] if seq else "N"
                nodes[node_id] = base

            elif record_type == "L":
                source = fields[1]
                target = fields[3]
                edges.add((source, target))

            elif record_type == "P":
                path_name = fields[1]
                raw_nodes = fields[2].split(",")

                path_nodes = []
                for item in raw_nodes:
                    node = item.rstrip("+-")
                    path_nodes.append(node)

                paths[path_name] = path_nodes

    return nodes, edges, paths


def count_edge_support(edges, paths):
    """
    Count how many paths support each edge.

    If an edge appears in many P paths, it is considered strong.
    """
    support = Counter()

    for path_name, path_nodes in paths.items():
        path_edges = zip(path_nodes[:-1], path_nodes[1:])

        for edge in path_edges:
            support[edge] += 1

    for edge in edges:
        support.setdefault(edge, 0)

    return support


def topological_x_layout(nodes, edges):
    """
    Assign x coordinate to each node according to graph order.

    x = longest distance from graph start.
    """
    graph = defaultdict(list)
    indegree = defaultdict(int)

    for node in nodes:
        indegree[node] = 0

    for source, target in edges:
        if source in nodes and target in nodes:
            graph[source].append(target)
            indegree[target] += 1

    queue = deque([node for node in nodes if indegree[node] == 0])
    x_pos = {node: 0 for node in nodes}

    visited = 0

    while queue:
        node = queue.popleft()
        visited += 1

        for nxt in graph[node]:
            x_pos[nxt] = max(x_pos.get(nxt, 0), x_pos[node] + 1)
            indegree[nxt] -= 1

            if indegree[nxt] == 0:
                queue.append(nxt)

    if visited < len(nodes):
        remaining = [n for n in nodes if n not in x_pos]
        start_x = max(x_pos.values()) + 1 if x_pos else 0

        for i, node in enumerate(remaining):
            x_pos[node] = start_x + i

    return x_pos


def get_node_positions(nodes, edges):
    """
    Get x/y coordinates.

    y is fixed by base:
        A, T, C, G
    """
    x_raw = topological_x_layout(nodes, edges)

    same_slot_count = defaultdict(int)
    pos = {}

    for node in sorted(nodes, key=lambda x: int(x) if x.isdigit() else x):
        base = nodes[node]

        if base not in BASE_Y:
            continue

        x = x_raw[node]
        y = BASE_Y[base]

        key = (x, y)
        offset_index = same_slot_count[key]
        same_slot_count[key] += 1

        jitter = (offset_index - 0.5) * 0.12 if offset_index > 0 else 0

        pos[node] = (x + jitter, y)

    return pos


def plot_gfa(nodes, edges, paths, output_file):
    support = count_edge_support(edges, paths)
    pos = get_node_positions(nodes, edges)

    max_support = max(support.values()) if support else 1
    if max_support == 0:
        max_support = 1

    max_x = max([p[0] for p in pos.values()]) if pos else 1
    fig_width = max(12, min(max_x * 0.18, 80))
    fig_height = 5

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Low support = blue, high support = red
    cmap = mpl.cm.coolwarm
    norm = mpl.colors.Normalize(vmin=0, vmax=max_support)

    # Draw only base row labels on the left.
    # No horizontal dashed guide lines.
    for base, y in BASE_Y.items():
        ax.text(
            -1.5,
            y,
            base,
            fontsize=14,
            fontweight="bold",
            color="black",
            ha="right",
            va="center"
        )

    # Draw edges
    for source, target in edges:
        if source not in pos or target not in pos:
            continue

        x1, y1 = pos[source]
        x2, y2 = pos[target]

        edge_support = support.get((source, target), 0)
        color = cmap(norm(edge_support))

        # Stronger support -> thicker line
        linewidth = 0.4 + 4.0 * edge_support / max_support

        ax.annotate(
            "",
            xy=(x2, y2),
            xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle="->",
                color=color,
                linewidth=linewidth,
                alpha=0.9,
                shrinkA=2,
                shrinkB=2,
                mutation_scale=8
            ),
            zorder=1
        )

    # Draw base letters only, no circles
    for node, (x, y) in pos.items():
        base = nodes[node]

        ax.text(
            x,
            y,
            base,
            fontsize=9,
            fontweight="bold",
            color="black",
            ha="center",
            va="center",
            zorder=4
        )

    # Colorbar for edge support
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = plt.colorbar(
        sm,
        ax=ax,
        orientation="horizontal",
        pad=0.18,
        fraction=0.06
    )
    cbar.set_label("Edge support: number of paths containing this connection")

    ax.set_ylim(-0.8, 3.8)
    ax.set_xlim(-2, max_x + 2)

    # Keep only y-axis labels A/T/C/G, remove y-axis line and tick marks
    ax.set_yticks([BASE_Y[b] for b in ["G", "C", "T", "A"]])
    ax.set_yticklabels(["G", "C", "T", "A"], fontweight="bold", color="black")
    ax.tick_params(axis="y", length=0)

    # Remove y-axis spine and other unnecessary spines
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_xlabel("Graph position")
    ax.set_title("Partial order alignment graph from GFA")

    plt.subplots_adjust(left=0.04, right=0.98, top=0.88, bottom=0.22)

    plt.savefig(
        output_file,
        dpi=200,
        bbox_inches="tight",
        pad_inches=0.3
    )

    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Plot a GFA partial order alignment graph with four base rows: A/T/C/G."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input GFA file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output figure file, e.g. graph.png, graph.pdf, graph.svg"
    )

    args = parser.parse_args()

    nodes, edges, paths = parse_gfa(args.input)

    if not nodes:
        raise ValueError("No S records found in the GFA file.")

    if not edges:
        raise ValueError("No L records found in the GFA file.")

    plot_gfa(nodes, edges, paths, args.output)


if __name__ == "__main__":
    main()
