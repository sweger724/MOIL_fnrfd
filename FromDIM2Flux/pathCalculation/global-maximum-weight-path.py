#!/usr/bin/env python

"""Compute maximum weight path.

This program computes global maximum weight paths using flux-space
graphs with flux-based edge weights (see [1]).

There is a modification to the algorithm: we extend the graph by
adding two extra nodes.  One node that acts as a source and another
node that acts as a sink.  The source is connected to the reactant
milestones through edges weighted according to the fluxes of the
reactants.  The product milestones are connected to the sink via edges
with zero weight.

Reference:

[1] Viswanath, S., Kreuzer, S. M., Cardenas, A. E., & Elber,
R. (2013). Analyzing milestoning networks for molecular kinetics:
definitions, algorithms, and examples. The Journal of Chemical
Physics, 139(17), 174105. doi:10.1063/1.4827495

"""

from __future__ import print_function

import networkx as nx


inf = float('inf')

# Meta-reactant and meta-product
meta_reactant = 0
meta_product = inf


def read_special_nodes(milestones_file):
    """Reads file indicating the indices of reactants and products.

    The file must contain lines of the form:

      milestone-index reactant
      ...
      milestone-index product
      ...
    """

    reactants = []
    products = []

    with file(milestones_file, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
        for line in lines:
            info = line.split()
            idx = int(info[0])
            if line.endswith('reactant'):
                reactants.append(idx)
            elif line.endswith('product'):
                products.append(idx)

    return reactants, products


def read_connectivity_graph(matrix_file, flux_file,
                            reactants, products):
    """Creates flux-based graph.

    The matrix file is a text file containing tuples of the form:

      number-of-entries
      index-i index-j value
      index-i index-j value
      index-i index-j value
      ...

    where the third column, value, is strictly positive.

    The flux file is a text file containing the values of the
    stationary flux vector at each milestone:

      q-1
      q-2
      ...

    For each triplet in the matrix file, we create an edge in the
    graph.  The weight of the edge is determined by the value of the
    stationary flux corresponding to the milestone labelled index-j.

    """

    G = nx.MultiDiGraph()

    weights = []
    with file(flux_file, 'r') as qf:
        weights = [float(line.strip()) for line in qf.readlines()]

    # Add meta nodes to the graph.
    for j in reactants:
        G.add_weighted_edges_from([(meta_reactant, j, weights[j-1])])
    for i in products:
        G.add_weighted_edges_from([(i, meta_product, 0)])

    with file(matrix_file, 'r') as mf:
        num_milestones = int(mf.readline().strip())
        assert(num_milestones == len(weights))

        lines = [line.strip() for line in mf.readlines()]
        for line in lines:
            a, b, c = line.split()
            i = int(a)
            j = int(b)
            if i in products:
                continue
            G.add_weighted_edges_from([(i, j, weights[j-1])])

    return G


def compute_max_weight(G, start=meta_reactant):
    """Compute maximum weight paths.

    This function annotates each node of the graph with its maximum
    weight from the starting node and a link to its predecessor.

    """
    def extract_max(labels):
        max_weight = -1
        max_label = None
        for label in labels:
            if G.node[label]['max_weight'] > max_weight:
                max_weight = G.node[label]['max_weight']
                max_label = label

        if max_label is not None:
            labels.remove(max_label)

        return max_label, max_weight

    labels = []
    G.node[start]['max_weight'] = inf
    G.node[start]['predecessor'] = None
    G.node[start]['bottleneck'] = None
    labels.append(start)
    for label in G.nodes_iter():
        if label is not start:
            G.node[label]['max_weight'] = -1
            G.node[label]['predecessor'] = None
            G.node[label]['bottleneck'] = None
            labels.append(label)

    while len(labels) > 0:
        # Find out which node has maximal weight.
        current, weight = extract_max(labels)

        for neighbor in G[current]:
            m = min(weight, G.edge[current][neighbor][0]['weight'])
            if G.node[neighbor]['max_weight'] < m:
                # We can reach the neighbor node from the current node
                # through a path with a wider bottleneck that the one
                # we (perhaps) already knew.
                G.node[neighbor]['max_weight'] = m
                G.node[neighbor]['predecessor'] = current

                # Update the location of the bottleneck.
                if weight < m:
                    G.node[neighbor]['bottleneck'] = \
                      G.node[current]['bottleneck']
                else:
                    G.node[neighbor]['bottleneck'] = (current, neighbor)


def compute_global_max_weight_path(G, start=meta_reactant, stop=meta_product):
    """Compute global maximum weight path.

    """
    if (start == stop):
        return []

    compute_max_weight(G, start)

    bottleneck = G.node[stop]['bottleneck']

    return compute_global_max_weight_path(G, start, bottleneck[0]) + \
      [(bottleneck[0], bottleneck[1])] + \
      compute_global_max_weight_path(G, stop, bottleneck[1])


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 4:
        print('Usage: %s matrix-file flux-file milestones-file'
              % sys.argv[0], file=sys.stderr)
        sys.exit(-1)

    matrix_file = sys.argv[1]
    flux_file = sys.argv[2]
    milestones_file = sys.argv[3]

    reactants, products = read_special_nodes(milestones_file)

    G = read_connectivity_graph(matrix_file, flux_file,
                                reactants, products)

    print('Computing global maximum weight path...', file=sys.stderr)
    path = compute_global_max_weight_path(G)
    print('Done.', file=sys.stderr)
    for p in path[1:]:
        print(p[0])
