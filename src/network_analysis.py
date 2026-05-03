import networkx as nx
from pyvis.network import Network
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def compute_candidate_proteins(data_dict: dict, query_proteins: list, n_stds: int = 2):
    data = data_dict.copy()
    for e in query_proteins:
        data.pop(e)
    return [node for node, v in data.items() if v > np.mean(list(data.values()))+(n_stds*np.std(list(data.values())))] # values more than 2 std

def compute_and_save_diagnostics(net: nx.classes.graph, ccoef: list, dc: list, cc: list, bc: list, n_stds: int = 2):
    global output_folder

    metrics_dict = {
        "clustering_coefficients" : np.array(list(ccoef.values())),
        "degree_centralities" : np.array(list(dc.values())),
        "closeness_centralities" : np.array(list(cc.values())),
        "betweenness_centralities" : np.array(list(bc.values()))
    }

    ### Compute degree distribution ###
    degrees = np.array([d for n, d in net.degree()])
    plt.figure(figsize=(20, 11))
    sns.histplot(degrees, color="#4C72B0", alpha=0.6, edgecolor='black', bins=50)

    plt.title("Degrees distribution", fontsize=20, fontweight='bold')
    plt.ylabel("Count")
    plt.xlabel("Degree")
    plt.savefig(output_folder+"/degrees_distribution.png", dpi=150, bbox_inches="tight")
    plt.clf()


    ### Plot topology metrics distribution ###
    fig, axes = plt.subplots(2, 2, figsize=(20, 11))
    axes = axes.flatten() # flatten vector for axes

    # Plot colors
    colors = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]

    # Iteration over each metric
    for ax, (col, values), color in zip(axes, metrics_dict.items(), colors):
        sns.histplot(values, ax=ax, color=color, alpha=0.6, edgecolor="black", bins=50)

        mean, std = values.mean(), values.std() # computed over np.array

        # add mean and positive std to check the metric relevance in computing candiate hubs
        ax.axvline(mean, color="black", linestyle='--', linewidth=1.8, label=f"Mean = {mean:.2f}")
        ax.axvline(mean+(n_stds*std), color="gray", linewidth=1.2, linestyle="--", label=f"Std  = {std:.2f}")

        ax.set_title(col, fontsize=13, fontweight="bold")
        ax.set_ylabel("Count")
        ax.legend(fontsize=9)

    plt.suptitle("Topology metrics distributions", fontsize=15, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_folder+"/topology_metrics_plots.png", dpi=150, bbox_inches="tight")
    plt.clf()


if __name__=='__main__':
    ### Collect cli args ###
    try:
        network_path = sys.argv[1]
        output_folder = sys.argv[2]
    except IndexError:
        raise FileNotFoundError("PPI file data is not valid. Only JSON files are accepted.")
    
    # Optional params handling
    params = {
        "-stds" : "n_stds",
    }
    params_values = {params[p]:float(sys.argv[sys.argv.index(p)+1]) for p in params.keys() if p in sys.argv}
    
    ### Load PPI network ###
    net = nx.read_gml(network_path)


    ### Compute clustering coefficients and centralities ###
    clustering_coeffs = nx.clustering(net, net.nodes)
    degree_centralities = nx.degree_centrality(net)
    closeness_centralities = nx.closeness_centrality(net)
    betweenness_centralities = nx.betweenness_centrality(net)

    # Compute diagnostic plots
    compute_and_save_diagnostics(net, clustering_coeffs, degree_centralities, closeness_centralities, betweenness_centralities, **params_values)


    query_proteins = [p[0] for p in net.nodes(data='query_protein') if p[1]] # exctract original query proteins
    hubs_per_topology_metrics = {
        'ccoef' : compute_candidate_proteins(clustering_coeffs, query_proteins, **params_values),
        'dc' : compute_candidate_proteins(degree_centralities, query_proteins, **params_values),
        'cc' : compute_candidate_proteins(closeness_centralities, query_proteins, **params_values),
        'bc' : compute_candidate_proteins(betweenness_centralities, query_proteins, **params_values)
    }


    # Eval optional cli parameters
    if '--all_metrics' not in sys.argv:
        hubs = []
        for metric in hubs_per_topology_metrics.keys():
            if f'-{metric}' in sys.argv:
                hubs.extend(hubs_per_topology_metrics[metric])
        hubs = list(set(hubs))
    else:
        hubs = [protein for v in hubs_per_topology_metrics.values() for protein in v]
        hubs = list(set(hubs))
    

    # Change color to candidate hubs
    nx.set_node_attributes(net, {p:"red" for p in hubs}, name='color')

    # Add candidate_hub flag to all nodes
    nx.set_node_attributes(net, {p:True for p in hubs}, name='candidate_hub')
    nx.set_node_attributes(net, {n:False for n in list(net.nodes) if n not in hubs}, name='candidate_hub')

    ### Visualization network ###
    inet = Network(notebook=True, cdn_resources="in_line", select_menu=True)
    nx.set_edge_attributes(net, {e:"#869BC4" for e in list(net.edges)}, name='color') # fix inconsistent edges color
    nx.set_edge_attributes(net, nx.get_edge_attributes(net, 'score'), name='value') # add edges width scaling based on STRING score
    inet.from_nx(net)
    inet.show(output_folder+"/interactive_network.html")

    ### Save Network ###
    nx.write_gml(net, output_folder+"/network.gml")