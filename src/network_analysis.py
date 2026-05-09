import networkx as nx
from pyvis.network import Network
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys, json, os, datetime, random


def compute_candidate_proteins(data_dict: dict, query_proteins: list, n_stds: int = 2):
    """
    Identifies candidate hub proteins.

    Calculates the mean and standard deviation of metrics for all proteins 
    excluding the provided query proteins. A protein is considered a 
    candidate if its metric value exceeds the threshold: 
    mean + (n_stds * standard_deviation).

    Args:
        data_dict (dict): A dictionary where keys are protein names (str) 
            and values are numerical topology metrics (float).
        query_proteins (list): A list of protein names (str) to be excluded 
            from the statistical baseline calculation.
        n_stds (int, optional): The number of standard deviations above the 
            mean to use as the threshold. Defaults to 2.

    Returns:
        list: A list of protein names (strings) that meet the outlier criteria, 
            sorted by their metric value in descending order.
    """
    data = data_dict.copy()
    for e in query_proteins:
        data.pop(e)
    data_mean = np.mean(list(data.values()))
    data_std = np.std(list(data.values()))
    candidates = [{'gene':node, 'value':v} for node, v in data.items() if v > data_mean+(n_stds*data_std)] # values more than 2 std
    return [e['gene'] for e in sorted(candidates, key=lambda x: x['value'], reverse=True)]

def borda_count_ranking(list_of_ranks: list[list]):
    """
    Performs Borda count aggregation on multiple ranked lists.

    Args:
        list_of_ranks (list[list]): A list containing multiple ranked lists, 
            where each sub-list represents the order of elements from 
            highest to lowest preference/score.

    Returns:
        list: A list of tuples containing (element, aggregated_score), 
            sorted by the aggregated score.
    """
    random.seed(284957) # for consistency

    # Set scrambled (unbiased) list of gene
    scrambled_reference_list = list_of_ranks[0].copy()
    random.shuffle(scrambled_reference_list)

    if not list_of_ranks:
        return []
    scores = {g:0 for g in scrambled_reference_list}
    n = len(list_of_ranks[0])
    for rank in list_of_ranks:
        for i, gene in enumerate(rank):
            points = (n - 1) - i
            scores[gene] = scores[gene] + points
            
    # Sort genes by score in descending order
    ranking = sorted(scores.items(), key=lambda item: item[1], reverse=True)
    return ranking
    

def compute_and_save_diagnostics(net: nx.classes.graph, ccoef: list, dc: list, cc: list, bc: list, n_stds: int = 2):
    """
    Generates and saves diagnostic visualization plots for network topology and edge scores.

    This function creates a 'diagnostic_plots' subdirectory within the global 
    `output_folder` and generates three specific types of plots:
    1. Distribution of STRING scores (based on edge attributes).
    2. Degree distribution of the network nodes.
    3. Distribution of topological metrics (clustering, degree, closeness, 
       and betweenness centralities) with mean and threshold lines.

    Args:
        net (nx.Graph): A NetworkX graph object containing edge data (e.g., STRING scores).
        ccoef (dict): A dictionary mapping node identifiers to their clustering 
            coefficients.
        dc (dict): A dictionary mapping node identifiers to their degree 
            centralities.
        cc (dict): A dictionary mapping node identifiers to their closeness 
            centralities.
        bc (dict): A dictionary mapping node identifiers to their betweenness 
            centralities.
        n_stds (int, optional): The selected number of standard deviations
            to compute the thresholds. Defaults to 2.

    Returns:
        None

    Note:
        - Creates a directory named 'diagnostic_plots' inside `output_folder` 
          if it does not already exist.
        - Saves 'scores_distribution.png' to the diagnostic folder.
        - Saves 'degrees_distribution.png' to the diagnostic folder.
        - Saves 'topology_metrics_plots.png' to the diagnostic folder.
        - Relies on the existence of a global variable `output_folder` 
          defining the base path.
    """
    global output_folder

    # create plot output foder
    try:
        os.mkdir(output_folder+"/diagnostic_plots")
    except:
        print("diagnostic_plots folder already exists")

    metrics_dict = {
        "clustering_coefficients" : np.array(list(ccoef.values())),
        "degree_centralities" : np.array(list(dc.values())),
        "closeness_centralities" : np.array(list(cc.values())),
        "betweenness_centralities" : np.array(list(bc.values()))
    }

    ### Plot STRING scores data ###
    df_scores = pd.DataFrame([e[2] for e in list(net.edges(data=True))])
    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(20, 11))
    axes = axes.flatten()

    for ax, col in zip(axes, df_scores.columns):
        col_mean = df_scores[col].mean()
        sns.histplot(data=df_scores, x=col, ax=ax, color="#4C72B0", alpha=0.6, edgecolor='black', bins=50)
        ax.axvline(col_mean, color="black", linestyle='--', linewidth=1.8, label=f"Mean = {col_mean:.2f}") # add mean line
        ax.set_title(col)
        ax.set_xlabel("")

    plt.suptitle("Distribution of STRING scores", fontsize=20, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_folder+"/diagnostic_plots/scores_distribution.png", dpi=150, bbox_inches="tight")
    plt.clf()


    ### Compute degree distribution ###
    degrees = np.array([d for n, d in net.degree()])
    plt.figure(figsize=(20, 11))
    sns.histplot(degrees, color="#4C72B0", alpha=0.6, edgecolor='black', bins=50)

    plt.title("Degrees distribution", fontsize=20, fontweight='bold')
    plt.ylabel("Count")
    plt.xlabel("Degree")
    plt.savefig(output_folder+"/diagnostic_plots/degrees_distribution.png", dpi=150, bbox_inches="tight")
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
    plt.savefig(output_folder+"/diagnostic_plots/topology_metrics_plots.png", dpi=150, bbox_inches="tight")
    plt.clf()


if __name__=='__main__':
    ### Collect cli args ###
    try:
        network_path = sys.argv[1]
        output_folder = sys.argv[2]
    except IndexError:
        raise FileNotFoundError("Wrong use of command line arguments.")
    
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
    topological_measures = {
        'ccoef' : clustering_coeffs,
        'dc' : degree_centralities,
        'cc' : closeness_centralities,
        'bc' : betweenness_centralities
    }

    # Compute diagnostic plots
    compute_and_save_diagnostics(net, clustering_coeffs, degree_centralities, closeness_centralities, betweenness_centralities, **params_values)


    query_proteins = [p[0] for p in net.nodes(data='query_protein') if p[1]] # exctract original query proteins
    hubs_per_topology_metrics = {k:compute_candidate_proteins(v, query_proteins, **params_values) for k,v in topological_measures.items()}

    # Eval optional cli metric parameters
    used_topology_metrics = []
    metrics_cands = {}
    if '--all_metrics' not in sys.argv:
        for metric in hubs_per_topology_metrics.keys():
            if f'-{metric}' in sys.argv:
                metrics_cands[metric] = hubs_per_topology_metrics[metric]
                used_topology_metrics.append(metric)
    else:
        metrics_cands = hubs_per_topology_metrics
        used_topology_metrics.extend(list(hubs_per_topology_metrics.keys()))

    # Get all candidate hubs
    cand_hubs = sorted(set([e for k,v in metrics_cands.items() for e in v]))

    # Get hubs topological measures
    hubs_metrics = {hub:{metric:dict(topological_measures[metric])[hub] for metric in used_topology_metrics} for hub in cand_hubs}

    # Eval full list of cand hubs for each used metric
    full_rank = []
    for metric in used_topology_metrics:
        full_rank.append([gene for gene, vals in dict(sorted(hubs_metrics.items(), key=lambda item: item[1][metric], reverse=True)).items()])

    # Aggregation of ranks
    ranking = borda_count_ranking(full_rank)
    hubs = [e[0] for e in ranking] # extract only ranked_gene symbols

    # Save Borda ranking data
    with open(output_folder+"/borda_ranking.json", 'w') as jsonfile:
        json.dump(ranking, jsonfile)



    # Change color to candidate hubs
    nx.set_node_attributes(net, {p:"red" for p in hubs}, name='color')

    # Add candidate_hub flag to all nodes
    nx.set_node_attributes(net, {p:True for p in hubs}, name='candidate_hub')
    nx.set_node_attributes(net, {n:False for n in list(net.nodes) if n not in hubs}, name='candidate_hub')

    ### Visualization network ###
    inet = Network(height='90vh', width='100%', notebook=True, cdn_resources="in_line", select_menu=True) # do not use % in height argument
    nx.set_edge_attributes(net, {e:"#869BC4" for e in list(net.edges)}, name='color') # fix inconsistent edges color
    nx.set_edge_attributes(net, nx.get_edge_attributes(net, 'score'), name='value') # add edges width scaling based on STRING score

    # Remove nodeds with degree == 1 for to remove visualization (if requested)
    if '--remove_dead_ends' in sys.argv:
        simplified_net = net.copy()
        simplified_net.remove_nodes_from([n for n, d in dict(simplified_net.degree()).items() if d == 1])
        inet.from_nx(simplified_net)
    else:
        inet.from_nx(net)
    inet.show(output_folder+"/interactive_network.html")

    ### Save Network ###
    nx.write_gml(net, output_folder+"/network.gml")

    ### Save candidate hub proteins ###
    output_data = {'header':{'datetime':str(datetime.datetime.now()), 'topology_metrics':used_topology_metrics, 'std_threshold':params_values['n_stds'] if 'n_stds' in params_values.keys() else 2.0}, "Ranked":True,
                   'hubs':[{'gene_symbol':protein, 'id':net.nodes[protein]['ENS_id']} for protein in hubs]}

    with open(output_folder+"/hubs.cand.json", "w") as jsonfile:
        json.dump(output_data, jsonfile)