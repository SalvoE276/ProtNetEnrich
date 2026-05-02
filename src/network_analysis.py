import networkx as nx
from pyvis.network import Network
import numpy as np
import sys

def compute_candidate_proteins(data_dict: dict, query_proteins: list):
    data = data_dict.copy()
    for e in query_proteins:
        data.pop(e)
    return [node for node, v in data.items() if v > np.mean(list(data.values()))+(2*np.std(list(data.values())))] # values more than 2 std

if __name__=='__main__':
    ### Collect cli args ###
    try:
        network_path = sys.argv[1]
        output_folder = sys.argv[2]
    except IndexError:
        raise FileNotFoundError("PPI file data is not valid. Only JSON files are accepted.")
    
    ### Load PPI network ###
    net = nx.read_gml(network_path)


    ### Compute clustering coefficients and centralities ###
    clustering_coeffs = nx.clustering(net, net.nodes)
    degree_centralities = nx.degree_centrality(net)
    closeness_centralities = nx.closeness_centrality(net)
    betweenness_centralities = nx.betweenness_centrality(net)

    query_proteins = [p[0] for p in net.nodes(data='query_protein') if p[1]] # exctract original query proteins
    topology_metrics = {
        'ccoef' : compute_candidate_proteins(clustering_coeffs, query_proteins),
        'dc' : compute_candidate_proteins(degree_centralities, query_proteins),
        'cc' : compute_candidate_proteins(closeness_centralities, query_proteins),
        'bc' : compute_candidate_proteins(betweenness_centralities, query_proteins)
    }

    # Eval optional cli parameters
    if '--all_metrics' not in sys.argv:
        hubs = []
        for metric in topology_metrics.keys():
            if f'-{metric}' in sys.argv:
                hubs.extend(topology_metrics[metric])
        hubs = list(set(hubs))
    else:
        hubs = [protein for v in topology_metrics.values() for protein in v]
        hubs = list(set(hubs))
    
    # Change color to candidate hubs
    nx.set_node_attributes(net, {p:"red" for p in hubs}, name='color')


    ### Visualization network ###
    inet = Network(notebook=True, cdn_resources="in_line", select_menu=True)
    nx.set_edge_attributes(net, {e:"#869BC4" for e in list(net.edges)}, name='color') # fix inconsistent edges color
    inet.from_nx(net)
    inet.show(output_folder+"/interactive_network.html")

    ### Save Network ###
    nx.write_gml(net, output_folder+"/network.gml")