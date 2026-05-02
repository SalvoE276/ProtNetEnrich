import networkx as nx
import sys


if __name__=='__main__':
    ### Collect cli args ###
    try:
        network_path = sys.argv[1]
        #output_folder = sys.argv[2]
    except IndexError:
        raise FileNotFoundError("PPI file data is not valid. Only JSON files are accepted.")
    
    ### Load PPI network ###
    net = nx.read_gml(network_path)

    print(net.nodes(data='id'))
    print(net.edges(data=True))