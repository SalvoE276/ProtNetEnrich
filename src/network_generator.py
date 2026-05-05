import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import sys


if __name__=="__main__":
    ### Collect cli args ###
    try:
        jsonfile_path = sys.argv[1]
        output_folder = sys.argv[2]
        original_source_file = sys.argv[3]
    except IndexError:
        raise FileNotFoundError("PPI file data is not valid. Only JSON files are accepted.")

    ### Load data from json file ###
    df = pd.read_json(jsonfile_path)
    taxon = df["ncbiTaxonId"].unique()[0]

    df_input = df[['preferredName_A', 'stringId_A']].rename(columns={'preferredName_A': 'gene_symbol', 'stringId_A': 'ENS_id'})
    df_interactors = df[['preferredName_B', 'stringId_B']].rename(columns={'preferredName_B': 'gene_symbol', 'stringId_B': 'ENS_id'})

    all_proteins_df = pd.concat([df_input, df_interactors]).drop_duplicates().sort_values('gene_symbol').reset_index(drop=True)
    all_proteins_df['ENS_id'] = all_proteins_df['ENS_id'].str.replace(f"{taxon}.", "", regex=False)

    ### Build network ###
    net = nx.Graph()
    net.add_nodes_from((row.gene_symbol, {'ENS_id': row.ENS_id}) for row in all_proteins_df.itertuples())

    scores = [c for c in list(df.columns) if "scor" in c] # Select revelant metadata
    net.add_edges_from((row.preferredName_A, row.preferredName_B, {k:v for k,v in row._asdict().items() if k not in ['preferredName_A', 'preferredName_B', 'Index']}) 
                       for row in df[scores + ['preferredName_A', 'preferredName_B']].itertuples())

    ### Add signature and color for query proteins ###
    with open(original_source_file,"r") as textfile:
        query_proteins = textfile.readlines()
        query_proteins = [p.strip() for p in query_proteins]
    nx.set_node_attributes(net, {k:True for k in query_proteins}, name="query_protein")
    nx.set_node_attributes(net, {k:'green' for k in query_proteins}, name="color")
    nx.set_node_attributes(net, {k:False for k in list(net.nodes) if k not in query_proteins}, name="query_protein")
    nx.set_node_attributes(net, {k:'#b4cffb' for k in list(net.nodes) if k not in query_proteins}, name="color")
    
    ### Save Network ###
    nx.write_gml(net, output_folder+"/raw_network.gml")