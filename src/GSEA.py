import gseapy as gp
import json, sys


if __name__=='__main__':
    ### Collect cli args ###
    try:
        candidate_hubs_path = sys.argv[1]
        output_folder = sys.argv[2]
    except IndexError:
        raise FileNotFoundError("Wrong use of command line arguments.")
    
    with open(candidate_hubs_path, 'r') as jsonfile:
        cand_hubs = json.load(jsonfile)
        cand_hubs = [e['gene_symbol'] for e in cand_hubs['hubs']]
    
    #print(cand_hubs)