import gseapy as gp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json, sys, os


if __name__=='__main__':
    ### Collect cli args ###
    try:
        borda_ranking_path = sys.argv[1]
        output_folder = sys.argv[2]
    except IndexError:
        raise FileNotFoundError("Wrong use of command line arguments.")
    
    with open(borda_ranking_path, 'r') as jsonfile:
        ranked_hubs = json.load(jsonfile)
        ranked_hubs = pd.DataFrame(ranked_hubs)
    

    # Run pre-ranked GSEA
    gsea = gp.prerank(
        rnk=ranked_hubs,               # Your ranked Series or DataFrame
        gene_sets='KEGG_2021_Human',    # Gene set database (see options below)
        outdir=output_folder+'/GSEA/',          # Output directory
        permutation_num=1000,           # Number of permutations (≥1000 for publication)
        min_size=5,                    # Min genes in a gene set to test
        max_size=500,                   # Max genes in a gene set to test
        threads=20,
        seed=478536,
        verbose=False
    )

    # Access results
    results_df = gsea.res2d
    results_df = results_df.sort_values(by='NES', ascending=False)
    terms = list(results_df['Term'])
    print(results_df)
    axs = gsea.plot(terms=terms[0])
    axs.set_size_inches(20, 11)
    axs.savefig('test.png')