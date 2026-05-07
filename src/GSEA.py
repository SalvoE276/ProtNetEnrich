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
        rnk=ranked_hubs,
        gene_sets='KEGG_2021_Human',
        outdir=output_folder+'/GSEA/',
        permutation_num=1000,
        min_size=5,     # Min genes in tested gene set
        max_size=500,   # Max genes in tested gene set
        threads=os.cpu_count(),
        seed=478536,
        verbose=False,
        no_plot=True
    )

    # Print and save results
    results_df = gsea.res2d.sort_values(by='NES', ascending=False)
    terms = results_df['Term']
    print("### GSEA dataframe results ###")
    print(results_df)

    if "--save_plots" in sys.argv:
        print("Saving GSEA plots...")
        try:
            os.mkdir(output_folder+"/GSEA/plots")
        except:
            pass

        for t in terms:
            fig = gsea.plot(terms=t)
            fig.set_size_inches(20, 11)
            fig.savefig(f'{output_folder}/GSEA/plots/{t}.png')