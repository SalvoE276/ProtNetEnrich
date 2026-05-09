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
        selected_gene_set = sys.argv[3]
    except IndexError:
        raise ValueError("Wrong use of command line arguments.")
    
    with open(borda_ranking_path, 'r') as jsonfile:
        ranked_hubs = json.load(jsonfile)
        ranked_hubs = pd.DataFrame(ranked_hubs)

    ### Default gene sets ###
    default_gene_sets = {
        'GO_26' : ['GO_Biological_Process_2026', 'GO_Cellular_Component_2026', 'GO_Molecular_Function_2026'],
        'KEGG_26' : 'KEGG_2026',
        'Disease' : ['OMIM_Expanded', 'Orphanet_Augmented_2021']
    }
    if selected_gene_set not in default_gene_sets.keys() and selected_gene_set == '--custom_set':
        try:
            custom_gene_set = sys.argv[4]
        except IndexError:
            raise IndexError("Missing custom gene set specification")
        gene_set = custom_gene_set
    elif selected_gene_set not in default_gene_sets and selected_gene_set != '--custom_set':
        raise ValueError("Wrong use of command line arguments.")
    elif selected_gene_set in default_gene_sets.keys():
        gene_set = default_gene_sets[selected_gene_set]
    

    # Run pre-ranked GSEA
    gsea = gp.prerank(
        rnk=ranked_hubs,
        gene_sets=gene_set,
        outdir=output_folder+'/GSEA/',
        permutation_num=1000,
        min_size=5,     # Min genes in tested gene set
        max_size=500,   # Max genes in tested gene set
        threads=os.cpu_count(),
        seed=478536,
        verbose=False,
        no_plot=True
    )
 

    # Generation of GSEA report
    results_df = gsea.res2d.sort_values(by=['FDR q-val', 'NES'], ascending=[True, False])
    terms = results_df['Term']

    results_df = results_df[results_df['FDR q-val'] < 0.25]
    results_df = results_df.drop(columns=['Name'])
    results_df.to_csv(output_folder+'/GSEA/GSEA_report.csv', index=False)


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