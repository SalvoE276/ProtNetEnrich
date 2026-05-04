import requests as r
from gprofiler import GProfiler
import pandas as pd
import json, sys, os


# All valid category types
ALL_CATEGORIES = [
    "GeneOntologyMolecularFunction",
    "GeneOntologyBiologicalProcess",
    "GeneOntologyCellularComponent",
    "HumanPheno",
    "MousePheno",
    "Domain",
    "Pathway",
    "Pubmed",
    "Interaction",
    "Cytoband",
    "TFBS",
    "GeneFamily",
    "Coexpression",
    "CoexpressionAtlas",
    "ToppCell",
    "Computational",
    "MicroRNA",
    "Drug",
    "Disease",
]


def load_ENS_ids(filepath: str):
    with open(filepath, 'r') as hubfile:
        cand_hubs_file = json.load(hubfile)
    return [e['id'] for e in cand_hubs_file['hubs']]

def from_ENSP_to_ENTREZ_id(ids: list):
    gp = GProfiler(return_dataframe=True)
    gp_df = gp.convert(organism='hsapiens', query=ids, target_namespace='ENTREZGENE_ACC')
    return list(gp_df['converted'])




def run_enrichment(
    entrez_ids: list[int],
    categories: list[str],
    p_value_cutoff: float = 0.05,
    min_genes: int = 1,
    max_genes: int = 1500,
    max_results: int = 50,
    correction: str = "FDR"):

    ENRICH_URL = "https://toppgene.cchmc.org/API/enrich"

    # Generation of POST request
    payload = {
        "Genes": [int(g) for g in entrez_ids],
        "Categories": [
            {
                "Type":       cat,
                "PValue":     p_value_cutoff,
                "MinGenes":   min_genes,
                "MaxGenes":   max_genes,
                "MaxResults": max_results,
                "Correction": correction,
            }
            for cat in categories
        ],
    }

    # POST request
    response = r.post(
        ENRICH_URL,
        json=payload,
        headers={"Content-Type": "application/json"},
        timeout=120,
    )

    # Check server response
    response.raise_for_status()

    ### Parse and sort response ###
    annotations = response.json()['Annotations']

    # Reorder data into raw = {category: list[dict]}
    raw = {}
    for a in annotations:
        cat = a.get("Category", "Unknown")
        raw.setdefault(cat, []).append({
            "ID":                 a.get("ID"),
            "Name":               a.get("Name"),
            "PValue":             a.get("PValue"),
            "QValueFDRBH":        a.get("QValueFDRBH"),
            "QValueFDRBY":        a.get("QValueFDRBY"),
            "QValueBonferroni":   a.get("QValueBonferroni"),
            "TotalGenes":         a.get("TotalGenes"),
            "GenesInTerm":        a.get("GenesInTerm"),
            "GenesInQuery":       a.get("GenesInQuery"),
            "GenesInTermInQuery": a.get("GenesInTermInQuery"),
            "Source":             a.get("Source"),
            "URL":                a.get("URL"),
            "Genes":              [g["Symbol"] for g in a.get("Genes", [])],
        })

    # Generation of Dataframe for each category
    dfs = {cat: pd.DataFrame(dictdata).sort_values("QValueFDRBH") for cat, dictdata in raw.items()}
    return dfs


# Accepts only dictionary of pd.Dataframe [see run_enrichment()]
def print_summary(results: dict[str, pd.DataFrame], top_n: int = 5):
    for cat in results.keys():
        df = results[cat]
        if df is None or df.empty:
            continue
        print(f"\n\n{'-'*50} {cat} - {len(df)} terms {'-'*50}\n")
        relevant_cols = ["Name", "PValue", "QValueFDRBH", "GenesInTermInQuery", "Genes"]
        print(df[relevant_cols].head(top_n).to_string(index=False))


# Accepts only dictionary of pd.Dataframe [see run_enrichment()]
def save_results(results: dict[str, pd.DataFrame], output_folder):
    ### Generation of output folder ###
    savepath = f"{output_folder}/enrichment_analysis"
    try:
        os.mkdir(savepath)
    except:
        print(f"{savepath} folder already exists")

    ### Save xlsx ###
    with pd.ExcelWriter(savepath+'/toppgene_results.xlsx') as writer:
        for cat, df in results.items():
            df.to_excel(writer, sheet_name=cat,  index=False)

    ### Save csvs ###
    csv_path = f'{savepath}/csvs'
    try:
        os.mkdir(csv_path)
    except:
        print(f'{csv_path} folder already exists')
    
    for cat, df in results.items():
        df.to_csv(f'{csv_path}/{cat}.csv', index=False, sep=',')



if __name__ == "__main__":
    ### Collect cli args ###
    try:
        candidate_hubs_path = sys.argv[1]
        output_folder = sys.argv[2]
    except IndexError:
        raise FileNotFoundError("Wrong use of command line arguments.")


    ### Loading and converion of candiate hub ids ###
    gene_list = load_ENS_ids(candidate_hubs_path)
    gene_list = from_ENSP_to_ENTREZ_id(gene_list)

    # Handle bad conversion
    if len(gene_list) == 0:
        raise SystemExit("Enrichment gene list is empty.")

    ### Toppgene enrichment ###
    results = run_enrichment(gene_list, categories=ALL_CATEGORIES)

    ### Print summary and save data ###
    print_summary(results, top_n=5) # for rapid inspection
    save_results(results, output_folder)