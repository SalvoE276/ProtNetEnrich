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
    """
    Loads Ensembl IDs from a hubs.cands.json file.

    Args:
        filepath (str): The file path to the JSON file containing hub data.

    Returns:
        list[str]: A list of Ensembl IDs extracted from the 'hubs' field in the JSON.
    """
    with open(filepath, 'r') as hubfile:
        cand_hubs_file = json.load(hubfile)
    return [e['id'] for e in cand_hubs_file['hubs']]

def from_ENSP_to_ENTREZ_id(ids: list, organism: str = 'hsapiens'):
    """
    Converts a list of Ensembl IDs to Entrez Gene IDs using gProfiler.

    This function queries the gProfiler API for the specified organism. 
    If any IDs fail to convert (resulting in 'None'), it prints a warning 
    and filters them out of the returned list.

    Args:
        ids (list): A list of Ensembl identifiers to be converted.
        organism (str) : The selected organism for the conversion. 
            Default to 'hsapiens'.

    Returns:
        list: A list of converted Entrez Gene IDs.
    """
    gp = GProfiler(return_dataframe=True)
    gp_df = gp.convert(organism=organism, query=ids, target_namespace='ENTREZGENE_ACC')
    result = list(gp_df['converted'])
    if 'None' in result:
        wrong_conversion = [g for g, c in zip(ids, result) if c == 'None']
        print("Some genes id where not converted: ")
        for g in wrong_conversion:
            print(g)
        print("Skipping these genes in enrichment...")
        return [id for id in result if id != 'None']
    else:
        return result




def run_enrichment(
    entrez_ids: list[int],
    categories: list[str],
    p_value_cutoff: float = 0.05,
    min_genes: int = 1,
    max_genes: int = 1500,
    max_results: int = 50,
    correction: str = "FDR"):
    """
    Performs gene enrichment analysis using the ToppGene API.

    Args:
        entrez_ids (list[int]): A list of Entrez Gene IDs to analyze.
        categories (list[str]): A list of ToppGene category strings to query.
        p_value_cutoff (float): The p-value threshold for significance. Defaults to 0.05.
        min_genes (int): The minimum number of genes required in a term. Defaults to 1.
        max_genes (int): The maximum number of genes to consider in a term. Defaults to 1500.
        max_results (int): The maximum number of enrichment results to return per category. Defaults to 50.
        correction (str): The multiple testing correction method (e.g., 'FDR', 'Bonferroni'). Defaults to "FDR".

    Returns:
        dict[str, pd.DataFrame]: A dictionary where keys are the category names and 
            values are Pandas DataFrames containing the parsed enrichment results, 
            sorted by the FDR-corrected p-value.

    Raises:
        requests.exceptions.HTTPError: If the API request fails.
    """

    ENRICH_URL = "https://toppgene.cchmc.org/API/enrich"

    # Generation of POST request
    payload = {
        "Genes": [int(g) for g in entrez_ids if g != 'None'],
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
    """Prints a formatted summary of the enrichment results to the console.

    Iterates through the provided dictionary of results, displaying the top 
    rows for each category. It specifically looks for columns related to 
    term names, p-values, and gene lists.

    Args:
        results (dict[str, pd.Dataframe]): A dictionary where keys are 
            category names and values are pandas DataFrames containing 
            the enrichment data.

    Returns:
        None
    """
    for cat in results.keys():
        df = results[cat]
        if df is None or df.empty:
            continue
        print(f"\n\n{'-'*50} {cat} - {len(df)} terms {'-'*50}\n")
        relevant_cols = ["Name", "PValue", "QValueFDRBH", "GenesInTermInQuery", "Genes"]
        print(df[relevant_cols].head(top_n).to_string(index=False))


# Accepts only dictionary of pd.Dataframe [see run_enrichment()]
def save_results(results: dict[str, pd.DataFrame], output_folder: str):
    """Saves the enrichment results to the local file system.

    Creates an 'enrichment_analysis' subdirectory within the specified 
    output folder. The function saves a single consolidated Excel file 
    containing all categories as separate sheets, and a collection of 
    individual CSV files (one for each category).

    Args:
        results (dict[str, pd.Dataframe]): A dictionary where keys are 
            category names and values are pandas DataFrames containing 
            the enrichment data.
        output_folder (str): The base directory path where the results 
            folder will be created and files will be stored.

    Raises:
        OSError: If the directory cannot be created or files cannot be written.

    Note:
        The directory structure created will be:
        output_folder/
        └── enrichment_analysis/
            ├── [category_name].csv
            └── enrichment_analysis.xlsx
    """
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