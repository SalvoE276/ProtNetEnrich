import sys, json
import requests as r


def retrieve_STRING_ppi_data(proteins: list, interactors_limit: int = 30, taxon: int = 9606):
    """Retrieves protein-protein interaction (PPI) data from the STRING-DB API.

    This function queries the STRING-DB web service to find physical interaction 
    partners for a given list of protein identifiers. It handles the formatting 
    of the API request and manages potential API limitations regarding query size.

    Args:
        proteins (list): A list of protein (only HGNC symbols accepted) to query.
        interactors_limit (int, optional): The maximum number of interaction partners 
            to retrieve per protein. Defaults to 30.
        taxon (int, optional): The NCBI taxonomy ID representing the species to be 
            queried. Defaults to 9606 (Homo sapiens).

    Returns:
        list[dict]: A list of dictionaries containing the retrieved interaction partner 
            data in JSON format, as provided by the STRING-DB API.

    Raises:
        ValueError: If the `proteins` list is empty.
        ConnectionError: If the API returns a 400 error (invalid query) or 
            if there is a failure in communicating with the STRING-DB server.

    Note:
        The function enforces a limit of 600 proteins per request to comply 
        with STRING-DB's API constraints. If the input list exceeds this, 
        the function will truncate the list and print a warning.
    """
    max_proteins = 600 # Max single request GET STRING-DB.org
    if len(proteins) == 1:
        query = proteins[0]
    elif len(proteins) == 0:
        raise ValueError('Protein list is empty')
    else:
        if len(proteins) > max_proteins:
            query = "%0d".join(proteins[:700])
            print(f"Query protein max length exceeded. Computing only the first {max_proteins} proteins")
        else:
            query = "%0d".join(proteins)
    url = f'https://string-db.org/api/json/interaction_partners?identifiers={query}&limit={interactors_limit}&species={taxon}&network_type=physical'
    response = r.get(url)
    if not response.ok:
        if response.status_code == 400:
            raise ConnectionError("Query data not valid")
        else:
            raise ConnectionError(f'Error during connection to STRING-db.org. Status code: {response.status_code}')
    return response.json()

if __name__=='__main__':
    ### Collect cli args ###
    try:
        source_file = sys.argv[1]
        output_path = sys.argv[2]
    except IndexError:
        raise FileNotFoundError('Protein source/output file not valid')
    
    # Optional params handling
    params = {
        "-il" : "interactors_limit",
        "-tx" : "taxon"
    }
    params_values = {params[p]:sys.argv[sys.argv.index(p)+1] for p in params.keys() if p in sys.argv}

    
    ### Retrieve and save physical protein interaction data ###
    with open(source_file, 'r') as file:
        proteins = file.readlines()
    
    data = retrieve_STRING_ppi_data(proteins, **params_values)
    with open(output_path+"/ppi_data.json", 'w') as jsonfile:
        json.dump(data, jsonfile)