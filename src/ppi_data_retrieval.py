import sys, json
import requests as r


def retrieve_STRING_ppi_data(proteins: list, interactors_limit: int = 30, taxon: int = 9606):
    if len(proteins) == 1:
        query = proteins[0]
    elif len(proteins) == 0:
        raise ValueError('Protein list is empty')
    else:
        query = "%0d".join(proteins)
    url = f'https://string-db.org/api/json/interaction_partners?identifiers={query}&limit={interactors_limit}&species={taxon}&network_type=physical'
    response = r.get(url)
    if not response.ok:
        if response.status_code == 400:
            raise ConnectionError("Query data not valid")
        else:
            raise ConnectionError('Error during connection to STRING-db.org')
    return response.json()

if __name__=='__main__':
    ### Collect cli args ###
    try:
        source_file = sys.argv[1]
        output_file = sys.argv[2]
    except IndexError:
        raise FileNotFoundError('Protein source/output file not specified')
    
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
    with open(output_file, 'w') as jsonfile:
        json.dump(data, jsonfile)