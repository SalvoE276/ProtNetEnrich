import sys
import requests as r


def retrieve_STRING_ppi_data(proteins: list, interactors_limit: int = 30):
    if len(proteins) == 1:
        query = proteins[0]
    elif len(proteins) == 0:
        raise ValueError('Protein list is empty')
    else:
        query = "%0d".join(proteins)
    url = f'https://string-db.org/api/json/interaction_partners?identifiers={query}&limit={interactors_limit}'
    response = r.get(url)
    if not response.ok:
        raise ConnectionError('Error during connection to STRING.org')
    return response.json()

if __name__=='__main__':
    source_file = sys.argv[1]
    with open(source_file, 'r') as file:
        proteins = file.readlines()
    data = retrieve_STRING_ppi_data(proteins, interactors_limit=1)
    for e in data:
        print(e)