
import requests
from .response import response_to_df

REST_API = "https://rest.uniprot.org"

def query_protein_name(
    protein_name: str,
    fields = ["accession", "gene_names", "protein_name", "organism_name", "cc_function"],
    size=500
):
    base_url = f"{REST_API}/uniprotkb/search"
    params = {
        "query": f"protein_name:{protein_name}",
        "fields": fields,
        "size": size
    }
    headers = {"accept": "application/json"}
    r = requests.get(base_url, params=params, headers=headers)
    data = r.json()["results"]
    return data