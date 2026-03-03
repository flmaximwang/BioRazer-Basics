import requests
from .info import AFDBEntry


def uniprot_to_entries(uniprot_id: str):
    r = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}")
    entries = r.json()
    for i in range(len(entries)):
        entries[i] = AFDBEntry(data=entries[i])
    return entries
