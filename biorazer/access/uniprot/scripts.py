import pandas as pd

def list_response_1(r):
    '''
    This helper function lists the Accession ID, Organism, Protein names, and Functions of the proteins in the response.
    '''
    res = r.json()["results"]
    res_df = pd.DataFrame()
    for i, entry in enumerate(res):
        if "recommendedName" in entry["proteinDescription"]:
            res_df.loc[
                i,
                [
                    "Accessions",
                    "Organisms",
                    "Protein names",
                ]
            ] = [
                entry["primaryAccession"],
                entry["organism"]["scientificName"],
                entry["proteinDescription"]["recommendedName"]["fullName"]["value"],
            ]
        else:
            res_df.loc[
                i,
                [
                    "Accessions",
                    "Organisms",
                ]
            ] = [
                entry["primaryAccession"],
                entry["organism"]["scientificName"],
            ]
            res_df.loc[i, "Protein names"] = ""
            for name in entry["proteinDescription"]["submissionNames"]:
                c_names = res_df.loc[i, "Protein names"]
                if not c_names:
                    res_df.loc[i, "Protein names"] = name["fullName"]["value"]
                else:
                    res_df.loc[i, "Protein names"] = ";".join([c_names, name["fullName"]["value"]])
        res_df.loc[i, "Functions"] = ""
        if "comments" in entry:
            for comment in entry["comments"]:
                for text in comment["texts"]:
                    c_funcs = res_df.loc[i, "Functions"]
                    if not c_funcs:
                        res_df.loc[i, "Functions"] = text["value"]
                    else:
                        res_df.loc[i, "Functions"] = ";".join([c_funcs, text["value"]])
        else:
            continue
    return res_df