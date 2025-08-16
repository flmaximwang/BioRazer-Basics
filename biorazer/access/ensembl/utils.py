import ensembl_rest
import requests, sys, json
from ...genome_analyzer.utils import nt_utils

def parse_kwargs_to_ext(ext, **kwargs):
    if kwargs:
        ext += "?"
        for key, value in kwargs.items():
            ext += f"{key}={value}&"
        ext = ext[:-1]
    return ext

def get_ensembl_response(ext, headers, **kwargs):
    '''
    虽然你在 URL 上看到 headers 和 ext 中额外的参数似乎是一样的, 但他们应该分开传入
    '''
    ext = parse_kwargs_to_ext(ext, **kwargs)
    response = requests.get(
        "http://rest.ensembl.org" + ext,
        headers=headers
    )
    if not response.ok:
        response.raise_for_status()
    return response

def ensembl_lookup(id, headers = {"Content-Type": "application/json"}, **kwargs):
    """
    The function `ensembl_lookup` sends a request to the Ensembl API to lookup an ID and returns the
    response.
    
    :param id: The `id` parameter in the `ensembl_lookup` function is used to specify the identifier for
    which you want to perform a lookup in the Ensembl database.
    
    :param headers: The `headers` parameter in the `ensembl_lookup` function is a dictionary that
    contains the HTTP headers to be included in the request.
    
    :return: The function `ensembl_lookup` is returning the response from the `get_ensembl_response`
    function with the specified parameters, including the `id` provided, custom `headers`, and additional keyword arguments (**kwargs) passed to the function.
    """

    return get_ensembl_response(ext=f"/lookup/id/{id}", headers=headers, content_type="application/json", **kwargs)

def ensembl_sequence_id(id, headers = {"Content-Type": "application/json"}, **kwargs):
    '''
    :param headers:
        - Content-Type: application/json: 同时获取序列以及元信息\n
        - Content-Type: text/x-fasta: 只获取序列
    '''
    return get_ensembl_response(ext=f"/sequence/id/{id}", headers=headers, **kwargs)

def get_simplified_exons_from_exon_list(exons):
    simplified_exons = []
    for exon in exons:
        simplified_exons.append((exon["start"], exon["end"]))
    return simplified_exons

def get_mRNA_from_ensembl_transcript(transcript_id):
    '''
    Joining exons from a transcript sequence to get the mRNA sequence.
    
    :param transcript_seq: The transcript sequence from which the mRNA sequence 
    will be extracted.
    
    :param exons: The exons of the transcript sequence. Every exon marked with a tuple, 
    where the first element is the start position and the second element is the end position.
    
    :return: The mRNA sequence.
    '''
    r = ensembl_sequence_id(transcript_id, headers={"Content-Type": "application/json"})
    transcript_seq = r.json()["seq"]
    r = ensembl_lookup(transcript_id, headers={"Content-Type": "application/json"}, expand="1")
    transcript_info = r.json()
    transcript_start = transcript_info["start"]
    exons = transcript_info["Exon"]
    simplified_exons = []
    for exon in exons:
        simplified_exons.append((exon["start"] - transcript_start, exon["end"] - transcript_start))
    mRNA_seq = nt_utils.get_mRNA_from_transcript(transcript_seq, simplified_exons)
    return mRNA_seq

def ensembl_transcript_haplotypes_get(id, species, headers = {"Content-Type": "application/json"}, **kwargs):
    return get_ensembl_response(ext=f"/transcript_haplotypes/{species}/{id}", headers=headers, **kwargs)