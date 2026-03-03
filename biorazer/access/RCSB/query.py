from rcsbapi.search import AttributeQuery, NestedAttributeQuery


def query_uniprot(uniprot_id, return_type="entry"):
    """
    This function queries all RCSB codes that contain the target protein denoted by the Uniprot accession ID
    Args:
    - return_type: str, one of ["entry", "polymer_entity", "nonpolymer_entity", "branched_entity"]
    """
    assert return_type in [
        "entry",
        "polymer_entity",
    ]
    q_uniprot = NestedAttributeQuery(
        AttributeQuery(
            attribute="rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
            operator="exact_match",
            value=uniprot_id,
        ),
        AttributeQuery(
            attribute="rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
            operator="exact_match",
            value="UniProt",
        ),
    )
    res = list(q_uniprot(return_type))
    return res
