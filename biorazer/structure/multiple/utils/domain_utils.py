import numpy as np
import biotite.structure as struc

def generate_domain_selector(struc: struc.AtomArray, domain_identifier: dict):
    res_selector = np.ones(struc.array_length(), dtype=bool)
    for kkey in domain_identifier:
        if isinstance(domain_identifier[kkey], list):
            step_selector = list(map(lambda x: x in domain_identifier[kkey], struc.get_annotation(kkey)))
        elif iter(domain_identifier[kkey]):
            identifier_list = list(domain_identifier[kkey])
            step_selector = list(map(lambda x: x in identifier_list, struc.get_annotation(kkey)))
        else:
            step_selector = struc.get_annotation(kkey) == domain_identifier[kkey]
        res_selector = np.logical_and(res_selector, step_selector)
    return res_selector

def select_domain(struc: struc.AtomArray, domain_identifier: dict):
    res_selector = generate_domain_selector(struc, domain_identifier)
    return struc[res_selector]

def assign_domains(struc: struc.AtomArray, domains):
    '''
    domain: dict of domain_identifiers
    '''
    struc.del_annotation("domain")
    struc.add_annotation("domain", 'U3')
    for domain in domains:
        domain_identifier = domains[domain]
        domain_selector = generate_domain_selector(struc, domain_identifier)
        for i, flag in enumerate(domain_selector):
            if flag:
                struc.domain[i] = domain

def get_indices_for_domains(struc: struc.AtomArray, domains):
    final_selector = np.ones(struc.array_length(), dtype=bool)
    domain_indices = {}
    for domain in domains:
        domain_identifier = domains[domain]
        domain_selector = generate_domain_selector(struc, domain_identifier)
        final_selector = np.logical_and(final_selector, domain_selector)
    return np.where(final_selector), domain_indices

def get_binder_tdomain_complex_indices(struc: struc.AtomArray, domain_identifier):
    tdomain_selector = generate_domain_selector(struc, DOMAINS_IN_PREDICTED[domain_identifier])
    binder_selector = generate_domain_selector(struc, DOMAINS_IN_PREDICTED["Bdr"])
    complex_selector = np.logical_and(tdomain_selector, binder_selector)
    return np.where(complex_selector)