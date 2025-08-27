import biotite.structure as struc

def read_macro(macro_str):
    '''
    A macro looks like this:
    /model/segi/chain/resn`resi/name
    model: model name
    segi: segment identifier
    chain: chain identifier
    resn: residue name
    resi: residue number
    name: atom name
    '''
    import re
    macro_pattern = re.compile(r'/(?P<model>\w*)/(?P<segi>\w*)/(?P<chain>\w)/(?P<resn>\w+)`(?P<resi>\d+)/(?P<name>\w+)')
    match = macro_pattern.match(macro_str)
    if match:
        return match.groupdict()
    else:
        raise ValueError("Invalid macro string: %s" % macro_str)

def to_macro(model, segi, chain, resn, resi, name):
    return f"/{model}/{segi}/{chain}/{resn}`{resi}/{name}"

def atom_to_macro(atom: struc.Atom):
    return to_macro(
        "",
        "",
        atom.chain_id,
        atom.res_name,
        atom.res_id,
        atom.atom_name
    )

def chain_priority(chain):
    return sum([ord(c) for c in chain])

def name_priority(name):
    return sum([ord(c) for c in name])

def atom_priority(atom_macro):
    atom_macro_dict = read_macro(atom_macro)
    return (
        chain_priority(atom_macro_dict['chain']),
        int(atom_macro_dict['resi']),
        name_priority(atom_macro_dict['name'])
    )

