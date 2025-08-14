def pocket_atom_specifiers(poc_file_path):
    '''
    Parse .poc file from CASTP web server.  Lines are PDB format with pocket number at end:

    ATOM    165  NH2 ARG A  45      -9.395 -35.779   0.675  1.00 32.68   2  POC
    ATOM    194  O   VAL B  27     -12.428  -9.493   5.152  1.00 13.76   1  POC
    ATOM    196  CG1 VAL B  27     -11.788  -7.554   2.677  1.00 14.96   1  POC
    This script is modified from https://rbvi.github.io/chimerax-recipes/castp/castp.html
    '''

    with open(poc_file_path, 'r') as f:
        lines = f.readlines()

    pockets: dict[int: list[dict]] = {}
    for line in lines:
        if line.startswith('ATOM  '):
            atom_name, chain_id, res_num = line[12:16].strip(), line[21:22], int(line[22:26])
            atom_dict = {'chain_id': chain_id, 'res_num': res_num, 'atom_name': atom_name}
            pocket_num = int(line[66:71])
            if pocket_num in pockets:
                pockets[pocket_num].append(atom_dict)
            else:
                pockets[pocket_num] = [atom_dict]

    return pockets