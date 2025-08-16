from Bio.PDB import MMCIFParser, PDBIO
from pathlib import Path

def convert_cif_to_pdb(cif_file, pdb_file):
    sstruc = MMCIFParser().get_structure("", str(cif_file))
    io = PDBIO()
    io.set_structure(sstruc)
    io.save(str(pdb_file))

def batch_convert(src_dir, src_file_format, target_dir, convert_func):
    src_dir = Path(src_dir)
    target_dir = Path(target_dir)
    if not src_dir.is_dir():
        raise ValueError(f"Source directory {src_dir} does not exist or is not a directory.")
    if not target_dir.exists():
        target_dir.mkdir(parents=True, exist_ok=True)

    for src_file in src_dir.glob(f"*.{src_file_format}"):
        target_file = target_dir / (src_file.stem + ".pdb")
        convert_func(src_file, target_file)
