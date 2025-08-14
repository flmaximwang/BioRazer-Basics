from pathlib import Path
import pymol

def align_models_to_ref(
    dir_of_models,
    ref_pdb_path,
    model_pymol_selector='all',
    ref_pymol_selector='all',
    target_dir=None,
    verbose=True
):
    dir_of_models = Path(dir_of_models)
    ref_pdb_path = Path(ref_pdb_path)
    if not target_dir:
        target_dir = dir_of_models.parent / (dir_of_models.name + "_aligned")
    if not target_dir.exists():
        target_dir.mkdir()
    for model in dir_of_models.iterdir():
        if not model.suffix in [".pdb", ".cif"]:
            continue
        if verbose:
            print(f">>> Aligning {model} to {ref_pdb_path}")
        pymol.cmd.reinitialize()
        pymol.cmd.load(model, "model_to_move")
        pymol.cmd.load(ref_pdb_path, "ref_model")
        pymol.cmd.align(f"model_to_move and ({model_pymol_selector})", f"ref_model and ({ref_pymol_selector})")
        pymol.cmd.save(target_dir / model.name, "model_to_move")