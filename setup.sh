mamba create -p /Users/maxim/Applications/BioRazer/conda-env python=3.10
mamba activate /Users/maxim/Applications/BioRazer/conda-env
mamba install \
    biopython \
    numpy \
    pandas \
    openpyxl \
    matplotlib \
    jupyter \
    pymol-open-source \
    biotite \
    mdtraj
cd /Users/maxim/Applications/PyRosetta/PyRosetta4.Release.python310.m1.release-388/setup
python setup.py install
pip install \
    impedance \
    umap
