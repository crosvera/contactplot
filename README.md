# contactplot.py: Solvent Accessible Surface Area plot software.

`contactplot.py` is a python script used to plot  the `dr-sasa`'s Solvent
Accessible Surface Area results.

## dr-sasa ##

`dr-sasa` is a solvent accessible surface area calculation software for 
biological molecules, that supports proteins, DNA and RNA inputs. The input 
files can be in the PDB or MOL2 format. PDB format files will use a NACCESS
compatible VdW radii lookup table, while MOL2 formats will use the same VdW
radii used in Chimera.


If you use this software in your research, please acknowledge it by citing the
following:

    Ribeiro, J., Ríos-Vera, C., Melo, F. and Schüller, A. (2018) 
    “Calculation of accurate contact surface areas between atoms for the 
    quantitative analysis of non-bonded molecular interactions”. 
    Bioinformatics (submitted).

## Install ##

### Requeriments 
- Python >= 2.7
- Matplotlib >= 2.1
- Pandas 
- Seaborn
- Doctopt
- Git (Optional, only if you are installing via pip)

### how to install

#### Using PIP:
Open your terminal (`cmd` in Windows) and run the following command:

    pip install git+https://github.com/crosvera/contactplot.git

### Usage

`coontactplot.py` has three operational modes:
- `residue`: Plots the Buried Surface Area (BSA or Delta SASA) matrix files
  calculated for residues (`*.by_res.tsv`) as a heatmap.
- `atom`: Plots the Buried Surface Area matrix file calculated for atoms 
  (`*._by_atom.tsv`) as a heatmap.
- `protein-ligand`: Plots the Buried Surface Area produced between a 
  protein-ligand interface using the BSA file calculated for atoms (`*._by_atom.tsv`)
  as heatmap.

Usage:

    contactplot.py residue <tablefile> <atmasafile> [--skip-non-contact <outputfile>]
    contactplot.py atom  <tablefile> <atmasafile> [--skip-non-contact <outputfile>]
    contactplot.py protein-ligand  <tablefile> <atmasafile> [--skip-non-contact <outputfile>]

Where:
- `--skip-non-contact` option is used to avoid ploting atoms or residues with
  none contact information.
- `<outputfile>` parameter if you want to change the output filename.

Example:

    #Plot using residue mode
    contactplot.py residue 3f3e_complex.LIGAND_vs_PROTEIN.by_res.tsv 3f3e_complex.atmasa --skip-non-contact

    #Plot using atom mode
    contactplot.py atom 3f3e_complex.LIGAND_vs_PROTEIN.by_atom.tsv 3f3e_complex.atmasa --skip-non-contact

    #Plot using protein-ligand mode
    contactplot.py atom 3f3e_complex.LIGAND_vs_PROTEIN.by_atom.tsv 3f3e_complex.atmasa --skip-non-contact

