# contactplot.py: plot interface residues delta sasa and burial percentage. #

## Requeriments ##
- Gnuplot
    * Microsoft Windows [download here](https://sourceforge.net/projects/gnuplot/files/gnuplot/5.2.2/)
    * Debian based Linux distributions: `sudo apt-get install gnuplot`
    * ArchLinux: `pacman -S gnuplot`
- Python 2.7 or 3.x
    * Install python dependencies: `pip install docopt biopython pandas`

## Usage ##
    contactplot.py <tablefile> <asafile> [--skip-none-contact --size <widthxheight> <outputfile>]

Where `<tablefile>` and `<asafile>` are the `.tsv` and `.asa` file respectively.
You may choose to plo only those residues which BSA > 0 using the 
`--skip-none-contact` option. By default `--size` option is 1920x1920
