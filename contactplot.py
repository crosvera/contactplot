"""Usage:
    contactplot.py <tablefile> <asafile> [--skip-none-contact --size <widthxheight> <outputfile>]

Options:
    --size <widthxheight>       Output image size in pixels.
                                [Default: 1920x1920]
"""
from docopt import docopt
from Bio.PDB import PDBParser
import pandas as pd
import subprocess
import tempfile
from collections import OrderedDict
import os

gnuplot_heatmap = """
set terminal pngcairo size {size[0]},{size[1]} enhanced font 'Verdana,12'
set output '{output}'
unset key

set view map
#set pm3d interpolate 0,0
#set pm3d map
#set palette negative grey
#set palette rgb -23,-30, -1

set datafile separator ","
set datafile missing "nan"


# line styles
set style line 1 lt 1 lc rgb '#FFFFCC' # very light yellow-orange-red
set style line 2 lt 1 lc rgb '#FFEDA0' # 
set style line 3 lt 1 lc rgb '#FED976' # light yellow-orange-red
set style line 4 lt 1 lc rgb '#FEB24C' # 
set style line 5 lt 1 lc rgb '#FD8D3C' # 
set style line 6 lt 1 lc rgb '#FC4E2A' # medium yellow-orange-red
set style line 7 lt 1 lc rgb '#E31A1C' #
set style line 8 lt 1 lc rgb '#B10026' # dark yellow-orange-red

# palette
set palette defined ( 0 '#FFFFCC',\\
                      1 '#FFEDA0',\\
                      2 '#FED976',\\
                      3 '#FEB24C',\\
                      4 '#FD8D3C',\\
                      5 '#FC4E2A',\\
                      6 '#E31A1C',\\
                      7 '#B10026' )





set border linewidth 1
set grid

#set lmargin screen 0.1
#set rmargin screen 0.9
#set tmargin screen 0.9
#set bmargin screen 0.1

set xtics nomirror rotate by -90 font "Verdana,10" offset character 0,-1,0 autojustify
set ytics nomirror font "Verdana,10" offset character 0,0,0 autojustify

#set cbrange [0:1]
#unset colorbox

#set xlabel "\\n\\nChain A"
#set ylabel "Chain B\\n"

set title "{structure} buried surface heatmap"



splot "{input}" using 1:2:3:xtic(4):ytic(5) with image
set xr [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
replot
"""


gnuplot_bsaasa = """
set terminal pngcairo size {size[0]},{size[1]} enhanced font 'Verdana,10'
set output '{output}'
unset key

set multiplot title "Contact map from {structure}.\\nBuried surface in residues of Chain {chains[0]} caused by residues of Chain {chains[1]}" font "Verdana,18"

set view map
set pm3d

set datafile separator ","
set datafile missing "nan"


# line styles
set style line 1 lt 1 lc rgb '#FFFFCC' # very light yellow-orange-red
set style line 2 lt 1 lc rgb '#FFEDA0' # 
set style line 3 lt 1 lc rgb '#FED976' # light yellow-orange-red
set style line 4 lt 1 lc rgb '#FEB24C' # 
set style line 5 lt 1 lc rgb '#FD8D3C' # 
set style line 6 lt 1 lc rgb '#FC4E2A' # medium yellow-orange-red
set style line 7 lt 1 lc rgb '#E31A1C' #
set style line 8 lt 1 lc rgb '#B10026' # dark yellow-orange-red

# palette
set palette defined ( 0 '#FFFFCC',\
                      1 '#FFEDA0',\
                      2 '#FED976',\
                      3 '#FEB24C',\
                      4 '#FD8D3C',\
                      5 '#FC4E2A',\
                      6 '#E31A1C',\
                      7 '#B10026' )




set border linewidth 0.1
set grid

set lmargin at screen 0.1305
set rmargin at screen 0.95
set bmargin at screen 0.151
set tmargin at screen 0.95


set colorbox user origin .96, .160 size .015, .790
set cbtics font "Verdana, 14"

unset xtics
unset ytics



splot "{contacts}" using 1:2:3:xtic(4):ytic(5) with image, \
    "{contacts}" using 1:2:3:xtic(4):ytic(5) with dots

set xr [-0.5:GPVAL_DATA_X_MAX+0.5]
rangex = GPVAL_DATA_X_MAX
set yrange [-0.5: GPVAL_DATA_Y_MAX-0.5]
replot

unset pm3d
unset key
set xrange [*:*]
set style fill transparent solid 0.5 

set lmargin at screen 0.08
set rmargin at screen 0.13

set ylabel "Chain {chains[1]}"
set ylabel font "Verdana,14"
set ytics nomirror font "Verdana,10" offset character 0,0,0 autojustify

#set x2tics 0,0.5,1 offset -1
set x2tics 0,50,100 offset -1 format "%1.0f%%"

plot "{bsaasa_b}" using ($3*0.5*100):1:($3*0.5*100):(0.4):ytic(2) with boxxyerrorbars lc rgb 0x5599 axes x2y1

set lmargin at screen 0.13
set rmargin at screen 0.95
set bmargin at screen 0.10
set tmargin at screen 0.15

set xtics nomirror rotate by -90 font "Verdana,10" offset character 0,0,0 autojustify
unset ytics
unset ylabel
set yrange [*:*]
unset arrow 1
#set arrow 2 from -0.5,0.95 to GPVAL_DATA_X_MAX-0.5,0.95 nohead

set key at -1, 0.5

set xlabel "Chain {chains[0]}"
set xlabel font "Verdana,14"
set xrange [-0.5: rangex+0.5]
unset x2tics
#set y2tics 0,0.5,1 offset -1
set y2tics 0,50,100 format "%1.0f%%"
plot "{bsaasa_a}" using 1:($3*100):xtic(2) with boxes lc rgb 0x5599 title "Burial percentage\\n(BSA/ASA)" axes x1y2

"""


def contact_table_heatmap(tablefile, outputfile, skip_non_contact=False, size=(1024,900)):
    with open(tablefile) as f:
        l = f.readline()
        columns = l.strip().split("\t")
    df = pd.read_csv(tablefile, sep="\t", index_col=0, skipinitialspace=True,
                     skip_blank_lines=True, usecols=columns)

    if skip_non_contact:
        columns = df.columns
        df = df.drop([col for col in columns if df[col].sum() == 0.0], axis=1)
        indexes = df.index
        df = df.drop([idx for idx in indexes if df.loc[idx].sum() == 0.0])

    columns = df.columns
    indexes = df.index
    inputf = tempfile.NamedTemporaryFile('w')
    scriptf = tempfile.NamedTemporaryFile('w')
    for i, col in enumerate(columns):
        for j, idx in enumerate(indexes):
            l = "%d, %d, %f, %s, %s \n" % (i, j, float(df[col][idx]), col, idx)
            inputf.file.write(l)
    inputf.file.close()
    
    script_parameters = dict(input=inputf.name, output=output, size=size,
                             structure=tablefile.split('.')[0])
    scriptf.file.write(gnuplot_heatmap.format(**script_parameters))
    scriptf.file.close()
    p = subprocess.Popen('gnuplot %s'%scriptf.name, shell=True)
    os.waitpid(p.pid, 0)
    inputf.close()
    scriptf.close()


def get_res_contact_area(tablefile, skip_none_contact=False):
    with open(tablefile) as f:
        l = f.readline()
        columns = l.strip().split("\t")
    df = pd.read_csv(tablefile, sep="\t", index_col=0, skipinitialspace=True,
                     skip_blank_lines=True, usecols=columns)

    if skip_none_contact:
        columns = df.columns
        df = df.drop([col for col in columns if df[col].sum() == 0.0], axis=1)
        indexes = df.index
        df = df.drop([idx for idx in indexes if df.loc[idx].sum() == 0.0])

    return df

def get_res_asa(asafile):
    parser = PDBParser()
    structure_name = 'structure'
    structure = parser.get_structure(structure_name, asafile)
    model = structure[0]
    
    res_asa = {}
    for res in model.get_residues():
        #print((res.get_resname(), res.parent.get_id(), res.get_id()[1]))
        resname = '%s/%s/%d' % (res.get_resname(), res.parent.get_id(), res.get_id()[1])
        res_asa[resname] = sum([atom.get_bfactor() for atom in res.get_atoms()])
        
    return res_asa


def get_res_bsa_vs_asa(tablefile, asafile, skip_none_contact=True):
    contacts_df = get_res_contact_area(tablefile, skip_none_contact)
    res_bsa_a = OrderedDict()
    res_bsa_b = OrderedDict()
    for res in contacts_df.columns:
        res_bsa_a[res] = contacts_df[res].sum()
    for res in contacts_df.index:
        res_bsa_b[res] = contacts_df[:res].sum().sum()

    res_asa = get_res_asa(asafile)
    total_res_asa = {}
    for res in res_bsa_a:
        total_res_asa[res] = res_asa[res] + res_bsa_a[res]
    for res in res_bsa_b:
        total_res_asa[res] = res_asa[res] + res_bsa_b[res]

    res_bsa_asa_a = OrderedDict((res, float(res_bsa_a[res]/total_res_asa[res]))\
                                for res in res_bsa_a)
    res_bsa_asa_b = OrderedDict((res, float(res_bsa_b[res]/total_res_asa[res]))\
                                for res in res_bsa_b)
    return res_bsa_asa_a, res_bsa_asa_b


def get_chains(contacts_df):
    index_name = contacts_df.index.name
    if '->' in index_name:
        chainb, chaina = index_name.split('->')
    elif '<-' in index_name:
        chaina, chainb = index_name.split('<-')
    else:
        chaina = chainb = ''
    return chaina, chainb


def plot_contact_bsaasa(tablefile, asafile, output, 
                        skip_none_contact=True, size=(1024,900)):
    res_bsa_asa_a, res_bsa_asa_b = get_res_bsa_vs_asa(tablefile, asafile,
                                                      skip_none_contact)
    chaina_input = tempfile.NamedTemporaryFile('w')
    for i, res in enumerate(res_bsa_asa_a):
        l = '%d, %s, %f\n' % (i, res, res_bsa_asa_a[res])
        chaina_input.file.write(l)
    chaina_input.file.close()

    chainb_input = tempfile.NamedTemporaryFile('w')
    for i, res in enumerate(res_bsa_asa_b):
        l = '%d, %s, %f\n' % (i, res, res_bsa_asa_b[res])
        chainb_input.file.write(l)
    chainb_input.file.close()

    df = get_res_contact_area(tablefile, skip_none_contact)
    columns = df.columns
    indexes = df.index
    inputf = tempfile.NamedTemporaryFile('w')
    scriptf = tempfile.NamedTemporaryFile('w')
    for i, col in enumerate(columns):
        for j, idx in enumerate(indexes):
            l = "%d, %d, %f, %s, %s \n" % (i, j, float(df[col][idx]), col, idx)
            inputf.file.write(l)
    inputf.file.close()

    structurename = tablefile.split(os.sep)[-1].split('.')[0]
    
    chaina, chainb = get_chains(df)
    script_parameters = dict(contacts=inputf.name, output=output, size=size,
                             structure=structurename,
                             bsaasa_a=chaina_input.name,
                             bsaasa_b=chainb_input.name,
                             chains = (chaina, chainb)
                             )
    scriptf.file.write(gnuplot_bsaasa.format(**script_parameters))
    scriptf.file.close()
    p = subprocess.Popen('gnuplot %s'%scriptf.name, shell=True)
    os.waitpid(p.pid, 0)
    inputf.close()
    scriptf.close()



if __name__ == '__main__':
    args = docopt(__doc__)
    print(args)
    size = tuple(int(s) for s in args['--size'].split('x'))
    skip_none_contact = args['--skip-none-contact']
    tablefile = args['<tablefile>']
    asafile = args['<asafile>']
    output = args['<outputfile>']
    if not output:
        output = os.path.splitext(tablefile)[0] + '.png'

    plot_contact_bsaasa(tablefile, asafile, output, skip_none_contact, size)
