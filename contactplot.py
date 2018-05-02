"""Usage:
    contactplot.py <contact_table> [--skip-non-contact --size <widthxheight> <outputfile>]

Options:
    --size <widthxheight>       Output image size in pixels.
                                [Default: 1024x900]
"""
from docopt import docopt
import pandas as pd
import subprocess
import tempfile
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


if __name__ == '__main__':
    args = docopt(__doc__)
    print(args)
    size = tuple(int(s) for s in args['--size'].split('x'))
    skip_non_contact = args['--skip-non-contact']
    table = args['<contact_table>']
    output = args['<outputfile>']
    if not output:
        output = os.path.splitext(table)[0] + '.png'

    contact_table_heatmap(table, output, skip_non_contact, size)
