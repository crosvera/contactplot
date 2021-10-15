#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Usage:
    contactplot.py residue <tablefile> <atmasafile> [--skip-non-contact <outputfile> --exclude-below=<below> --exclude-above=<above>]
    contactplot.py atom  <tablefile> <atmasafile> [--skip-non-contact <outputfile> --exclude-below=<below> --exclude-above=<above>]
    contactplot.py protein-ligand  <tablefile> <atmasafile> [--skip-non-contact <outputfile> --exclude-below=<below> --exclude-above=<above>]

"""
from docopt import docopt
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import pandas as pd
import seaborn as sns
import tempfile
from collections import OrderedDict
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from math import ceil
import gc


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


def get_atom_contact_area(tablefile, skip_none_contact=False):
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

    rename = lambda a: "{}/{}/{}/{}/{}".format(*a.split("/"))
    df.columns = [rename(c) for c in df.columns]
    iname = df.index.name
    df.index = [rename(i) for i in df.index]
    df.index.name = iname
    
    return df



def get_ligand_protein_contact_area(tablefile, skip_none_contact=False):
    with open(tablefile) as f:
        l = f.readline()
        columns = l.strip().split("\t")
    df = pd.read_csv(tablefile, sep="\t", index_col=0, skipinitialspace=True,
                     skip_blank_lines=True, usecols=columns)

    #check if protein is at column or index
    data = OrderedDict()
    iname = df.index.name
    if 'PROTEIN' in df.columns[0]:
        pivot = 'column'
        columns = []
        for col in df.columns:
            colname = '/'.join(col.split("/")[1:4])
            try:
                data[colname] += df[col]
            except KeyError:
                data[colname] = df[col]
                columns.append(colname)
        df = pd.DataFrame(data)[columns]
        df.index.name = iname

    elif 'PROTEIN' in df.index[0]:
        pivot = 'index'
        columns = df.columns
        indexes = []
        for idx in df.index:
            idxname = '/'.join(idx.split('/')[1:4])
            try:
                data[idxname] += df.loc[idx]
            except KeyError:
                data[idxname] = df.loc[idx]
                indexes.append(idxname)
        df = pd.DataFrame.from_dict(data, orient='index')[columns]
        #indexes.sort(key = lambda e: (e.split('/')[1], e.split('/')[2]))
        df = df.loc[indexes]
        df.index.name = iname


    if skip_none_contact:
        columns = df.columns
        df = df.drop([col for col in columns if df[col].sum() == 0.0], axis=1)
        indexes = df.index
        df = df.drop([idx for idx in indexes if df.loc[idx].sum() == 0.0])

    rename = lambda a: "{}/{}/{}/{}/{}".format(*a.split("/"))
    if pivot != 'column':
        df.columns = [rename(c) for c in df.columns]

    if pivot != 'index':
        iname = df.index.name
        df.index = [rename(i) for i in df.index]
        df.index.name = iname
    
    return df, pivot


def get_ligand_protein_bsa_vs_asa(tablefile, atmasafile, skip_none_contact=True):
    contacts_df, pivot = get_ligand_protein_contact_area(tablefile, skip_none_contact)
    ligand_bsa_asa = OrderedDict()
    protein_bsa_asa = OrderedDict()
    
    atmasa = pd.read_csv(atmasafile, sep='\t')
    atmasa['resnum'] = atmasa['resnum'].astype(str)
    atmasa['chain'] = atmasa['chain'].astype(str)
    atmasa_by_id = atmasa.set_index('ID')
    atom_total_asa = {str(atnum): atmasa_by_id['total_ASA'][atnum]\
                      for atnum in atmasa_by_id.index}
    atmasa_by_resnum = atmasa.groupby(['chain', 'resnum']).sum()
    res_total_asa = {resnum: atmasa_by_resnum['total_ASA'][resnum]\
                     for resnum in atmasa_by_resnum.index}
    
    protein_bsa = {}
    ligand_bsa = {}
    if pivot == 'column':
        for r in contacts_df.columns:
            rr = r.split('/')
            res = (rr[1], rr[2])
            protein_bsa[res] = contacts_df[r].sum()
        for a in contacts_df.index:
            aa = a.split('/')
            ligand_bsa[aa[-1]] = contacts_df.loc[a].sum()
            
    elif pivot == 'index':
        for r in contacts_df.index:
            rr = r.split('/')
            res = (rr[1], rr[2])
            protein_bsa[res] = contacts_df.loc[r].sum()
        for a in contacts_df.columns:
            aa = a.split('/')
            ligand_bsa[aa[-1]] = contacts_df[a].sum()
         
    ligand_total_asa = {atm: ligand_bsa[atm]+atom_total_asa[atm]\
                        for atm in ligand_bsa}
    protein_total_asa = {res: protein_bsa[res]+res_total_asa[res]\
                         for res in protein_bsa}
    
    for atm in ligand_total_asa:
        ligand_bsa_asa[atm] = ligand_bsa[atm] / ligand_total_asa[atm]
    for resnum in protein_total_asa:
        protein_bsa_asa[resnum] = protein_bsa[resnum] / protein_total_asa[resnum]
        
    return ligand_bsa_asa, protein_bsa_asa



def get_size(df, dpi):
    N = max(len(df.columns), len(df.index))
    W = max(max([len(c) for c in df.columns]),
            max([len(i) for i in df.index]))

    #print (N)
    if N>= 1000:
        sns.set_context('notebook', font_scale=0.5)
        factor = 0.65
    elif N > 700:
        sns.set_context('paper', font_scale=0.8)
        factor = 0.7
    elif N > 100:
        sns.set_context('talk')
        factor = 0.8
    else:
        sns.set_context('paper')
        factor = 1.8

    dsize = 18
    fontsizes = {
        'xx-small': (3/5.) * dsize,
        'x-small': (3/4.) * dsize,
        'small': (8/9.) * dsize,
        'medium': 1 * dsize,
        'large': (6/5.) * dsize,
        'x-large': (3/2.) * dsize,
        'xx-large': 2 * dsize,
    }


    #adjust size
    try:
        fontsize_pt = fontsizes[plt.rcParams['ytick.labelsize']]
        #print(1, fontsize_pt)
    except KeyError:
        fontsize_pt = plt.rcParams['ytick.labelsize']
        #print(2, fontsize_pt)
    matrix_height_pt = fontsize_pt * (N + W)
    matrix_height_in = matrix_height_pt / dpi

    matrix_weight_pt = fontsize_pt * (N + W)
    matrix_weight_in = matrix_weight_pt / dpi

    top_margin = 0.3 #if N >= 100 else 0.4
    button_margin = 0.3
    left_margin = 0.3
    right_margin = 0.2
    fig_height_in = matrix_height_in / (1 - top_margin - button_margin) * factor
    fig_weight_in = matrix_weight_in / (1 - left_margin - right_margin) * factor

    size = max(fig_weight_in, fig_height_in) *1.1
    inchsize = (int(ceil(size)), int(ceil(size)))
    #print(inchsize)
    return inchsize


def filter_out_area(df, below=None, above=None):
    if below is not None:
        df = df.apply(lambda x: np.where(x < below, np.nan, x))
    if above is not None:
        df = df.apply(lambda x: np.where(x > above, np.nan, x))

    return df
    


def plot_contact_ligand_protein(tablefile, atmasafile, output,
                                skip_none_contact=True, size=(1920,1920), dpi=72, filter_below=None, filter_above=None):
    df, pivot = get_ligand_protein_contact_area(tablefile, skip_none_contact)
    df = filter_out_area(df, filter_below, filter_above)

    ligand_bsa_asa, protein_bsa_asa = get_ligand_protein_bsa_vs_asa(tablefile, atmasafile, skip_none_contact)
 
    columns = df.columns
    indexes = df.index
    N = max(len(columns), len(indexes))

    structurename = tablefile.split(os.sep)[-1].split('.')[0]
    
    inchsize = get_size(df, dpi)
    fig = plt.figure(figsize=inchsize)
    size = inchsize[0]
    gc.collect()



    #cbar formatter:
    cbar_fmt = mtick.FuncFormatter(lambda x, pos: "{} $\AA^2$".format(x))


    #ax1 = plt.subplot2grid(inchsize, (0,1), colspan=inchsize[0]-2, rowspan=inchsize[1]-1)
    #ax2 = plt.subplot2grid(inchsize, (inchsize[0]-1,1), colspan=inchsize[0]-2, rowspan=1)
    ##ax3 = plt.subplot2grid(inchsize, (0,0), colspan=1, rowspan=inchsize[1]-1)
    #ax4 = plt.subplot2grid(inchsize, (0,inchsize[1]-1), colspan=1, rowspan=inchsize[1]-1)
    ax1 = plt.subplot2grid((20,20), (0,1), colspan=18, rowspan=19)
    ax2 = plt.subplot2grid((20,20), (19,1), colspan=18, rowspan=1)
    ax4 = plt.subplot2grid((20,20), (0,19), colspan=1, rowspan=19)

    sns.heatmap(df, ax=ax1, annot=False, xticklabels=False, yticklabels=True,
                cmap='YlOrRd', cbar_ax=ax4, cbar_kws=dict(format=cbar_fmt))
    ax1.set_yticklabels(df.index, rotation=0, size='x-large')
    ax4.tick_params(labelsize='x-large')

    if pivot == 'index':
        sns.barplot(x=list(df.columns), y=list(ligand_bsa_asa[e.split('/')[-1]]*100 for e in df.columns),
                    color='#005599', ax=ax2, label="BSA/ASA %")
    elif pivot == 'column':
        cols = []
        for c in df.columns:
            cc = c.split('/')
            cols.append( (cc[1], cc[2]) )
        sns.barplot(x=list(df.columns), y=list(protein_bsa_asa[c]*100 for c in cols),
                    color='#005599', ax=ax2, label="BSA/ASA %")

    ax2.set_xticklabels(list(df.columns), rotation=90, size='x-large')
    ax2.set_ylim([0, 100])
    ax2.set_yticks([])
    ax2t = ax2.twinx()
    ax2t.set_yticks([0, 100])
    ax2t.set_ylim([0, 100])
    ax2t.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f%%'))



    #set cbar outline 
    ax4.set_frame_on(True)

    X, Y = np.meshgrid(np.arange(0.5, len(columns)),
                       np.arange(0.5, len(indexes)))
    ax1.scatter(X, Y, color='gray', s=3)
    
    #legend = fig.legend(loc='lower left')
    #loc = (size*0.17/N, size*0.21/N)
    p = ax2.get_position()
    legend_size = 100.0 if N < 100 else 140.0
    x0 = ((p.x0 * size)  - legend_size/dpi)/ size
    loc = (x0, p.y0*1.11)
    legend = fig.legend(loc=loc)

    ax1.set_ylabel("")
    iname = df.index.name.split('/')[0]
    plt.suptitle("Buried surface area in structure %s.\n%s"%(structurename, iname), size='xx-large', y=0.92)



    fig.savefig(output, dpi=dpi)


def get_res_bsa_vs_asa(tablefile, atmasafile, skip_none_contact=True):
    contacts_df = get_res_contact_area(tablefile, skip_none_contact)
    res_bsa_asa_a = {}
    res_bsa_asa_b = {}
    


    atmasa = pd.read_csv(atmasafile, sep='\t')
    atmasa['resnum'] = atmasa['resnum'].astype(str)
    atmasa['chain'] = atmasa['chain'].astype(str)
    atmasa_by_id = atmasa.set_index('ID')
    atom_total_asa = {str(atnum): atmasa_by_id['total_ASA'][atnum]\
                      for atnum in atmasa_by_id.index}
    atmasa_by_resnum = atmasa.groupby(['chain', 'resnum']).sum()
    res_total_asa = {resnum: atmasa_by_resnum['total_ASA'][resnum]\
                     for resnum in atmasa_by_resnum.index}

    res_bsa_a = {}
    res_bsa_b = {}

    for r in contacts_df.columns:
        rr = r.split('/')
        res = (rr[1], rr[2])
        res_bsa_a[res] = contacts_df[r].sum()
    
    for r in contacts_df.index:
        rr = r.split('/')
        res = (rr[1], rr[2])
        res_bsa_b[res] = contacts_df.loc[r].sum()

    res_total_asa_a = {res: res_bsa_a[res]+res_total_asa[res] for res in res_bsa_a}
    res_total_asa_b = {res: res_bsa_b[res]+res_total_asa[res] for res in res_bsa_b}
    
    for res in res_total_asa_a:
        res_bsa_asa_a[res] = res_bsa_a[res] / res_total_asa_a[res]
    for res in res_total_asa_b:
        res_bsa_asa_b[res] = res_bsa_b[res] / res_total_asa_b[res]
        
    return res_bsa_asa_a, res_bsa_asa_b


def plot_contact_res_bsaasa(tablefile, atmasafile, output,
                            skip_none_contact=True, size=(1440,1400), dpi=72, filter_below=None, filter_above=None):
    res_bsa_asa_a, res_bsa_asa_b = get_res_bsa_vs_asa(tablefile, atmasafile,
                                                      skip_none_contact)

    df = get_res_contact_area(tablefile, skip_none_contact)
    df = filter_out_area(df, filter_below, filter_above)
    columns = df.columns
    indexes = df.index
    N = max(len(columns), len(indexes))

    structurename = tablefile.split(os.sep)[-1].split('.')[0]

    
    inchsize = get_size(df, dpi)
    fig = plt.figure(figsize=inchsize)
    size = inchsize[0]
    gc.collect()

    #cbar formatter:
    cbar_fmt = mtick.FuncFormatter(lambda x, pos: "{} $\AA^2$".format(x))


    #ax1 = plt.subplot2grid(inchsize, (0,1), colspan=inchsize[0]-2, rowspan=inchsize[1]-1)
    #ax2 = plt.subplot2grid(inchsize, (inchsize[0]-1,1), colspan=inchsize[0]-2, rowspan=1)
    ##ax3 = plt.subplot2grid(inchsize, (0,0), colspan=1, rowspan=inchsize[1]-1)
    #ax4 = plt.subplot2grid(inchsize, (0,inchsize[1]-1), colspan=1, rowspan=inchsize[1]-1)
    ax1 = plt.subplot2grid((20,20), (0,1), colspan=18, rowspan=19)
    ax2 = plt.subplot2grid((20,20), (19,1), colspan=18, rowspan=1)
    ax4 = plt.subplot2grid((20,20), (0,19), colspan=1, rowspan=19)

    sns.heatmap(df, ax=ax1, annot=False, xticklabels=False, yticklabels=True,
                cmap='YlOrRd', cbar_ax=ax4, cbar_kws=dict(format=cbar_fmt))
    ax1.set_yticklabels(df.index, rotation=0, size='x-large')
    ax4.tick_params(labelsize='x-large')

    cols = []
    for c in df.columns:
        cc = c.split('/')
        cols.append( (cc[1], cc[2]) )
    sns.barplot(x=list(df.columns), y=[res_bsa_asa_a[c]*100 for c in cols],
                color='#005599', ax=ax2, label="BSA/ASA %")

    ax2.set_xticklabels(list(df.columns), rotation=90, size='x-large')
    ax2.set_ylim([0, 100])
    ax2.set_yticks([])
    ax2t = ax2.twinx()
    ax2t.set_yticks([0, 100])
    ax2t.set_ylim([0, 100])
    ax2t.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f%%'))



    #set cbar outline 
    ax4.set_frame_on(True)

    X, Y = np.meshgrid(np.arange(0.5, len(columns)),
                       np.arange(0.5, len(indexes)))
    ax1.scatter(X, Y, color='gray', s=3)
    
    #legend = fig.legend(loc='lower left')
    #loc = (size*0.23/N, size*0.23/N)
    p = ax2.get_position()
    legend_size = 100.0 if N < 100 else 140.0
    x0 = ((p.x0 * size)  - legend_size/dpi)/ size
    loc = (x0, p.y0*1.11)
    legend = fig.legend(loc=loc)

    ax1.set_ylabel("")
    iname = df.index.name.split('/')[0]
    plt.suptitle("Buried surface area in structure %s.\n%s"%(structurename, iname), size='xx-large', y=0.92)
    fig.savefig(output, dpi=dpi)


def get_atom_bsa_vs_asa(tablefile, atmasafile, skip_none_contact=True):
    contacts_df = get_atom_contact_area(tablefile, skip_none_contact)
    atom_bsa_asa_a = {}
    atom_bsa_asa_b = {}
    
    atmasa = pd.read_csv(atmasafile, sep='\t')
    atmasa_by_id = atmasa.set_index('ID')
    atom_total_asa = {str(atnum): atmasa_by_id['total_ASA'][atnum] for atnum in atmasa_by_id.index}
    atmasa_by_resnum = atmasa.groupby(['resnum']).sum()
    res_total_asa = {str(resnum): atmasa_by_resnum['total_ASA'][resnum] for resnum in atmasa_by_resnum.index}
    
    atom_bsa_a = {atom.split('/')[-1]:contacts_df[atom].sum()\
                 for atom in contacts_df.columns}
    atom_bsa_b = {atom.split('/')[-1]:contacts_df.loc[atom].sum()\
                 for atom in contacts_df.index}
    
    atom_total_asa_a = {atom: atom_bsa_a[atom]+atom_total_asa[atom] for atom in atom_bsa_a}
    atom_total_asa_b = {atom: atom_bsa_b[atom]+atom_total_asa[atom] for atom in atom_bsa_b}
    
    for atom in atom_total_asa_a:
        atom_bsa_asa_a[atom] = atom_bsa_a[atom] / atom_total_asa_a[atom]
    for atom in atom_total_asa_b:
        atom_bsa_asa_b[atom] = atom_bsa_b[atom] / atom_total_asa_b[atom]
        
    return atom_bsa_asa_a, atom_bsa_asa_b


def plot_contact_atom_bsaasa(tablefile, atmasafile, output,
                             skip_none_contact=True, size=(1920,1920), dpi=72, filter_below=None, filter_above=None):
    atom_bsa_asa_a, atom_bsa_asa_b = get_atom_bsa_vs_asa(tablefile, atmasafile,
                                                         skip_none_contact)

    df = get_atom_contact_area(tablefile, skip_none_contact)
    df = filter_out_area(df, filter_below, filter_above)
    columns = df.columns
    indexes = df.index
    N = max(len(columns), len(indexes))

    structurename = tablefile.split(os.sep)[-1].split('.')[0]

    inchsize = get_size(df, dpi)
    fig = plt.figure(figsize=inchsize)
    size = inchsize[0]
    gc.collect()

    #cbar formatter:
    cbar_fmt = mtick.FuncFormatter(lambda x, pos: "{} $\AA^2$".format(x))


    #ax1 = plt.subplot2grid(inchsize, (0,1), colspan=inchsize[0]-2, rowspan=inchsize[1]-1)
    #ax2 = plt.subplot2grid(inchsize, (inchsize[0]-1,1), colspan=inchsize[0]-2, rowspan=1)
    ##ax3 = plt.subplot2grid(inchsize, (0,0), colspan=1, rowspan=inchsize[1]-1)
    #ax4 = plt.subplot2grid(inchsize, (0,inchsize[1]-1), colspan=1, rowspan=inchsize[1]-1)
    ax1 = plt.subplot2grid((20,20), (0,1), colspan=18, rowspan=19)
    ax2 = plt.subplot2grid((20,20), (19,1), colspan=18, rowspan=1)
    ax4 = plt.subplot2grid((20,20), (0,19), colspan=1, rowspan=19)

    sns.heatmap(df, ax=ax1, annot=False, xticklabels=False, yticklabels=True,
                cmap='YlOrRd', cbar_ax=ax4, cbar_kws=dict(format=cbar_fmt))
    ax1.set_yticklabels(df.index, rotation=0, size='x-large')
    ax4.tick_params(labelsize='x-large')

    sns.barplot(x=list(df.columns), y=list(atom_bsa_asa_a[e.split('/')[-1]]*100 for e in df.columns),
                color='#005599', ax=ax2, label="BSA/ASA %")

    ax2.set_xticklabels(list(df.columns), rotation=90, size='x-large')
    ax2.set_ylim([0, 100])
    ax2.set_yticks([])
    ax2t = ax2.twinx()
    ax2t.set_yticks([0, 100])
    ax2t.set_ylim([0, 100])
    ax2t.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f%%'))



    #set cbar outline 
    ax4.set_frame_on(True)

    X, Y = np.meshgrid(np.arange(0.5, len(columns)),
                       np.arange(0.5, len(indexes)))
    ax1.scatter(X, Y, color='gray', s=3)
    
    p = ax2.get_position()
    legend_size = 100.0 if N < 100 else 140.0
    x0 = ((p.x0 * size)  - legend_size/dpi)/ size
    loc = (x0, p.y0*1.11)
    legend = fig.legend(loc=loc)

    ax1.set_ylabel("")
    iname = df.index.name.split('/')[0]
    plt.suptitle("Buried surface area in structure %s.\n%s"%(structurename, iname), size='xx-large', y=0.92)
    fig.savefig(output, dpi=dpi)


def check_tablefile_type(tablefile):
    file_split = tablefile.split('.')
    table_elements = file_split[-2]
    table_type = file_split[-4]

    return table_elements, table_type




def contactplot():
    args = docopt(__doc__)
    skip_none_contact = args['--skip-non-contact']
    tablefile = args['<tablefile>']
    atmasafile = args['<atmasafile>']
    output = args['<outputfile>']
    below = float(args['--exclude-below']) if args['--exclude-below']  else None
    above = float(args['--exclude-above']) if args['--exclude-above']  else None
    if not output:
        if args['protein-ligand']:
            output = os.path.splitext(tablefile)[0] + '_protein_ligand.png'
        else:
            output = os.path.splitext(tablefile)[0] + '.png'

    if atmasafile.split('.')[-1] != 'atmasa':
        print("ERROR: {} is not the .atmasa file.".format(atmasafile))
        sys.exit(-1)

    if args['residue']:
        if check_tablefile_type(tablefile)[0] != 'by_res':
            print ("ERROR: Incorrect tablefile, you should use the *.by_res.tsv file.")
            sys.exit(-1)
        plot_contact_res_bsaasa(tablefile, atmasafile, output, skip_none_contact, filter_below=below, filter_above=above)
    elif args['atom']:
        if check_tablefile_type(tablefile)[0] != 'by_atom':
            print ("ERROR: Incorrect tablefile, you should use the *.by_atom.tsv file.")
            sys.exit(-1)
        plot_contact_atom_bsaasa(tablefile, atmasafile, output, skip_none_contact, filter_below=below, filter_above=above)
    elif args['protein-ligand']:
        if check_tablefile_type(tablefile)[0] != 'by_atom':
            print ("ERROR: Incorrect tablefile, you should use the *.by_atom.tsv file.")
            sys.exit(-1)
        if check_tablefile_type(tablefile)[1] == 'matrix':
            print ("ERROR: Table file type 'matrix' still not supported in this mode.")
            sys.exit(-1)
        plot_contact_ligand_protein(tablefile, atmasafile, output, skip_none_contact, filter_below=below, filter_above=above)


if __name__ == '__main__':
    contactplot()
