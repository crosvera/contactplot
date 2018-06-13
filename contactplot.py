#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Usage:
    contactplot.py residue <tablefile> <asafile> [--skip-none-contact --size <widthxheight> <outputfile>]
    contactplot.py atom  <tablefile> <asafile> [--skip-none-contact --size <widthxheight> <outputfile>]
    contactplot.py ligand-protein  <tablefile> <asafile> [--skip-none-contact --size <widthxheight> <outputfile>]

Options:
    --size <widthxheight>       Output image size in pixels.
                                [Default: 1920x1920]
"""
from docopt import docopt
from Bio.PDB import PDBParser
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import pandas as pd
import seaborn as sns
import tempfile
from collections import OrderedDict
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


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
        resname = '%s/%s/%d' % (res.get_resname().strip(), res.parent.get_id(), res.get_id()[1])
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


def plot_contact_res_bsaasa(tablefile, asafile, output,
                         skip_none_contact=True, size=(1440,1400), dpi=72):
    res_bsa_asa_a, res_bsa_asa_b = get_res_bsa_vs_asa(tablefile, asafile,
                                                      skip_none_contact)

    df = get_res_contact_area(tablefile, skip_none_contact)
    columns = df.columns
    indexes = df.index


    structurename = tablefile.split(os.sep)[-1].split('.')[0]
    
    inchsize = (int(size[0]/dpi), int(size[1]/dpi))
    fig = plt.figure(figsize=inchsize)
    ax1 = plt.subplot2grid(inchsize, (0,1), colspan=inchsize[0]-2, rowspan=inchsize[1]-1)
    ax2 = plt.subplot2grid(inchsize, (inchsize[0]-1,1), colspan=inchsize[0]-2, rowspan=1)
    ax3 = plt.subplot2grid(inchsize, (0,0), colspan=1, rowspan=inchsize[1]-1)
    ax4 = plt.subplot2grid(inchsize, (0,inchsize[1]-1), colspan=1, rowspan=inchsize[1]-1)
    
    #cbar formatter:
    cbar_fmt = mtick.FuncFormatter(lambda x, pos: "{} $\AA^2$".format(x))

    sns.heatmap(df, ax=ax1, annot=False, xticklabels=False, yticklabels=False,
                cmap='YlOrRd', cbar_ax=ax4, cbar_kws=dict(format=cbar_fmt))
    #set cbar outline 
    ax4.set_frame_on(True)

    X, Y = np.meshgrid(np.arange(0.5, len(columns)),
                       np.arange(0.5, len(indexes)))
    ax1.scatter(X, Y, color='gray', s=3)

    sns.barplot(x=list(df.columns), y=list(res_bsa_asa_a[e]*100 for e in df.columns),
                color='#005599', ax=ax2, label="BSA/ASA %")
    ax4.tick_params(labelsize='x-large')
    ax2.set_xticklabels(list(df.columns), rotation=90, size='x-large')
    ax2.set_ylim([0, 100])
    ax2.set_yticks([])
    ax2t = ax2.twinx()
    ax2t.set_yticks([0, 100])
    ax2t.set_ylim([0, 100])
    ax2t.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f%%'))
    legend = fig.legend(loc=(0.09,0.11))
    
    sns.barplot(y=list(df.index), x=list(res_bsa_asa_b[e]*100 for e in df.index), color='#005599', ax=ax3)
    ax3.set_yticklabels(df.index, size='x-large')
    ax3.set_xlim([0, 100])
    ax3.set_xticks([])
    ax3t = ax3.twiny()
    ax3t.set_xticks([0, 100])
    ax3t.set_xlim([0, 100])
    ax3t.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f%%'))


    ax1.set_ylabel("")
    iname = df.index.name.split('/')[0]
    plt.suptitle("Buried surface area in structure %s.\n%s"%(structurename, iname), size='xx-large', y=0.92)
    fig.savefig(output, dpi=dpi)



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
    data = {}
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
        indexes.sort(key = lambda e: (int(e.split('/')[2]), e.split('/')[1]))
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



def get_atom_asa(asafile):
    parser = PDBParser()
    structure_name = 'structure'
    structure = parser.get_structure(structure_name, asafile)
    model = structure[0]
    
    atom_asa = {}
    for a in model.get_atoms():
        aname = "%s/%s/%s/%d/%d" % (a.get_id(), a.parent.get_resname(), 
                                    a.parent.parent.get_id(), a.parent.get_id()[1], 
                                    a.get_serial_number())
        atom_asa[aname] = a.get_bfactor()
        
    return atom_asa

def get_atom_bsa_vs_asa(tablefile, asafile, skip_none_contact=True):
    contacts_df = get_atom_contact_area(tablefile, skip_none_contact)
    atom_bsa_a = OrderedDict()
    atom_bsa_b = OrderedDict()
    for atom in contacts_df.columns:
        atom_bsa_a[atom] = contacts_df[atom].sum()
    for atom in contacts_df.index:
        atom_bsa_b[atom] = contacts_df[:atom].sum().sum()
        
    atom_asa = get_atom_asa(asafile)
    total_atom_asa = {}
    for atom in atom_bsa_a:
        total_atom_asa[atom] = atom_asa[atom] + atom_bsa_a[atom]
    for atom in atom_bsa_b:
        total_atom_asa[atom] = atom_asa[atom] + atom_bsa_b[atom]
        
    atom_bsa_asa_a = OrderedDict((atom, float(atom_bsa_a[atom]/total_atom_asa[atom]))\
                                 for atom in atom_bsa_a)
    atom_bsa_asa_b = OrderedDict((atom, float(atom_bsa_b[atom]/total_atom_asa[atom]))\
                                 for atom in atom_bsa_b)
    
    return atom_bsa_asa_a, atom_bsa_asa_b


def get_ligand_protein_bsa_vs_asa(tablefile, asafile, skip_none_contact=True):
    contacts_df, pivot = get_ligand_protein_contact_area(tablefile, skip_none_contact)
    ligand_bsa = OrderedDict()
    protein_bsa = OrderedDict()
    total_ligand_asa = {}
    total_protein_asa = {}

    atom_asa = get_atom_asa(asafile)
    res_asa = get_res_asa(asafile)

    if pivot == 'column':
        for res in contacts_df.columns:
            protein_bsa[res] = contacts_df[res].sum()
        for atom in contacts_df.index:
            ligand_bsa[atom] = contacts_df.loc[atom].sum()
    elif pivot == 'index':
        for res in contacts_df.index:
            protein_bsa[res] = contacts_df.loc[res].sum()
        for atom in contacts_df.columns:
            ligand_bsa[atom] = contacts_df[atom].sum()

    for res in protein_bsa:
        total_protein_asa[res] = res_asa[res] + protein_bsa[res]
    for atom in ligand_bsa:
        total_ligand_asa[atom] = atom_asa[atom] + ligand_bsa[atom]

    protein_bsa_asa = OrderedDict((res, float(protein_bsa[res]/total_protein_asa[res]))\
                                  for res in protein_bsa)
    ligand_bsa_asa = OrderedDict((atom, float(ligand_bsa[atom]/total_ligand_asa[atom]))\
                                 for atom in ligand_bsa)

    return ligand_bsa_asa, protein_bsa_asa



def plot_contact_atom_bsaasa(tablefile, asafile, output,
                             skip_none_contact=True, size=(1920,1920), dpi=72):
    atom_bsa_asa_a, atom_bsa_asa_b = get_atom_bsa_vs_asa(tablefile, asafile,
                                                         skip_none_contact)

    df = get_atom_contact_area(tablefile, skip_none_contact)
    columns = df.columns
    indexes = df.index


    structurename = tablefile.split(os.sep)[-1].split('.')[0]
    
    inchsize = (int(size[0]/dpi), int(size[1]/dpi))
    fig = plt.figure(figsize=inchsize)
    ax1 = plt.subplot2grid(inchsize, (0,1), colspan=inchsize[0]-2, rowspan=inchsize[1]-1)
    ax2 = plt.subplot2grid(inchsize, (inchsize[0]-1,1), colspan=inchsize[0]-2, rowspan=1)
    ax3 = plt.subplot2grid(inchsize, (0,0), colspan=1, rowspan=inchsize[1]-1)
    ax4 = plt.subplot2grid(inchsize, (0,inchsize[1]-1), colspan=1, rowspan=inchsize[1]-1)

    #cbar formatter:
    cbar_fmt = mtick.FuncFormatter(lambda x, pos: "{} $\AA^2$".format(x))

    sns.heatmap(df, ax=ax1, annot=False, xticklabels=False, yticklabels=False,
                cmap='YlOrRd', cbar_ax=ax4, cbar_kws=dict(format=cbar_fmt))
    #set cbar outline 
    ax4.set_frame_on(True)
    ax4.tick_params(labelsize='x-large')

    X, Y = np.meshgrid(np.arange(0.5, len(columns)),
                       np.arange(0.5, len(indexes)))
    ax1.scatter(X, Y, color='gray', s=3)

    sns.barplot(x=list(df.columns), y=list(atom_bsa_asa_a[e]*100 for e in df.columns),
                color='#005599', ax=ax2, label="BSA/ASA %")
    ax2.set_xticklabels(list(df.columns), rotation=90, size='x-large')
    ax2.set_ylim([0, 100])
    ax2.set_yticks([])
    ax2t = ax2.twinx()
    ax2t.set_yticks([0, 100])
    ax2t.set_ylim([0, 100])
    ax2t.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f%%'))
    legend = fig.legend(loc=(0.09,0.11))
    
    sns.barplot(y=list(df.index), x=list(atom_bsa_asa_b[e]*100 for e in df.index), color='#005599', ax=ax3)
    ax3.set_yticklabels(df.index, size='x-large')
    ax3.set_xlim([0, 100])
    ax3.set_xticks([])
    ax3t = ax3.twiny()
    ax3t.set_xticks([0, 100])
    ax3t.set_xlim([0, 100])
    ax3t.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f%%'))


    ax1.set_ylabel("")
    iname = df.index.name.split('/')[0]
    plt.suptitle("Buried surface area in structure %s.\n%s"%(structurename, iname), fontsize='xx-large', y=0.92)
    fig.savefig(output, dpi=dpi)



def plot_contact_ligand_protein(tablefile, asafile, output,
                                skip_none_contact=True, size=(1920,1920), dpi=72):
    df, pivot = get_ligand_protein_contact_area(tablefile, skip_none_contact)

    ligand_bsa_asa, protein_bsa_asa = get_ligand_protein_bsa_vs_asa(tablefile,
                                                                    asafile,
                                                                    skip_none_contact)

    columns = df.columns
    indexes = df.index

    structurename = tablefile.split(os.sep)[-1].split('.')[0]
    
    inchsize = (int(size[0]/dpi), int(size[1]/dpi))
    fig = plt.figure(figsize=inchsize)



    #cbar formatter:
    cbar_fmt = mtick.FuncFormatter(lambda x, pos: "{} $\AA^2$".format(x))


    ax1 = plt.subplot2grid(inchsize, (0,1), colspan=inchsize[0]-2, rowspan=inchsize[1]-1)
    ax2 = plt.subplot2grid(inchsize, (inchsize[0]-1,1), colspan=inchsize[0]-2, rowspan=1)
    #ax3 = plt.subplot2grid(inchsize, (0,0), colspan=1, rowspan=inchsize[1]-1)
    ax4 = plt.subplot2grid(inchsize, (0,inchsize[1]-1), colspan=1, rowspan=inchsize[1]-1)

    sns.heatmap(df, ax=ax1, annot=False, xticklabels=False, yticklabels=True,
                cmap='YlOrRd', cbar_ax=ax4, cbar_kws=dict(format=cbar_fmt))
    ax1.set_yticklabels(df.index, rotation=0, size='x-large')
    ax4.tick_params(labelsize='x-large')

    if pivot == 'index':
        sns.barplot(x=list(df.columns), y=list(ligand_bsa_asa[e]*100 for e in df.columns),
                    color='#005599', ax=ax2, label="BSA/ASA %")
    elif pivot == 'column':
        sns.barplot(x=list(df.columns), y=list(protein_bsa_asa[e]*100 for e in df.columns),
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
    
    legend = fig.legend(loc=(0.09,0.11))

    ax1.set_ylabel("")
    iname = df.index.name.split('/')[0]
    plt.suptitle("Buried surface area in structure %s.\n%s"%(structurename, iname), size='xx-large', y=0.92)
    fig.savefig(output, dpi=dpi)






def contactplot():
    args = docopt(__doc__)
    size = tuple(int(s) for s in args['--size'].split('x'))
    skip_none_contact = args['--skip-none-contact']
    tablefile = args['<tablefile>']
    asafile = args['<asafile>']
    output = args['<outputfile>']
    if not output:
        output = os.path.splitext(tablefile)[0] + '.png'

    if args['residue']:
        plot_contact_res_bsaasa(tablefile, asafile, output, skip_none_contact, size)
    elif args['atom']:
        plot_contact_atom_bsaasa(tablefile, asafile, output, skip_none_contact, size)
    elif args['ligand-protein']:
        plot_contact_ligand_protein(tablefile, asafile, output, skip_none_contact, size)


if __name__ == '__main__':
    contactplot()
