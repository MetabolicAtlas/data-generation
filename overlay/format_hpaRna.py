#!/usr/bin/env python

# Description: Format RNA HPA tissue gene data so that it is usable for dataOverlay in MetAtlas

# Raw data can be downloaded from https://www.proteinatlas.org/download/rna_tissue_hpa.tsv.zip

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
import yaml

def scalingfunc(x, U = 10):
    return min(1, np.log2(x + 1)/U)

def logfunc(x):
    return np.log2(x+1)

def calMax(qt):
    """Estimate the upper whisker of the data set. See https://en.wikipedia.org/wiki/Box_plot for details.
    Input: percentile [25, 50 ,75]
    """
    return qt[2] + 1.5*(qt[2]-qt[0])

def apply_scale(df, U):
    """Apply the scaling function to the data set.
    """
    dt_tpm_lg = df['TPM'].apply(scalingfunc, U=U)
    df_lg = pd.DataFrame({ 
        'TPM (lg, U=%g)'%(U) : dt_tpm_lg,
        })
    return df_lg

def getGeneIdSetFromModel(modelfile):
    """Retrieve set of gene Ids from the modelfile.
       Return None if failed to parse the model file.
    """
    geneIdList = []
    with open(modelfile, "r") as stream:
        try:
            modeldata = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print("Failed to parse model file %s"%(modelfile),  file = sys.stderr)
            return None
    for (key, content) in modeldata:
        if key == "genes":
            for idlist in content:
                for idtuple in idlist:
                    geneIdList.append(idtuple[1])
    geneIdSet = set(geneIdList)
    return geneIdSet

def plot(df_lg, U, figfile):
    """Plot the curve of the scaling function and the distribution of the data set after scaling.
    """
    fig = plt.figure(figsize=(15,8))
    gs = fig.add_gridspec(1,2)

    ax1 = fig.add_subplot(gs[0,0])
    vf = np.vectorize(scalingfunc)
    x = np.linspace(0,100,2000)
    y = vf(x, U)
    ax1.set_xlabel('Raw data')
    ax1.set_ylabel('Scaled')
    ax1.plot(x, y, 'r')
    ax1.set_ylim(0, 1)

    ax2 = fig.add_subplot(gs[0,1])
    ax2.autoscale()
    ax2.set_xlabel('Scaled')
    ax2.set_ylabel('Percentage')
    df_lg.hist(column = df_lg.columns[0], bins = 50, ax=ax2, density=1)

    bsname = figfile
    fig.savefig('%s.png'%(bsname))
    fig.savefig('%s.pdf'%(bsname))
    plt.clf()
    plt.close()

def formatHpaRna(df, U, geneIdSet, outfile):
    """Output the scaled TPM values from the HPA RNA tissue data.
    Genes that do not exist in the Human-GEM model are filtered.
    """
    try:
        fpout = open(outfile, 'w')
        tissueList = df['Tissue'].unique().tolist()
        fpout.write("\t".join(["id"]+tissueList)+"\n")
        grouped_df = df.groupby(['Gene'])
        for key, item in grouped_df:
            if geneIdSet is None or key in geneIdSet:
                df1 = grouped_df.get_group(key).reset_index()
                values = df1['TPM'].apply(scalingfunc, U)
                dt = {}
                for i in range(len(df1)):
                    dt[df1.iloc[i]['Tissue']] = values[i]
                li = []
                li.append(key)
                for tissue in tissueList:
                    li.append("%.3g"%(dt[tissue]))
                fpout.write("\t".join(li)+"\n")
        fpout.close()
    except:
        sys.stderr.write("Failed to write to file %s"%(outfile))
        return 1

def outputDistFig(df, U, figfile):
    df_lg = apply_scale(df, U)
    plot(df_lg, U, figfile)


def main():
    parser = argparse.ArgumentParser(
            description='Format hpaRna TPM data to that is suitable for MetAtlas dataOverlay',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='''\
Example:
  format_hpaRna.py -i rna_tissue_hpa.tsv -o hpaRna.tsv -m Human-GEM_v1.8.0.yml -ofig ana_distribution

Created 2021-09-24, updated 2021-09-28, Nanjiang Shu
''')
    parser.add_argument('-i', metavar='FILE', dest='hpaRnaFile', required=True,
            help='provide the hpaRna data file')
    parser.add_argument('-m', metavar='FILE', dest='modelfile', required=True,
            help='provide the Human-GEM yaml file')
    parser.add_argument('-o' , metavar='OUTFILE', dest='outfile', required=True,
            help='output the result to outfile')
    parser.add_argument('-ofig' , metavar='OUTFILE', dest='figfile', required=False,
            help='output the distribution figure')

    args = parser.parse_args()
    hpaRnaFile = args.hpaRnaFile
    modelfile = args.modelfile
    outfile = args.outfile
    figfile = args.figfile

    geneIdSet = getGeneIdSetFromModel(modelfile)


    df = pd.read_csv(hpaRnaFile, sep='\t')

    # calculate parameter U
    qt_TPM = np.percentile(df['TPM'], [25, 50, 75])
    vf_log = np.vectorize(logfunc)
    U = logfunc(calMax(qt_TPM))

    formatHpaRna(df, U, geneIdSet, outfile)

    if figfile != None and figfile != "":
        outputDistFig(df, U, figfile)

    return 0

if __name__ == '__main__' :
    main()
