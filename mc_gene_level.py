#!/usr/bin/python3
#
# Read in allc tables and generate summarized mc levels within gene bodies
#

import pandas as pd
# import mypy
import tabix
import multiprocessing as mp
import numpy as np
import os
import snmcseq_utils as mypy
import argparse

def mc_gene_level(samples,
    genebody='/cndd/projects/Public_Datasets/references/hg19/transcriptome/gencode.v19.annotation_genes_mypy.tsv',
    outdir="./genebody"):
    """
    samples = list of file names of allc files
    genebody = BED file with gene body or other annotation features
    """

    chromosomes = mypy.get_human_chromosomes()

    df_gtf = pd.read_csv(genebody, sep="\t")
    df_gtf.chr = df_gtf.chr.apply(lambda x: x[3:])
    df_gtf = df_gtf.loc[df_gtf.chr.isin(chromosomes)]

    print("mc_gene_level processing: "+sample)
    
    outfilename = os.path.join(outdir,sample+"_mCH_genebody.txt")
    if os.path.isfile(outfilename):
        print "File exists "+outfilename+", skipping..."
        continue

    outfile_CH = open(, "w")
    outfile_CH.write("id\tname\tchr\tstart\tend\tstrand\tmc\tc\n")

    for i,row in df_gtf.iterrows():

        allc = tabix.open(sample+'_bismark/allc_'+sample+'_'+row['chr']+'.tsv.gz')

        # Gene body CH
        records = allc.query(row['chr'], row['start'], row['end'])
        mc, c = mypy.tabix_summary(records, context="CH", cap=2)
        outfile_CH.write(row['gene_id'] + "\t" + row['name'] + "\t" + row['chr'] + "\t" + str(row['start']) + "\t" + 
           str(row['end']) + "\t" + row['strand'] + "\t" + str(mc) + "\t" + str(c) + "\n")

        return True
        
        procs = min(len(samples), 16)
        p = mp.Pool(processes=procs)
        split_samples = np.array_split(samples,procs)
        pool_results = p.map(process, split_samples)
        p.close()
        p.join()

def mc_gene_level(samples,
    genebody='/cndd/projects/Public_Datasets/references/hg19/transcriptome/gencode.v19.annotation_genes_mypy.tsv',
    outdir="./genebody"):
