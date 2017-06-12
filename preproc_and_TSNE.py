import mypy
import numpy as np
import pandas as pd
from collections import OrderedDict
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import ipdb
import random
import math
import numpy as np
import sys, getopt
import argparse
import warnings


warnings.filterwarnings("ignore")

###################################
# Parse command line arguments
###################################
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description="Outputs the TSNE coordinates for a set of input cells. Input is the methylated and total base calls for each cell across bins/genes in tsv format. "+
                                 "DETAILS: Load the data, filter out any samples not in the QC list, impute mCH at bins/gene where coverage is poor, normalize the data, "+
                                 "run PCA (50 components) and TSNE (2 components) and save the output to file.")

parser.add_argument("-i", "--input", help="input file name. Must have columns named as sample1_mc, sample1_c, sample2_mc, etc.", 
                    required=True)
parser.add_argument("-o", "--output", help="output file name", required=True)
parser.add_argument("-s", "--species", help="mouse or human", default="mouse")
parser.add_argument("-n", "--normalize", help="normalize the data before running PCA and TSNE", default=True)
parser.add_argument("-p", "--perplexity", type=int, help="TSNE perplexity", default=25)
parser.add_argument("-d", "--seed", type=int, help="TSNE seed", default=1)
parser.add_argument("-b", "--base_call_cutoff", type=int, help="minimum base calls for a bin to not be imputed.", default=100)
parser.add_argument("-m", "--mdata", help="path to metdata file", required=True)
args = parser.parse_args()

species = args.species
normalize = args.normalize
perplexity = args.perplexity
infile = args.input
outfile = args.output
base_call_cutoff = args.base_call_cutoff
seed = args.seed
mdata = args.mdata




##################################
# Load data
##################################
print("Loading data.")


# Load input and metadata
df_gene_level_mCH = pd.read_csv(infile, sep="\t")
metadata = pd.read_csv(mdata, sep="\t")

df = df_gene_level_mCH
samples = df.samples.tolist()

print("Computing mCH levels.")
df = df.loc[(df.filter(regex='_c') > base_call_cutoff).sum(axis=1) >= .995*len(samples)]
df_mc = df[[x+'_mc' for x in samples]]
df_c = df[[x+'_c' for x in samples]]
df_c[df_c < base_call_cutoff] = np.nan
df_c.columns = [x+'_mcc' for x in samples]
df_mc.columns = [x+'_mcc' for x in samples]
df = df_mc/df_c

print('Imputing data.')
df = df.loc[df.count(axis=1) > 0]
df.reset_index(inplace=True, drop=True)
means = df.mean(axis=1)
fill_value = pd.DataFrame({col: means for col in df.columns})
df.fillna(fill_value, inplace=True)

if normalize:
    print('Normalizing.')
    for i,row in metadata.iterrows():
        samp = row['Sample']
        df[samp+'_mcc'] = (df[samp+'_mcc'] / (row['mCH/CH']+.01))



#######################
# RUNNING PCA
#######################
print("Running PCA.")
pca = PCA(n_components=50)
sklearn_transf = pca.fit_transform(df.T)
# print(pca.explained_variance_ratio_)
sklearn_transf_PCA = sklearn_transf



#######################
# RUNNING TSNE
#######################
print("Running TSNE.")
num_components = 2
tsne = TSNE(n_components=num_components, init='pca', random_state=seed, perplexity=perplexity, verbose=3)
sklearn_transf = tsne.fit_transform(sklearn_transf_PCA)

print("Saving output to file.")
df_tsne = pd.DataFrame(sklearn_transf, columns=['tsne_x','tsne_y'])
df_tsne['cells'] = [sample[:-4] for sample in df.columns.tolist()]
df_tsne.to_csv(outfile, sep="\t", index=False)


# plt.plot(sklearn_transf[:,0], sklearn_transf[:,1], '.r')

