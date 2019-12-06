import os
import multiprocessing as mp
import numpy as np



################
#  CREATE THE BINNED DATA FILES FROM ALLC FILES
################

######### INPUT PARAMETERS #################
species = 'mouse' # or human
samples = [] # LIST OF SAMPLE NAMES
allc_dir = '' # DIRECTORY WITH ALLC FILES IN IT
# FILE NAME TO ALLC FILES WILL BE CONSTRUCTED AS: allc_dir + "/allc_" + sample[n] + "_" + chromosome + ".tsv"
# OUTPUT FILE NAME: "binc_" + sample[n] + "_" + str(bin_size) + "_" + chromosome + ".tsv",



def read_allc(fname, position_as_index=True, compressed=False):
    # Read methylation data stored in an "allc" file
    if compressed:
        os.system("bgzip -cd " + fname + ".gz > " + fname)

    if position_as_index == True:
        df = pd.read_csv(fname, sep="\t", index_col=1, skiprows=1,
                         names=['chr','pos','strand','context','mc','c','methylated'])
    else:
        df = pd.read_csv(fname, sep="\t", skiprows=1,
                         names=['chr','pos','strand','context','mc','c','methylated'])

    if compressed:
        os.remove(fname)
    return df


def get_mCH_contexts():
    contexts = []
    for base1 in ['A','C','T']:
        for base2 in ['A','C','G','T']:
            contexts.append('C' + base1 + base2)
    return contexts

def get_mCG_contexts():
    return ['CGA','CGC','CGG','CGT']

def get_mouse_chromosomes(include_x=True):
    chromosomes = [str(x) for x in range(1,20)]
    if include_x:
        chromosomes.append('X')
    return chromosomes

def get_human_chromosomes(include_x=True):
    chromosomes = [str(x) for x in range(1,23)]
    if include_x:
        chromosomes.append('X')
    return chromosomes


### FUNCTION FOR BINNING THE ALLC FILES
def bin_allc(sample, path='.', bin_sizes=[1000,5000,10000,25000,100000], chromosomes=None, outpath = '.', species='mouse', compressed=False):

    if chromosomes == None:
        if species == 'human':
            chromosomes = get_human_chromosomes()
        else:
            chromosomes = get_mouse_chromosomes()

    for chromosome in chromosomes:

        fname = path + "/allc_" + sample + "_" + chromosome + ".tsv"

        if compressed:
            os.system("bgzip -cd " + fname + ".gz > " + fname)

        if not os.path.isfile(fname):
            print("bin_allc: " + fname + " does not exist.")
            return

        df = read_allc(fname)
        if compressed:
            os.remove(fname)

        for bin_size in bin_sizes:

            if species == 'human':
                bins = np.arange(0, get_chrom_lengths_human()[chromosome], bin_size)
            else:
                bins = np.arange(0, get_chrom_lengths_mouse()[chromosome], bin_size)

            # mCG
            df_CG = df.loc[df.context.isin(get_mCG_contexts())]
            groups = df_CG.groupby(pd.cut(df_CG.index, bins))
            mCG = groups.sum().mc.fillna(0)
            CG = groups.sum().c.fillna(0)

            # mCH
            df_CH = df.loc[df.context.isin(get_mCH_contexts())]
            groups = df_CH.groupby(pd.cut(df_CH.index, bins))
            mCH = groups.sum().mc.fillna(0)
            CH = groups.sum().c.fillna(0)

            data = np.array([bins[:len(bins)-1], mCG.values, CG.values, mCH.values, CH.values]).astype(int)
            binned_allc = pd.DataFrame(data.transpose(), columns=['bin','mCG','CG','mCH','CH'])
            binned_allc['chr'] = chromosome
            binned_allc = binned_allc[['chr','bin','mCG','CG','mCH','CH']]

            binned_allc.to_csv(outpath + "/binc_" + sample + "_" + str(bin_size) + "_" + chromosome + ".tsv",
                               sep="\t", header=False, index=False)


### PROCESS THE WILL CALL THE BINNING FUNCTION IN A PARALLEL FASHION
def process(samples):
    for sample in samples:
        print(sample)
        if not os.path.exists(sample):
            os.makedirs(sample)
        mypy.bin_allc(sample, path=allc_dir+sample, outpath=sample, species=species, bin_sizes=[100000], compressed=False)
    return True


procs = min(16, len(samples))

p = mp.Pool(processes=procs)
split_samples = np.array_split(samples,procs)
pool_results = p.map(process, split_samples)
p.close()
p.join()

