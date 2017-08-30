# Library with utilities from Chris Keown's "mypy"

def get_human_chromosomes(include_x=True):
    chromosomes = [str(x) for x in range(1,23)]
    if include_x:
        chromosomes.append('X')
    return chromosomes


def read_gencode_human(version='v19', pc=False):
	# pc = protein coding
    prefix = '/cndd/projects/Public_Datasets/references/hg19/transcriptome/'
    if pc:
        fname= prefix+'gencode.'+version+'.annotation_genes_pc_mypy.tsv'
    else:
        fname= prefix+'gencode.'+version+'.annotation_genes_mypy.tsv'
    return pd.read_csv(fname, sep="\t")
