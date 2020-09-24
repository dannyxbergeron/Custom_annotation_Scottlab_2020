from snakemake.shell import shell

import pandas as pd
import numpy as np

from pybedtools import BedTool as bt

input_file = snakemake.input.gene_only
out_file = snakemake.output.tok

def load_df(file):

    biotypes = ['tRNA', 'snoRNA', 'rRNA', 'scaRNA', 'intergenic_cluster',
                'intronic_cluster', 'pre-tRNA', 'tRNA_fragment']

    df = pd.read_csv(file, sep='\t')
    df = df.loc[df.gene_biotype.isin(biotypes)]
    df.rename(columns={'seqname': 'chr'}, inplace=True)
    df = df[['chr', 'start', 'end', 'gene_name', 'gene_biotype', 'strand']]
    return df


def bedtools_intersect(df1, df2):
    df1_bt = bt.from_dataframe(df1)
    df2_bt = bt.from_dataframe(df2)
    intersect = df1_bt.intersect(df2_bt, wo=True, s=True, sorted=False)
    new_cols = ['chr1', 'start1', 'end1', 'gene_name1', 'gene_biotype1', 'strand',
                'chr2', 'start2', 'end2', 'gene_name2', 'gene_biotype2', 'strand2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    intersect_df = intersect_df.loc[intersect_df.gene_name1 != intersect_df.gene_name2]
    return intersect_df

def main():

    df = load_df(input_file)

    intersect_df = bedtools_intersect(df, df)
    print(len(intersect_df))
    print(intersect_df[['start1', 'start2', 'end1', 'end2', 'gene_name1',
                        'gene_name2', 'gene_biotype1', 'gene_biotype2']])

    shell(f'touch {out_file}')

main()
