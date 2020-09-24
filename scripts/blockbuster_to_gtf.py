from snakemake.shell import shell

import pandas as pd
import numpy as np

from pybedtools import BedTool as bt

current_custom_file = snakemake.input.current_custom
parsed_custom_file = snakemake.input.parsed_custom
blockbuster_file = snakemake.input.blockbuster_table
tmp = snakemake.params.tmp
fetch_script = snakemake.params.fetch_script
blockbuster_bed_file = snakemake.output.blockbuster_bed
out_file = snakemake.output.blockbuster_gtf


def load_parsed(file):

    df = pd.read_csv(file, sep='\t')
    df = df.loc[(df.source == 'blockbuster_SK1') & (df.feature == 'gene')]
    df = df[['seqname', 'start', 'end', 'gene_id', 'gene_biotype', 'strand']]
    df.sort_values(['seqname', 'start', 'end'], inplace=True)
    return df


def load_blockbuster(file):

    df = pd.read_csv(file, sep='\t')
    df = df[['chr', 'start', 'end', 'Cluster ID', 'Final Biotype', 'strand']]
    df.sort_values(['chr', 'start', 'end'], inplace=True)
    return df


def bedtools_intersect(df1, df2):
    df1_bt = bt.from_dataframe(df1)
    df2_bt = bt.from_dataframe(df2)
    intersect = df1_bt.intersect(df2_bt, wo=True, s=True, sorted=False)
    new_cols = ['chr1', 'start1', 'end1', 'gene_id1', 'gene_biotype1', 'strand',
                'chr2', 'start2', 'end2', 'gene_id2', 'gene_biotype2', 'strand2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    # print(intersect_df[['chr1', 'start1', 'end1', 'gene_id1', 'gene_biotype2']])
    return intersect_df

def create_blockbuster_gtf(df, tmp, custom_gtf_file, out_file):

    with open(tmp, 'w') as f:
        for id, biotype  in zip(df.gene_id1, df.gene_biotype2):
            if 'snoRNA' in biotype:
                biotype = 'snoRNA'
            elif biotype == 'tRNA fragment':
                biotype = 'tRNA_fragment'
            f.write(f'gene_id "{id}"\t{biotype}')
            f.write('\n')

    cmd = f'{fetch_script} {tmp} {custom_gtf_file} > {out_file} && rm {tmp}'
    shell(cmd)


def main():

    parsed_custom_df = load_parsed(parsed_custom_file)

    blockbuster_table = load_blockbuster(blockbuster_file)
    # File needed for the removing of overlaping snodb snoRNA
    blockbuster_table.to_csv(blockbuster_bed_file, sep='\t', index=False)

    intersect_df = bedtools_intersect(parsed_custom_df, blockbuster_table)

    create_blockbuster_gtf(intersect_df, tmp, current_custom_file, out_file)


main()
