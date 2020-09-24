import pandas as pd
import numpy as np

parsed_gtf = snakemake.input.gene_only
start_snodb = snakemake.input.start_snodb
out_file = snakemake.output.associated_snoDB

def load_gtf(file):

    df = pd.read_csv(file, sep='\t', dtype=str)
    df = df.loc[df.gene_biotype.isin(['snoRNA', 'scaRNA'])]
    df.rename(columns={'seqname': 'chr'}, inplace=True)
    df = df[['chr', 'start', 'end', 'gene_name', 'gene_id', 'strand']]
    df['merged_id'] = df['chr'] + '_' + df['start'] + '_' + df['end']
    return df

def load_snodb(file):
    df = pd.read_csv(file, sep='\t', dtype=str)
    df['merged_id'] = df['chr'] + '_' + df['st'] + '_' + df['nd']
    return df

def get_gene_id(snodb_df_, gtf_df):

    mId_geneId_dict = dict(zip(gtf_df.merged_id, gtf_df.gene_id))
    mId_geneNamedict = dict(zip(gtf_df.merged_id, gtf_df.gene_name))

    snodb_df = snodb_df_.copy(deep=True)
    snodb_df['gene_id_annot2020'] = snodb_df['merged_id'].map(mId_geneId_dict)
    snodb_df['gene_name_annot2020'] = snodb_df['merged_id'].map(mId_geneNamedict)
    snodb_df.drop(columns='merged_id', inplace=True)

    snodb_df = snodb_df[['id', 'gene_id_annot2020', 'gene_name_annot2020',
                            *list(snodb_df.columns)[1:-2]]]

    snodb_df.to_csv(out_file, sep='\t', index=False)

def main():

    gtf_df = load_gtf(parsed_gtf)
    snodb_df = load_snodb(start_snodb)

    final_df = get_gene_id(snodb_df, gtf_df)


main()
