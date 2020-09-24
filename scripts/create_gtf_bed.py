import numpy as np
import pandas as pd

parsed_gtf = snakemake.input.parsed_gtf
out_file = snakemake.output.gene_bed

df = pd.read_csv(parsed_gtf, sep='\t', dtype={'seqname': str})


def transform_refseq(df_):
    """ Deal with the Refseq nonsense chromosome notation and the ids that
        are only in the exon... """
    df = df_.copy(deep=True)
    df = df.loc[df.feature.isin(['gene'])]
    # keep the exon for the gene_id of the snoRNA
    df.rename(columns={'gene': 'gene_name'}, inplace=True)

    # Change the chr naming
    seqname = []
    for chr in df['seqname']:
        chr = chr.split('.')[0].replace('NC_00000', 'chr').replace('NC_0000', 'chr')
        if chr == 'chr23': chr = 'chrX'
        elif chr == 'chr24': chr = 'chrY'
        seqname.append(chr)
    df['seqname'] = seqname
    df = df.loc[df.seqname.str.startswith('chr')]

    # Add the real gene_id to each gene because IT'S NOT THERE !!!!!
    exon_df = df_.loc[df_.feature == 'exon']
    gene_id_dict = {
        exon_df.at[i, 'gene_id']: exon_df.at[i, 'transcript_id'].split('.')[0]
        for i in exon_df.index
    }
    df['refseq_id'] = df['gene_id'].map(gene_id_dict)
    df['refseq_id'] = np.where(df.gene_id.isna(), df['gene_id'], df['refseq_id'])
    df = df[['seqname', 'start', 'end', 'gene_id', 'gene_name', 'strand',
             'gene_biotype', 'refseq_id']]
    return df


def transform_ensembl(df_):
    """ Simply add chr in front of chromosome for further manipulation """
    df = df_.copy(deep=True)
    df = df.loc[df.feature.isin(['gene'])]
    df['seqname'] = 'chr' + df['seqname']
    df = df[['seqname', 'start', 'end', 'gene_id', 'gene_name', 'strand',
             'gene_biotype']]
    return df


if 'Refseq' in parsed_gtf:
    df = transform_refseq(df)
else:
    df = transform_ensembl(df)


df.rename(columns={'seqname': 'chr'}, inplace=True)

df.to_csv(out_file, sep='\t', index=False)
