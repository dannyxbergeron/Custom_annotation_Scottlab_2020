import pandas as pd
import numpy as np

from pybedtools import BedTool as bt

ensembl_bed_file = snakemake.input.ensembl_bed
refseq_bed_file = snakemake.input.refseq_bed
snodb_file = snakemake.input.snodb
blockbuster_file = snakemake.input.blockbuster_bed
to_remove_file = snakemake.output.snoRNA_to_delete
out_file = snakemake.output.sno_gtf


def load_df(file):
    df = pd.read_csv(file, sep='\t')
    df = df.loc[(df.gene_biotype == 'snoRNA')
                | (df.gene_biotype == 'scaRNA')
                | (df.gene_biotype == 'guide_RNA')]
    id_refseqId = {}
    if 'Refseq' in file:
        id_refseqId = dict(zip(df.gene_id, df.refseq_id))

    df = df[['chr', 'start', 'end', 'gene_id', 'gene_name', 'strand']]
    print('initial num of snoRNA: ', len(df))
    return df, id_refseqId

def load_snodb(file):

    chromosomes = [f'chr{x}' for x in range(1,23)] + ['chrX', 'chrY']

    df = pd.read_csv(file, sep='\t')
    df.chr = 'chr' + df.chr
    df = df.loc[df.chr.isin(chromosomes)]
    df.rename(columns={'st': 'start', 'nd': 'end', 'symbol': 'gene_name', 'id': 'gene_id'}, inplace=True)
    df.gene_id = [f'snoDB{x}' for x in df.gene_id]
    df_bed = df[['chr', 'start', 'end', 'gene_id', 'gene_name', 'strand']]
    return df, df_bed

def load_blockbuster(file):
    df = pd.read_csv(file, sep='\t')
    df.chr = 'chr' + df.chr
    return df

def bedtools_intersect(df1, df2):
    df1_bt = bt.from_dataframe(df1)
    df2_bt = bt.from_dataframe(df2)
    intersect = df1_bt.intersect(df2_bt, wo=True, s=True, sorted=False)
    new_cols = ['chr1', 'start1', 'end1', 'gene_id1', 'gene_name1', 'strand',
                'chr2', 'start2', 'end2', 'gene_id2', 'gene_name2', 'strand2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    intersect_df = intersect_df[['start1', 'end1', 'start2', 'end2', 'gene_id1',
                                'gene_id2', 'gene_name1', 'gene_name2', 'overlap']]
    # print(intersect_df[['chr1', 'start1', 'end1', 'gene_id1', 'gene_biotype2']])
    return intersect_df

def remove_duplicates(df_):

    df = df_.copy(deep=True)
    intersect_df = bedtools_intersect(df, df)
    intersect_df = intersect_df.loc[intersect_df.gene_id1 != intersect_df.gene_id2]

    intersect_df['merge_id'] = [
        '_'.join(sorted([x,y]))
        for (x,y) in zip(intersect_df.gene_id1, intersect_df.gene_id2)
    ]
    intersect_df.drop_duplicates('merge_id', inplace=True)
    intersect_df.drop(columns=['merge_id'], inplace=True)
    intersect_df['length_1'] = intersect_df.end1 - intersect_df.start1
    intersect_df['length_2'] = intersect_df.end2 - intersect_df.start2

    to_remove = []
    idx_to_remove = -1
    for i in intersect_df.index:
        gene_name1 = intersect_df.at[i, 'gene_name1']
        gene_name2 = intersect_df.at[i, 'gene_name2']
        length_1 = intersect_df.at[i, 'length_1']
        length_2 = intersect_df.at[i, 'length_2']
        sno_bool1 = 'SNOR' in gene_name1 or gene_name1.startswith('U')
        sno_bool2 = 'SNOR' in gene_name2 or gene_name2.startswith('U')

        if (gene_name1 == gene_name2) or \
            (sno_bool1 and sno_bool2) or \
            (not sno_bool1 and not sno_bool2):
            idx_to_remove = 1 if length_1 < length_2 else 2
        elif sno_bool1:
            idx_to_remove = 2
        else:
            idx_to_remove = 1

        to_remove.append(intersect_df.at[i, f'gene_id{idx_to_remove}'])

    df = df.loc[~(df.gene_id.isin(to_remove))]
    print('num of snoRNA after removing the duplicates: ', len(df))
    return df, to_remove

def write_to_remove(to_remove):

    # Write the snoRNA to remove from the base ensembl annotation
    with open(to_remove_file, 'w') as f:
        for id in to_remove:
            f.write(id)
            f.write('\n')


def get_refseq_snoRNA(ensembl_df, refseq_df_, id_refseqId):

    refseq_df = refseq_df_.copy(deep=True)

    refseq_initial_ids = set(refseq_df.gene_id)
    intersect_df = bedtools_intersect(ensembl_df, refseq_df)
    overlapping = set(intersect_df.gene_id2)

    missing = refseq_initial_ids - overlapping
    print('Number of refseq missing snoRNA: ', len(missing))
    sno_refseq_missing_df = refseq_df.loc[refseq_df.gene_id.isin(missing)]
    sno_refseq_missing_df['gene_id'] = sno_refseq_missing_df['gene_id'].map(id_refseqId)

    return sno_refseq_missing_df

def get_snodb_missing(refseq_ensembl_snoRNAs_df, snodb_bed_df_):

    snodb_bed_df = snodb_bed_df_.copy(deep=True)

    initial_ids = set(snodb_bed_df.gene_id)
    intersect_df = bedtools_intersect(refseq_ensembl_snoRNAs_df, snodb_bed_df)
    overlapping = set(intersect_df.gene_id2)
    missing = initial_ids - overlapping

    snodb_bed_df = snodb_bed_df.loc[snodb_bed_df.gene_id.isin(missing)]
    snodb_bed_df = snodb_bed_df.loc[snodb_bed_df.gene_name != 'TERC']
    # print(snodb_bed_df)
    # print(len(missing))

    return snodb_bed_df


def create_snoRNA_gtf(df, out_file):

    score, frame = '.', '.'
    gene_version = 'gene_version "1";'
    transcript_version = 'transcript_version "1";'
    exon_version = 'exon_version "1";'
    exon_number = 'exon_number "1";'
    tag = 'tag "basic";'
    tsl = 'transcript_support_level "NA";'

    master_list = []
    for chr, start, end, gene_id, gene_name, strand, source, biotype in df.values:

        id = f'gene_id "{gene_id}";'
        transcript_id = f'transcript_id "{gene_id}";'
        exon_id = f'exon_id "{gene_id}";'
        name = f'gene_name "{gene_name}";'
        transcript_name = f'transcript_name "{gene_name}";'
        gene_source = f'gene_source "{source}";'
        transcript_source = f'transcript_source "{source}";'
        gene_biotype = f'gene_biotype "{biotype}";'
        transcript_biotype = f'transcript_biotype "{biotype}";'

        # For the gene
        att_list = [id, gene_version, name, gene_source, gene_biotype]
        attributes = ' '.join(att_list)
        merge = [chr, source, 'gene', start, end, score, strand, frame, attributes]
        master_list.append(merge)

        # For the transcript
        att_list = [id, gene_version, transcript_id, transcript_version,
                    name, gene_source, gene_biotype,
                    transcript_name, transcript_source, transcript_biotype, tag, tsl]
        attributes = ' '.join(att_list)
        merge = [chr, source, 'transcript', start, end, score, strand, frame, attributes]
        master_list.append(merge)

        # For the exon
        att_list = [id, gene_version, transcript_id, transcript_version,
                    exon_number,
                    name, gene_source, gene_biotype,
                    transcript_name, transcript_source, transcript_biotype,
                    exon_id, exon_version,
                    tag, tsl]
        attributes = ' '.join(att_list)
        merge = [chr, source, 'exon', start, end, score, strand, frame, attributes]
        master_list.append(merge)

    # For some reasons, the to_csv was adding nonsense " in the attribute field...
    with open(out_file, 'w') as f:
        for row in master_list:
            for val in row:
                val = str(val)
                f.write(val)
                f.write('\t')
            f.write('\n')


def main():

    snodb_df, snodb_bed_df = load_snodb(snodb_file)
    blockbuster_df = load_blockbuster(blockbuster_file)

    print('Ensembl----------------------------------------')
    ensembl_df, _ = load_df(ensembl_bed_file)
    ensembl_df, to_remove = remove_duplicates(ensembl_df)
    write_to_remove(to_remove)

    print('Refseq----------------------------------------')
    refseq_df, id_refseqId = load_df(refseq_bed_file)
    refseq_df, _  = remove_duplicates(refseq_df)

    # Get the refseq missing snoRNA
    refseq_missing_df = get_refseq_snoRNA(ensembl_df, refseq_df, id_refseqId)

    # Get the snoDB missing snoRNA
    blockbuster_df.columns = ensembl_df.columns # Just to be able to merge the dfs
    refseq_ensembl_snoRNAs_df = pd.concat([ensembl_df, refseq_missing_df, blockbuster_df])
    snodb_missing = get_snodb_missing(refseq_ensembl_snoRNAs_df, snodb_bed_df)

    # Remove the duplicates from the snodb missing
    snodb_missing, _ = remove_duplicates(snodb_missing)

    refseq_missing_df['source'] = 'BestRefSeq'
    snodb_missing['source'] = 'snoDB'

    full_missing = pd.concat([refseq_missing_df, snodb_missing])
    full_missing.chr = full_missing.chr.str.replace('chr', '')
    full_missing['biotype'] = np.where(full_missing.gene_name.str.startswith('SCA'), 'scaRNA', 'snoRNA')

    create_snoRNA_gtf(full_missing, out_file)




main()
