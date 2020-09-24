import pandas as pd
import numpy as np

from pybedtools import BedTool as bt

rRNA_file = snakemake.input.refseq_parsed_file
refseq_bed = snakemake.input.refseq_bed
ens_bed = snakemake.input.ensembl_bed
out_file = snakemake.output.rRNA_gtf
to_remove_out = snakemake.output.rRNA_to_delete

all_refseq_df = pd.read_csv(rRNA_file, sep='\t')
all_refseq_df = all_refseq_df.loc[all_refseq_df.gene_biotype == 'rRNA']
# all_refseq_df = all_refseq_df[['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
#                                 'gene_biotype', 'gene_id', 'gene']]

def load_df(file):

    df = pd.read_csv(file, sep='\t')
    df = df.loc[df.gene_biotype == 'rRNA']
    df_bed = df[['chr', 'start', 'end', 'gene_id', 'gene_name', 'strand']]
    gene_set = set(df.gene_id)
    return df_bed, gene_set, df

def bedtools_intersect(df1, df2):
    df1_bt = bt.from_dataframe(df1)
    df2_bt = bt.from_dataframe(df2)
    intersect = df1_bt.intersect(df2_bt, wo=True, s=True, sorted=False)
    new_cols = ['chr1', 'start1', 'end1', 'gene_id1', 'gene_name1', 'strand',
                'chr2', 'start2', 'end2', 'gene_id2', 'gene_name2', 'strand2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    return intersect_df

def get_missing(intersect_df, ref_gene_set, ensembl_bed_df, ref_bed_df):

    # The missing genes are RNA18S, RNA28S and 45S
    missing = ref_gene_set - set(intersect_df.gene_id2)
    refseq_missing_df = ref_bed_df.loc[(ref_bed_df.gene_id.isin(missing))
                                | (ref_bed_df.gene_id.str.contains('45S'))]
    # Remove the ensembl rRNA overlapping with the 45S
    ensembl_45S_overlap_set = set(intersect_df.loc[intersect_df.gene_id2.str.contains('45S')].gene_id1)
    ensembl_45S_overlap_df = ensembl_bed_df.loc[ensembl_bed_df.gene_id.isin(ensembl_45S_overlap_set)]
    to_remove = ensembl_45S_overlap_df[['gene_id']]
    # print(to_remove)
    to_remove.to_csv(to_remove_out, sep='\t', index=False, header=False)

    # Get the layers to pur the 18S, 5.8S and 28S in exons
    gene45S_df = ref_bed_df.loc[ref_bed_df.gene_id.str.contains('45S')]
    gene_exon_df = ref_bed_df.loc[~(ref_bed_df.gene_id.str.contains('45S'))]
    all_missing_df = bedtools_intersect(gene45S_df, gene_exon_df)
    return all_missing_df

def create_to_add_file(missing_df, refid_geneId_dict):

    source, score, frame = 'BestRefSeq', '.', '.'
    gene_version = 'gene_version "1";'
    transcript_version = 'transcript_version "1";'
    exon_version = 'exon_version "1";'
    gene_source = 'gene_source "BestRefSeq";'
    transcript_source = 'transcript_source "BestRefSeq";'
    gene_biotype = 'gene_biotype "rRNA";'
    transcript_biotype = 'transcript_biotype "rRNA";'
    tag = 'tag "basic";'
    tsl = 'transcript_support_level "NA";'

    master_list = []
    for long_rRNA in set(missing_df.gene_id1):
        tmp = missing_df.loc[missing_df.gene_id1 == long_rRNA]
        tmp.reset_index(drop=True, inplace=True)
        # print(tmp)
        matrix = tmp.values

        chr, start, end, gene_id, gene_name, strand = matrix[0,:6]
        chr = chr.replace('chr', '')
        refseq_id = refid_geneId_dict[gene_id]
        id = f'gene_id "{refseq_id}";'
        transcript_id = f'transcript_id "{refseq_id}";'
        name = f'gene_name "{gene_name}";'
        transcript_name = f'transcript_name "{gene_name}";'

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

        # For the exons (18S, 28S and 5.8S)
        exons = matrix[:,6:-1]
        count = 1
        for _, exon_start, exon_end, exon_id, exon_name, _ in exons:
            refseq_id_exon = refid_geneId_dict[exon_id]
            exon_id = f'exon_id "{refseq_id_exon}";'
            exon_number = f'exon_number "{count}";'

            att_list = [id, gene_version, transcript_id, transcript_version,
                        exon_number,
                        name, gene_source, gene_biotype,
                        transcript_name, transcript_source, transcript_biotype,
                        exon_id, exon_version,
                        tag, tsl]
            attributes = ' '.join(att_list)
            merge = [chr, source, 'exon', exon_start, exon_end, score, strand, frame, attributes]
            master_list.append(merge)
            count += 1

    # For some reasons, the to_csv was adding nonsense " in the attribute field...
    with open(out_file, 'w') as f:
        for row in master_list:
            for val in row:
                val = str(val)
                f.write(val)
                f.write('\t')
            f.write('\n')




def main():

    ref_bed_df, ref_gene_set, ref_original = load_df(refseq_bed)
    ensembl_bed_df, ens_gene_set, ensembl_original = load_df(ens_bed)

    intersect_df = bedtools_intersect(ensembl_bed_df, ref_bed_df)

    all_missing_df = get_missing(intersect_df, ref_gene_set,
                                 ensembl_bed_df, ref_bed_df)

    refid_geneId_dict = dict(zip(ref_original.gene_id, ref_original.refseq_id))
    create_to_add_file(all_missing_df, refid_geneId_dict)


main()
