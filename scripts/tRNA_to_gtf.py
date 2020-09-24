import pandas as pd

tRNA_file = snakemake.input.tRNA_file
out_file = snakemake.output.tRNA_gtf


def load_df(tRNA_file):
    colnames = ['chr', 'start', 'end', 'gene_id', 'scorecrap',
                'strand', 'st', 'nd', 'crap1', 'crap2', 'crap3', 'crap4']
    df = pd.read_csv(tRNA_file, sep='\t', names=colnames)
    df = df[['chr', 'start', 'end', 'gene_id', 'strand']]
    df['chr'] = df['chr'].str.replace('chr', '')
    return df


def create_tRNA_gtf(df):
    # attributes
    gene_version = 'gene_version "1";'
    gene_source = 'gene_source "gtrnadb";'
    gene_biotype = 'gene_biotype "tRNA";'
    transcript_version = 'transcript_version "1";'
    transcript_source = 'transcript_source "gtrnadb";'
    transcript_biotype = 'transcript_biotype "tRNA";'
    exon_number = 'exon_number "1";'
    exon_version = 'exon_version "1";'
    tag = 'tag "basic";'
    tsl = 'transcript_support_level "NA";'

    master_list = []
    for chr, start, end, id, strand in df.values:
        # for gene
        gene_id = f'gene_id "{id}";'
        gene_name = f'gene_name "{id}";'
        gene_attributes = ' '. join([gene_id, gene_version, gene_name, gene_source, gene_biotype])
        merge = [chr, 'gtrnadb', 'gene', start, end, '.', strand, '.', gene_attributes]
        master_list.append(merge)

        # for transcript
        transcript_id = f'transcript_id "{id}";'
        transcript_name = f'transcript_name "{id}";'
        att_list = [gene_id, gene_version, transcript_id, transcript_version,
                    gene_name, gene_source, gene_biotype, transcript_name,
                    transcript_source, transcript_biotype, tag, tsl]
        gene_attributes = ' '. join(att_list)
        merge = [chr, 'gtrnadb', 'transcript', start, end, '.', strand, '.', gene_attributes]
        master_list.append(merge)

        # for exon
        exon_id = f'exon_id "{id}";'
        att_list = [gene_id, gene_version, transcript_id, transcript_version,
                    exon_number, gene_name, gene_source, gene_biotype, transcript_name,
                    transcript_source, transcript_biotype, exon_id, exon_version, tag, tsl]
        gene_attributes = ' '. join(att_list)
        merge = [chr, 'gtrnadb', 'exon', start, end, '.', strand, '.', gene_attributes]
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

    tRNA_df = load_df(tRNA_file)

    create_tRNA_gtf(tRNA_df)


main()
