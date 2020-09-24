from os.path import join

configfile: "config.json"

rule all:
    input:
        final_gtf = join(config['path']['final'],
                         config['gtf']['final'] + '.gtf'),
        gene_only = join(config['path']['final'],
                         config['database_file']['parsed_final_gene_only']),
        associated_snoDB = join(config['path']['final'],
                                config['database_file']['final_snodb']),
        # tok = join(config['path']['test'], 'overlap.tok')


rule download_all:
    """ Download all the data needed for the complete gtf file construction """
    output:
        start_ensembl = join(config['path']['ref'],
                             config['gtf']['ensembl'] + '.gtf'),
        start_refseq = join(config['path']['ref'],
                            config['gtf']['refseq'] + '.gtf'),
        start_snodb = join(config['path']['ref'],
                           config['database_file']['snodb'])
    params:
        link_ensembl = config['download']['ensembl'],
        link_refseq = config['download']['refseq'],
        link_snodb = config['download']['snodb']
    shell:
        "wget --quiet -O {output.start_ensembl}.gz {params.link_ensembl} && "
        "wget --quiet -O {output.start_refseq}.gz {params.link_refseq} && "
        "wget --quiet -O {output.start_snodb} {params.link_snodb} && "
        "gunzip data/references/*.gz"

rule parse_gtf:
    """ Parse the gtf files into tabular files """
    input:
        gtf = join(config['path']['ref'], "{gtf}.gtf")
    output:
        parsed_gtf = join(config['path']['parsed_gtf'], "{gtf}.tsv")
    shell:
        "echo {input.gtf} && "
        "scripts/gtfParser find_parse {input.gtf} > {output.parsed_gtf}"

rule create_beds:
    """ Create bed files with only genes to use with pybedtools """
    input:
        parsed_gtf = join(config['path']['parsed_gtf'], "{gtf}.tsv")
    output:
        gene_bed = join(config['path']['bed'], "{gtf}.bed")
    conda:
        "envs/python.yaml"
    script:
        "scripts/create_gtf_bed.py"

rule create_tRNA_annotation:
    """ Create the gtf part for the tRNAs """
    input:
        tRNA_file = join(config['path']['ref'],
                         config['database_file']['tRNAs']),
    output:
        tRNA_gtf = join(config['path']['to_add_gtf'],
                        config['database_file']['gtf_tRNAs'])
    conda:
        "envs/python.yaml"
    script:
        "scripts/tRNA_to_gtf.py"

rule create_rRNA_annotation:
    """ Create the gtf part for the rRNAs from RefSeq annotation """
    input:
        refseq_bed = join(config['path']['bed'],
                          config['gtf']['refseq'] + '.bed'),
        ensembl_bed = join(config['path']['bed'],
                           config['gtf']['ensembl'] + '.bed'),
        refseq_parsed_file = join(config['path']['parsed_gtf'],
                                  config['gtf']['refseq'] + '.tsv'),
    output:
        rRNA_gtf = join(config['path']['to_add_gtf'],
                        config['database_file']['gtf_rRNAs']),
        rRNA_to_delete = join(config['path']['to_remove_gtf'],
                        config['database_file']['rRNA_remove_gtf'])
    conda:
        "envs/python.yaml"
    script:
        "scripts/rRNA_to_gtf.py"

rule get_blockbuster_genes:
    """ Get the blockbuster genes from our old custom annotation and the table
        from the blockbuster article """
    input:
        current_custom = join(config['path']['ref'],
                                config['gtf']['current_custom'] + '.gtf'),
        parsed_custom = join(config['path']['parsed_gtf'],
                                config['gtf']['current_custom'] + '.tsv'),
        blockbuster_table = join(config['path']['ref'],
                                config['database_file']['blockbuster']),
    output:
        blockbuster_gtf = join(config['path']['to_add_gtf'],
                                config['database_file']['gtf_blockbuster']),
        blockbuster_bed = join(config['path']['bed'],
                                config['database_file']['blockbuster_bed'])
    params:
        tmp = 'data/tmp',
        fetch_script = 'scripts/fetch_gtf'
    conda:
        "envs/python.yaml"
    script:
        "scripts/blockbuster_to_gtf.py"

rule get_snoRNAs:
    """ Get the snoRNAs from RefSeq and snodb """
    input:
        ensembl_bed = join(config['path']['bed'],
                          config['gtf']['ensembl'] + '.bed'),
        refseq_bed = join(config['path']['bed'],
                          config['gtf']['refseq'] + '.bed'),
        snodb = join(config['path']['ref'],
                     config['database_file']['snodb']),
        blockbuster_bed = join(config['path']['bed'],
                               config['database_file']['blockbuster_bed'])
    output:
        sno_gtf = join(config['path']['to_add_gtf'],
                       config['database_file']['gtf_snoRNA']),
        snoRNA_to_delete = join(config['path']['to_remove_gtf'],
                                config['database_file']['snoRNA_remove_gtf'])
    conda:
        "envs/python.yaml"
    script:
        "scripts/snoRNA_to_gtf.py"

rule create_final_annotation:
    """ Remove the duplicates genes in the annotation and add the missing
        genes. It also adds the snoDB version in the header """
    input:
        start_ensembl = join(config['path']['ref'],
                             config['gtf']['ensembl'] + '.gtf'),
        rRNA_to_delete = join(config['path']['to_remove_gtf'],
                              config['database_file']['rRNA_remove_gtf']),
        snoRNA_to_delete = join(config['path']['to_remove_gtf'],
                                config['database_file']['snoRNA_remove_gtf']),
        tRNA_gtf = join(config['path']['to_add_gtf'],
                        config['database_file']['gtf_tRNAs']),
        rRNA_gtf = join(config['path']['to_add_gtf'],
                        config['database_file']['gtf_rRNAs']),
        blockbuster_gtf = join(config['path']['to_add_gtf'],
                               config['database_file']['gtf_blockbuster']),
        sno_gtf = join(config['path']['to_add_gtf'],
                       config['database_file']['gtf_snoRNA'])
    output:
        final_gtf = join(config['path']['final'],
                         config['gtf']['final'] + '.gtf')
    shell:
        "cat {input.rRNA_to_delete} {input.snoRNA_to_delete} > data/tmp && "
        "grep -v -f data/tmp {input.start_ensembl} > {output.final_gtf} && "
        "grep '^#!' {output.final_gtf} > data/tmp && "
        "echo '#!snoDB-version 1.1.1' >> data/tmp && "
        "grep -v '^#!' {output.final_gtf} >> data/tmp && "
        "mv data/tmp {output.final_gtf} && "
        "cat {input.tRNA_gtf} {input.rRNA_gtf} "
        "{input.blockbuster_gtf} {input.sno_gtf} >> {output.final_gtf}"

rule parse_final_gtf:
    """ Parse the final gtf to get a tabular full version and another one just
        for the genes (feature == 'gene') """
    input:
        final_gtf = join(config['path']['final'],
                         config['gtf']['final'] + '.gtf')
    output:
        final_parsed = join(config['path']['final'],
                         config['database_file']['parsed_final']),
        gene_only = join(config['path']['final'],
                         config['database_file']['parsed_final_gene_only'])
    shell:
        "scripts/gtfParser find_parse {input.final_gtf} > {output.final_parsed} && "
        "awk -F '\t' '$3 == \"gene\" || NR == 1' {output.final_parsed} > {output.gene_only} "

rule create_associated_snoDB:
    """ Add a gene_id column corresponding to the gene_id in the gtf """
    input:
        gene_only = join(config['path']['final'],
                         config['database_file']['parsed_final_gene_only']),
        start_snodb = join(config['path']['ref'],
                           config['database_file']['snodb'])
    output:
        associated_snoDB = join(config['path']['final'],
                                config['database_file']['final_snodb'])
    conda:
        "envs/python.yaml"
    script:
        "scripts/create_final_snoDB.py"


include: "rules/test.smk"
