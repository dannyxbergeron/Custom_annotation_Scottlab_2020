rule check_overlapping_genes:
    input:
        gene_only = join(config['path']['final'],
                         config['database_file']['parsed_final_gene_only'])
    output:
        tok = join(config['path']['test'], 'overlap.tok')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/check_overlap_test.py"
