# Performing RER test

1. Run **core_genes_phylo.nf**, which extracts SNVs from core genes, to
create alignment fasta files and uses FastTree to make phylogenetic
reconstruction.

2. Run **snvs_to_rerconverge.nf**. Performs RERconverge test
    * First runs **get_all_gene_alns.r**

3. **create_permuted_maps.r**
