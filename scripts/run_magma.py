#!/usr/bin/env python

#import libraries
import subprocess

#generate variables for input files
magma_input = "/projectnb/bs859/students/amurrill/asd_proj/common-vs-rare-variants-ASD/data/processed/common_MAGMA.txt" #SNP file and pval file
geneloc_file = "/projectnb/bs859/students/amurrill/asd_proj/common-vs-rare-variants-ASD/ref/NCBI37.3.gene.loc"
anno_prefix = "/projectnb/bs859/students/amurrill/asd_proj/common-vs-rare-variants-ASD/data/processed/common_MAGMA"
gene_prefix = "/projectnb/bs859/students/amurrill/asd_proj/common-vs-rare-variants-ASD/data/processed/common_MAGMA_genes"
ref_panel = "/projectnb/bs859/data/1000G/plinkformat/1000G_EUR"

#perform annotation
anno_result = subprocess.run(
    ["magma", 
    "--annotate",
    "--snp-loc", magma_input,
    "--gene-loc", geneloc_file,
    "--out", anno_prefix],
    check = True
)
#do analysis
magma_result = subprocess.run(
    ["magma", 
    "--bfile", ref_panel,
    "--pval", magma_input, "ncol=NOBS", #change N to newly calculated NOBS
    "--gene-annot", anno_prefix + ".genes.annot",
    "--out", gene_prefix],
    check = True
)

