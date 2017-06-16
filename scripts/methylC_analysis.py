import os
import sys
import pandas as pd
import pybedtools as pbt

functionsfile = '../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

allc="allc/Fvesca_"+sys.argv[1]+"_allc_total.tsv"
fasta="../ref/"+sys.argv[1]+"/bowtie2/Fvesca_"+sys.argv[1]+".fa"
genome_file="../ref/"+sys.argv[1]+"/misc/Fvesca_"+sys.argv[1]+".genome"
genes_gff="../ref/"+sys.argv[1]+"/misc/Fvesca_"+sys.argv[1]+".gff"
filter_chr=['ChrL','ChrC','ChrM']
context=['CG','CHG','CHH','CAA','CAT','CAC','CAG','CTA','CTT','CTC',
         'CTG','CCA','CCT','CCC','CCG','CGA','CGT','CGC','CGG']

print("Finding mC distributions")
functions.mC_distribution(genome_file,allc,genes_gff,fasta,
                          windows=100000,stepsize=50000,cutoff=0,filter=0.5,feature='gene',
                          output_path="results/",save_windows=True,return_windows=False)

if os.path.exists(genes_gff):
    print("Gene metaplot")
    functions.gene_metaplot(allc,genes_gff,genome_file,output="results/gene_metaplot.tsv",
                            ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0,
                            first_feature='mRNA',second_feature='CDS',filter_chr=filter_chr,
                            remove_tmp=False)
else:
    print("No gene annotations found")
