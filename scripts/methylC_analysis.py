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
unique_genes="../ref/"+sys.argv[1]+"/misc/Fvesca_"+sys.argv[1]+"_unique.gff"
syntenic_genes="../ref/"+sys.argv[1]+"/misc/Fvesca_"+sys.argv[1]+"_syntelogs.gff"
repeats_gff="../ref/"+sys.argv[1]+"/misc/Fvesca_"+sys.argv[1]+"_filtered_repeats.gff"
filter_chr=['ChrL','ChrC','ChrM']
context=['CG','CHG','CHH','CAA','CAT','CAC','CAG','CTA','CTT','CTC',
         'CTG','CCA','CCT','CCC','CCG','CGA','CGT','CGC','CGG']

print("Finding Total Weighted Methylation")
functions.weighted_mC(allc,output="results/total_weighted_mC.tsv",cutoff=0)

print("Getting per-site methylation levels")
functions.per_site_mC(allc,output_path='results/',context=context)

print("Analyzing Subcontext Methylation")
functions.subcontext_methylation(allc,fasta,context=context,output='results/subcontext_methylation.tsv',filter_chr=filter_chr)

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

if os.path.exists(syntenic_genes):
    print("Syntenic genes metaplot")
    functions.feature_metaplot(allc,syntenic_genes,genome_file,output="results/syntenic_gene_metaplot.tsv",
                               ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0)
else:
    print("No syntentic genes found")

if os.path.exists(unique_genes):
    print("Unique genes metaplot")
    functions.feature_metaplot(allc,unique_genes,genome_file,output="results/unique_gene_metaplot.tsv",
                               ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0)
else:
    print("No unique genes found")

if os.path.exists(genes_gff):
    print("Gene methylation levels")
    functions.feature_mC_levels('CDS_allc.tmp',genes_gff,output="results/gene_methylation_levels.tsv",
                                cutoff=0,filter_features='gene',filter_chr=filter_chr)
else:
    print("No gene annotations found")
