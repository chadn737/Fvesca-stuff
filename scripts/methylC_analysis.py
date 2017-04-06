import os
import sys
import pandas as pd
import pybedtools as pbt

functionsfile = '../../../../scripts/functions.py'
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

print("Gene metaplot")
mC_bed = functions.allc2bed(allc)
bed = pbt.BedTool(genes_gff).filter(functions.feat_filter,'gene').filter(functions.chr_filter,filter_chr)
flank_bed = pbt.bedtool.BedTool.flank(bed,g=genome_file,l=2000,r=2000,s=True).saveas('f_tmp')
cds_bed = bed.filter(functions.feat_filter,'CDS').filter(functions.chr_filter,filter_chr).saveas('c_tmp')
bed = cds_bed.cat(flank_bed, postmerge=False)
mapping = pbt.bedtool.BedTool.intersect(mC_bed,bed,wa=True)
m = pd.read_table(mapping.fn, header=None, usecols = [0,1,5,6,7,8,9])
m.columns = ['chr','pos','strand','mc_class','mc_count','total','methylated']
m = m.drop_duplicates()
m.to_csv('CDS_allc.tmp', sep='\t', index=False)

functions.feature_metaplot('CDS_allc.tmp',genes_gff,genome_file,output="results/gene_metaplot.tsv",
                           ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0,
                           filter_features='gene',filter_chr=filter_chr)

for i in ['CDS_allc.tmp','c_tmp','f_tmp']:
    os.remove(i)

if os.path.exists(syntenic_genes):
    print("Syntenic genes metaplot")
    mC_bed = functions.allc2bed(allc)
    bed = pbt.BedTool(unique_genes).filter(functions.feat_filter,'gene').filter(functions.chr_filter,filter_chr)
    flank_bed = pbt.bedtool.BedTool.flank(bed,g=genome_file,l=2000,r=2000,s=True).saveas('f_tmp')
    cds_bed = bed.filter(functions.feat_filter,'CDS').filter(functions.chr_filter,filter_chr).saveas('c_tmp')
    bed = cds_bed.cat(flank_bed, postmerge=False)
    mapping = pbt.bedtool.BedTool.intersect(mC_bed,bed,wa=True)
    m = pd.read_table(mapping.fn, header=None, usecols = [0,1,5,6,7,8,9])
    m.columns = ['chr','pos','strand','mc_class','mc_count','total','methylated']
    m = m.drop_duplicates()
    m.to_csv('CDS_allc.tmp', sep='\t', index=False)

    functions.feature_metaplot('CDS_allc.tmp',syntenic_genes,genome_file,output="results/syntenic_gene_metaplot.tsv",
                               ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0,
                               filter_features='gene',filter_chr=filter_chr)

    for i in ['CDS_allc.tmp','c_tmp','f_tmp']:
        os.remove(i)

else:
    print("No syntentic genes found")

if os.path.exists(unique_genes):
    print("Unique genes metaplot")
    mC_bed = functions.allc2bed(allc)
    bed = pbt.BedTool(unique_genes).filter(functions.feat_filter,'mRNA').filter(functions.chr_filter,filter_chr)
    flank_bed = pbt.bedtool.BedTool.flank(bed,g=genome_file,l=2000,r=2000,s=True).saveas('f_tmp')
    cds_bed = bed.filter(functions.feat_filter,'CDS').filter(functions.chr_filter,filter_chr).saveas('c_tmp')
    bed = cds_bed.cat(flank_bed, postmerge=False)
    mapping = pbt.bedtool.BedTool.intersect(mC_bed,bed,wa=True)
    m = pd.read_table(mapping.fn, header=None, usecols = [0,1,5,6,7,8,9])
    m.columns = ['chr','pos','strand','mc_class','mc_count','total','methylated']
    m = m.drop_duplicates()
    m.to_csv('CDS_allc.tmp', sep='\t', index=False)

    functions.feature_metaplot('CDS_allc.tmp',unique_genes,genome_file,output="results/unique_gene_metaplot.tsv",
                               ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0,
                               filter_features='mRNA',filter_chr=filter_chr)

    for i in ['CDS_allc.tmp','c_tmp','f_tmp']:
        os.remove(i)

else:
    print("No unique genes found")
