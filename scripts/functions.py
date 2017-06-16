import os
import sys
import pybedtools as pbt
import pandas as pd
import itertools
import numpy as np
from scipy.stats.stats import pearsonr
from collections import Counter
from Bio import SeqIO

#interpret sequence context, taken from methylpy.utils
def expand_nucleotide_code(mc_type=["C"]):
    iub_dict = {"N":["A","C","G","T"],"H":["A","C","T"],"C":["C"],"G":["G"],"T":["T"],"A":["A"]}
    for type in mc_type[:]:
        type += "N" * (3 - len(type))
        mc_type.extend(["".join(i) for i in itertools.product(*[iub_dict[nuc] for nuc in type])])
    if "C" in mc_type:
        mc_type.extend(["CG", "CHG", "CHH","CNN"])
    if "CG" in mc_type:
        mc_type.extend(["CGN"])
    return mc_type

#filter allc file based on sequence context
def filter_context(allc,context=["C"]):
    a = pd.read_table(allc)
    a = a[a.mc_class.isin(expand_nucleotide_code(context))]
    return a

def allc2bed(allc,context=["C"],bed=True):
    a = filter_context(allc,context)
    a['pos2'] = a.pos
    a['name'] = a.index
    a['score'] = "."
    a = a[['chr','pos','pos2','name','score','strand','mc_class','mc_count','total','methylated']]
    if bed is True:
        a = pbt.BedTool.from_dataframe(a)
    return a

#simple function for filtering gff files based on feature (gene, exon, mRNA, etc)
def feat_filter(x,feature):
    if feature:
        return x[2] == feature
    else:
        return x

#simple function for filtering gff files based on strand
def strand_filter(x,strand):
    return x.strand == strand

#simple function for filtering gff files based on chromosome
def chr_filter(x,chr):
    return x.chrom not in chr

#generic function for calculating methylation levels in windows
def window_methylation_levels(m,cutoff=0,filter=0.5,nuc_bed=()):
    a = pd.DataFrame(columns=['window','mCG','mCHG','mCHH'])
    name = "none"
    CG = mCG = CHG = mCHG = CHH = mCHH = 0
    if nuc_bed:
        nuc = pd.read_table(nuc_bed.fn, usecols = [3,7,8])
        m = pd.merge(m,nuc,left_on=13,right_on='4_usercol')
    for c in m.itertuples():
      if name == "none":
        name = c[5]
        if nuc_bed:
            GC = int(c[7]) + int(c[8])
        if int(c[3]) >= int(cutoff):
          if c[1].startswith("CN") or c[1].endswith("N"):
            continue
          elif c[1].startswith("CG"):
            CG = CG + int(c[3])
            mCG = mCG + int(c[2])
          elif c[1].endswith("G"):
            CHG = CHG + int(c[3])
            mCHG = mCHG + int(c[2])
          else:
            CHH = CHH + int(c[3])
            mCHH = mCHH + int(c[2])
      elif c[5] != name:
        if nuc_bed:
            if ((CG + CHG + CHH)/GC) >= filter:
                a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
        else:
            a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
        name = c[5]
        if nuc_bed:
            GC = int(c[7]) + int(c[8])
        CG = mCG = CHG = mCHG = CHH = mCHH = 0
        if int(c[3]) >= int(cutoff):
            if c[1].startswith("CN") or c[1].endswith("N"):
              continue
            elif c[1].startswith("CG"):
              CG = CG + int(c[3])
              mCG = mCG + int(c[2])
            elif c[1].endswith("G"):
              CHG = CHG + int(c[3])
              mCHG = mCHG + int(c[2])
            else:
              CHH = CHH + int(c[3])
              mCHH = mCHH + int(c[2])
      elif c[5] == name:
        if int(c[3]) >= int(cutoff):
          if c[1].startswith("CN") or c[1].endswith("N"):
            continue
          elif c[1].startswith("CG"):
            CG = CG + int(c[3])
            mCG = mCG + int(c[2])
          elif c[1].endswith("G"):
            CHG = CHG + int(c[3])
            mCHG = mCHG + int(c[2])
          else:
            CHH = CHH + int(c[3])
            mCHH = mCHH + int(c[2])
    if nuc_bed:
        if ((CG + CHG + CHH)/GC) >= filter:
            a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
    else:
        a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
    return a

#get correlations of mC with genes
def mC_distribution(genome_file,allc,gff_file,fasta,windows=100000,stepsize=50000,cutoff=0,filter=0.5,feature='gene',output_path=(),save_windows=False,return_windows=False):
    w_bed = pbt.bedtool.BedTool.window_maker(pbt.BedTool(genome_file),g=genome_file,w=windows,s=stepsize,i='srcwinnum')
    nuc = pbt.BedTool.nucleotide_content(w_bed,fasta)
    mC_bed = allc2bed(allc)
    gene_gff = pbt.BedTool(gff_file).filter(feat_filter,feature=feature)
    gene_mapping = pbt.bedtool.BedTool.intersect(w_bed,gene_gff,wa=True)
    g = pd.DataFrame.from_dict(Counter([(f[3]) for f in gene_mapping]), orient='index').reset_index()
    g.columns=['window','genes']
    allc_mapping = pbt.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True)
    m = pd.read_table(allc_mapping.fn,header=None,usecols=[10,13,6,7,8])
    m = m.sort_values(by = 13,ascending=True)
    a = window_methylation_levels(m,cutoff,filter,nuc)
    df = pd.merge(a,g,on='window')
    correlations = pd.DataFrame(columns=['Context','r','p-value'])
    corCG = pearsonr(df['mCG'],df['genes'])
    correlations = correlations.append({'Context': 'mCG', 'r': corCG[0], 'p-value': corCG[1]}, ignore_index=True)
    corCHG = pearsonr(df['mCHG'],df['genes'])
    correlations = correlations.append({'Context': 'mCHG', 'r': corCHG[0], 'p-value': corCHG[1]}, ignore_index=True)
    corCHH = pearsonr(df['mCHH'],df['genes'])
    correlations = correlations.append({'Context': 'mCHH', 'r': corCHH[0], 'p-value': corCHH[1]}, ignore_index=True)
    if output_path:
        correlations.to_csv(output_path + "correlations.tsv", sep='\t', index=False)
    if save_windows:
        df.to_csv(output_path + "genome_windows_data.tsv", sep='\t', index=False)
    if return_windows:
        return df
    else:
        return correlations

#plot methylation level for features
def feature_metaplot(allc,features,genome_file,output=(),ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0,filter_features=(),filter_chr=[]):
    counter = 1
    mC_bed = allc2bed(allc)
    a = pbt.BedTool(features)
    if ignoreStrand:
        p_bed = a.filter(feat_filter,filter_features).filter(chr_filter,filter_chr).saveas('p_tmp')
    else:
        p_bed = a.filter(strand_filter,strand='+').filter(feat_filter,filter_features).filter(chr_filter,filter_chr).saveas('p_tmp')
        n_bed = a.filter(strand_filter,strand='-').filter(feat_filter,filter_features).filter(chr_filter,filter_chr).saveas('n_tmp')
    CG = mCG = CHG = mCHG = CHH = mCHH = 0
    metaplot = pd.DataFrame(columns=['Bin','mCG','mCHG','mCHH'])
    for y in ['u','f','d']:
        if y == 'u':
            if ignoreStrand:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
            else:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
                nf_bed = pbt.bedtool.BedTool.flank(n_bed,g=genome_file,l=updown_stream,r=0,s=True)
                pw_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=int(windows/3),i='srcwinnum')
                nw_bed = pbt.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=int(windows/3),i='srcwinnum',reverse=True)
        elif y == 'f':
            if ignoreStrand:
                w_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=int(windows/3),i='srcwinnum')
            else:
                pw_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=int(windows/3),i='srcwinnum')
                nw_bed = pbt.bedtool.BedTool.window_maker(n_bed,b=n_bed,n=int(windows/3),i='srcwinnum',reverse=True)
        elif y == 'd':
            if ignoreStrand:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
            else:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
                nf_bed = pbt.bedtool.BedTool.flank(n_bed,g=genome_file,l=0,r=updown_stream,s=True)
                pw_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=int(windows/3),i='srcwinnum')
                nw_bed = pbt.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=int(windows/3),i='srcwinnum',reverse=True)
        if ignoreStrand:
            w_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=int(windows/3),i='srcwinnum')
        else:
            w_bed = pw_bed.cat(nw_bed, postmerge=False)
        mapping = pbt.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True)
        m = pd.read_table(mapping.fn, header=None, usecols = [13,6,7,8])
        for x in list(range(1,int(windows/3)+1)):
            for c in m.itertuples():
                if c[4].endswith("_"+str(x)):
                    if int(c[3]) >= int(cutoff):
                        if c[1].startswith("CN") or c[1].endswith("N"):
                            continue
                        elif c[1].startswith("CG"):
                            CG = CG + int(c[3])
                            mCG = mCG + int(c[2])
                        elif c[1].endswith("G"):
                            CHG = CHG + int(c[3])
                            mCHG = mCHG + int(c[2])
                        else:
                            CHH = CHH + int(c[3])
                            mCHH = mCHH + int(c[2])
            metaplot = metaplot.append({'Bin': counter,'mCG': (np.float64(mCG)/np.float64(CG)),
                                        'mCHG': (np.float64(mCHG)/np.float64(CHG)),
                                        'mCHH': (np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
            counter = counter + 1
            CG = mCG = CHG = mCHG = CHH = mCHH = 0
    os.remove('p_tmp')
    os.remove('n_tmp')
    if output:
        metaplot.to_csv(output, sep='\t', index=False)
    else:
        return metaplot

#map methylation to features
def map2features(allc,features,genome_file,updown_stream=2000,first_feature=(),second_feature=(),filter_chr=[]):
    bed = pbt.BedTool(features).filter(feat_filter,first_feature).filter(chr_filter,filter_chr)
    flank_bed = pbt.bedtool.BedTool.flank(bed,g=genome_file,l=updown_stream,r=updown_stream,s=True).saveas('f_tmp')
    cds_bed = pbt.BedTool(features).filter(feat_filter,second_feature).filter(chr_filter,filter_chr).saveas('c_tmp')
    bed = cds_bed.cat(flank_bed, postmerge=False)
    mC_bed = allc2bed(allc)
    mapping = pbt.bedtool.BedTool.intersect(mC_bed,bed,wa=True)
    m = pd.read_table(mapping.fn, header=None, usecols = [0,1,5,6,7,8,9])
    m.columns = ['chr','pos','strand','mc_class','mc_count','total','methylated']
    m = m.drop_duplicates()
    m.to_csv('CDS_allc.tmp', sep='\t', index=False)

#plot methylation levels for genes
def gene_metaplot(allc,features,genome_file,output=(),ignoreStrand=False,windows=60,
                  updown_stream=2000,cutoff=0,first_feature=(),second_feature=(),
                  filter_chr=[],remove_tmp=True):
    map2features(allc,features,genome_file,updown_stream,first_feature,second_feature,filter_chr)
    feature_metaplot('CDS_allc.tmp',features,genome_file,output,ignoreStrand,
                     windows,updown_stream,cutoff,first_feature,filter_chr)
    if remove_tmp:
        for i in ['CDS_allc.tmp','c_tmp','f_tmp']:
            os.remove(i)
    else:
        for i in ['c_tmp','f_tmp']:
            os.remove(i)
