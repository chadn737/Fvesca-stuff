#PBS -S /bin/bash
#PBS -q batch
#PBS -N prepare_references
#PBS -l nodes=1:ppn=2:rjsnode
#PBS -l walltime=480:00:00
#PBS -l mem=10gb

echo "Starting"

#Load required modules
#Required only for use on GACRC Sapelo cluster, delete/modify elsewhere
cd $PBS_O_WORKDIR
module load bowtie2/2.2.9
module load bowtie/1.1.1
module load samtools/1.2
module load picard/2.4.1
module load tophat/2.0.13
module load python/2.7.8

#Unpack and format needed fasta files
cd ../misc
gunzip ChrL.fa.gz
cd ../data/ref/organelles
gunzip F_vesca_mito.fasta.gz
gunzip Fragaria_vesca_f._alba_H4_chloroplast.fasta.gz
sed s/F_vesca_mito/ChrM/g F_vesca_mito.fasta > Fvesca_ChrM.fa
sed s/Fragaria_vesca_f._alba_2095/ChrC/g Fragaria_vesca_f._alba_H4_chloroplast.fasta > Fvesca_ChrC.fa
cd ..

#Build v2 ref indexes
echo "Preparing v2 refs"
cd v2/bowtie2
gunzip Fragaria_vesca_v2.0.a1_pseudomolecules.fasta.gz
cat Fragaria_vesca_v2.0.a1_pseudomolecules.fasta ../../organelles/Fvesca_ChrC.fa ../../organelles/Fvesca_ChrM.fa > tmp
gzip Fragaria_vesca_v2.0.a1_pseudomolecules.fasta
python ../../../../scripts/fix_fasta.py -i tmp -o Fvesca_v2.fa
rm tmp
time samtools faidx Fvesca_v2.fa
cut -f1,2 Fvesca_v2.fa.fai > Fvesca_v2.genome
cat Fvesca_v2.fa ../../../../misc/ChrL.fa > ../methylCseq/tmp
cd ../

cd methylCseq
python ../../../../scripts/fix_fasta.py -i tmp -o Fvesca_v2.fa
rm tmp
time samtools faidx Fvesca_v2.fa
python ../../../../scripts/build_v2_Fvesca_methylCseq_index.py
cd ../

cd misc/
gunzip Fragaria_vesca_v2.0.a1.transcripts.gff3.gz
grep -v \# Fragaria_vesca_v2.0.a1.transcripts.gff3 > Fvesca_v2.gff
gzip Fragaria_vesca_v2.0.a1.transcripts.gff3.gz
cd ../../

#Build v4 ref indexes
echo "Preparing v4 refs"
cd v4/bowtie2
ls
gunzip F_vesca_V4.1.fasta.gz
cat F_vesca_V4.1.fasta ../../organelles/Fvesca_ChrC.fa ../../organelles/Fvesca_ChrM.fa > tmp
gzip F_vesca_V4.1.fasta
python ../../../../scripts/fix_fasta.py -i tmp -o Fvesca_v4.fa
rm tmp
time samtools faidx Fvesca_v4.fa
cut -f1,2 Fvesca_v4.fa.fai > Fvesca_v4.genome
cat Fvesca_v4.fa ../../../../misc/ChrL.fa > ../methylCseq/tmp
cd ../

cd methylCseq
python ../../../../scripts/fix_fasta.py -i tmp -o Fvesca_v4.fa
rm tmp
time samtools faidx Fvesca_v4.fa
python ../../../../scripts/build_v4_Fvesca_methylCseq_index.py
cd ../

cd misc/
gunzip F_vesca_v4.1.6_makerStandard_woTpases_genemodels_new_gene_ids.gff.gz
grep -v \# F_vesca_v4.1.6_makerStandard_woTpases_genemodels_new_gene_ids.gff.gz > Fvesca_v4.gff
gzip F_vesca_v4.1.6_makerStandard_woTpases_genemodels_new_gene_ids.gff.gz
cd ../../

#Wrap it up
cd organelles
for i in *fa *fasta
do
  gzip "$i"
done
cd ../../../misc
gzip ChrL.fa
cd ../data
echo "Finished"
