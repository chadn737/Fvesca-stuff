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

#Unpack needed fasta files
cd ../misc
gunzip ChrL.fa.gz
cd ../data

#Build old ref indexes
cd ref/old/misc
gunzip Fvesca_ChrC.fa.gz
cd ../bowtie2
gunzip Fvesca_226_v1.1.fa.gz
cat Fvesca_226_v1.1.fa ../misc/Fvesca_ChrC.fa > tmp
gzip Fvesca_226_v1.1.fa
python ../../../../scripts/fix_fasta.py -i tmp -o Fvesca_old.fa
rm tmp
time samtools faidx Fvesca_old.fa
cut -f1,2 Fvesca_old.fa.fai > Fvesca_old.genome
cat Fvesca_old.fa ../../../../misc/ChrL.fa.gz > ../methylCseq/tmp
cd ../misc
gzip Fvesca_ChrC.fa
cd ../

cd methylCseq
python ../../../../scripts/fix_fasta.py -i tmp -o Fvesca_old.fa
rm tmp
time samtools faidx Fvesca_old.fa
python ../../../../scripts/build_old_Fvesca_methylCseq_index.py
cd ../../

#Build new ref indexes
cd new/misc
gunzip Fvesca_ChrC.fa.gz
cd ../bowtie2
gunzip F_vesca_V4.1.fasta.gz
cat F_vesca_V4.1.fasta ../misc/Fvesca_ChrC.fa > tmp
gzip F_vesca_V4.1.fasta
python ../../../../scripts/fix_fasta.py -i tmp -o Fvesca_new.fa
rm tmp
time samtools faidx Fvesca_new.fa
cut -f1,2 Fvesca_new.fa.fai > Fvesca_new.genome
cat Fvesca_new.fa ../../../../misc/ChrL.fa.gz > ../methylCseq/tmp
cd ../misc
gzip Fvesca_ChrC.fa
cd ../

cd methylCseq
python ../../../../scripts/fix_fasta.py -i tmp -o Fvesca_new.fa
rm tmp
time samtools faidx Fvesca_new.fa
python ../../../../scripts/build_new_Fvesca_methylCseq_index.py
cd ../../

#Wrap it up
cd ../../misc
gzip ChrL.fa
cd ../data
echo "Finished"
