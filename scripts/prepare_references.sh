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

#Build old ref indexes
echo "Preparing old refs"
cd old/bowtie2
gunzip Fvesca_226_v1.1.fa.gz
cat Fvesca_226_v1.1.fa ../../organelles/Fvesca_ChrC.fa ../../organelles/Fvesca_ChrM.fa > tmp
gzip Fvesca_226_v1.1.fa
python ../../../../scripts/fix_fasta.py -i tmp -o Fvesca_old.fa
rm tmp
time samtools faidx Fvesca_old.fa
cut -f1,2 Fvesca_old.fa.fai > Fvesca_old.genome
cat Fvesca_old.fa ../../../../misc/ChrL.fa > ../methylCseq/tmp
cd ../

cd methylCseq
python ../../../../scripts/fix_fasta.py -i tmp -o Fvesca_old.fa
rm tmp
time samtools faidx Fvesca_old.fa
python ../../../../scripts/build_old_Fvesca_methylCseq_index.py
cd ../../
ls

#Build new ref indexes
echo "Preparing new refs"
cd new/bowtie2
ls
gunzip F_vesca_V4.1.fasta.gz
cat F_vesca_V4.1.fasta ../../organelles/Fvesca_ChrC.fa ../../organelles/Fvesca_ChrM.fa > tmp
gzip F_vesca_V4.1.fasta
python ../../../../scripts/fix_fasta.py -i tmp -o Fvesca_new.fa
rm tmp
time samtools faidx Fvesca_new.fa
cut -f1,2 Fvesca_new.fa.fai > Fvesca_new.genome
cat Fvesca_new.fa ../../../../misc/ChrL.fa > ../methylCseq/tmp
cd ../

cd methylCseq
python ../../../../scripts/fix_fasta.py -i tmp -o Fvesca_new.fa
rm tmp
time samtools faidx Fvesca_new.fa
python ../../../../scripts/build_new_Fvesca_methylCseq_index.py
cd ../../

#Wrap it up
cd organelles
for i in *fa *fasta
do
  gzip "$i"
done
cd ../../misc
gzip ChrL.fa
cd ../data
echo "Finished"
