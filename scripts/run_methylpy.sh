#PBS -S /bin/bash
#PBS -q batch
#PBS -N methylpy
#PBS -l nodes=1:ppn=12:rjsnode
#PBS -l walltime=480:00:00
#PBS -l mem=30gb

echo "Starting"
module load python/2.7.8
cd $PBS_O_WORKDIR
mkdir new old

#Uncompress fastq files
echo "Uncompressing fastq files"
cd fastq
gunzip AD22DYACXX_lane6-strawberry_151_Schmitz_Pool8_INDEX_mixed_R1.fastq.gz
cd ../

#Run methylpy
echo "run methylpy on old ref"
cd old/
mkdir allc reports
python ../../scripts/run_methylpy.py Fvesca_old "../fastq/*.fastq" "../ref/old/methylCseq/Fvesca_old" "10" "9" "ChrC" > reports/Fvesca_old_output.txt

#Format allc files
echo "Formatting allc files"
mv allc_* allc/
cd allc
mkdir tmp
head -1 allc_Fvesca_old_ChrL.tsv > tmp/header
for i in allc_Fvesca_old_*
do
sed '1d' "$i" > tmp/"$i"
done

tar -cjvf Fvesca_old_allc.tar.bz2 allc_Fvesca_old_*
rm allc_*
cd tmp
rm allc_Fvesca_old_ChrL.tsv allc_Fvesca_old_ChrC.tsv
cat header allc_* > ../Fvesca_old_allc_total.tsv
cd ..
rm -R tmp
tar -cjvf Fvesca_old_allc_total.tar.bz2 Fvesca_old_allc_total.tsv
cd ../

#Cleanup directory
echo "Cleaning up intermediate files"
rm *mpileup* *.bam *.bam.bai
cd ../

#Run methylpy
echo "run methylpy on new ref"
cd new/
mkdir allc reports
python ../../scripts/run_methylpy.py Fvesca_new "../fastq/*.fastq" "../ref/new/methylCseq/Fvesca_new" "10" "9" "ChrC" > reports/Fvesca_new_output.txt

#Format allc files
echo "Formatting allc files"
mv allc_* allc/
cd allc
mkdir tmp
head -1 allc_Fvesca_new_ChrL.tsv > tmp/header
for i in allc_Fvesca_new_*
do
sed '1d' "$i" > tmp/"$i"
done

tar -cjvf Fvesca_new_allc.tar.bz2 allc_Fvesca_new_*
rm allc_*
cd tmp
rm allc_Fvesca_new_ChrL.tsv allc_Fvesca_new_ChrC.tsv
cat header allc_* > ../Fvesca_new_allc_total.tsv
cd ..
rm -R tmp
tar -cjvf Fvesca_new_allc_total.tar.bz2 Fvesca_new_allc_total.tsv
cd ../

#Cleanup directory
echo "Cleaning up intermediate files"
rm *mpileup* *.bam *.bam.bai
cd ../

#Compress fastq files
echo "Compressing fastq files"
cd fastq
for i in *fastq
do
  gzip "$i"
done

echo "done"
