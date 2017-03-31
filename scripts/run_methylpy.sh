#PBS -S /bin/bash
#PBS -q batch
#PBS -N methylpy
#PBS -l nodes=1:ppn=12:rjsnode
#PBS -l walltime=480:00:00
#PBS -l mem=50gb

echo "Starting"
module load python/2.7.8
cd $PBS_O_WORKDIR
mkdir v4 v2

#Uncompress fastq files
echo "Uncompressing fastq files"
cd fastq
gunzip AD22DYACXX_lane6-strawberry_151_Schmitz_Pool8_INDEX_mixed_R1.fastq.gz
cd ../

#Run methylpy
echo "run methylpy on v2 ref"
cd v2/
mkdir allc reports
python ../../scripts/run_methylpy.py Fvesca_v2 "../fastq/*.fastq" "../ref/v2/methylCseq/Fvesca_v2" "10" "4" "ChrC" > reports/Fvesca_v2_output.txt

#Format allc files
echo "Formatting allc files"
mv allc_* allc/
cd allc
mkdir tmp
head -1 allc_Fvesca_v2_ChrL.tsv > tmp/header
for i in allc_Fvesca_v2_*
do
sed '1d' "$i" > tmp/"$i"
done

tar -cjvf Fvesca_v2_allc.tar.bz2 allc_Fvesca_v2_*
rm allc_*
cd tmp
rm allc_Fvesca_v2_ChrL.tsv allc_Fvesca_v2_ChrC.tsv allc_Fvesca_v2_ChrM.tsv
cat header allc_* > ../Fvesca_v2_allc_total.tsv
cd ..
rm -R tmp
tar -cjvf Fvesca_v2_allc_total.tar.bz2 Fvesca_v2_allc_total.tsv
cd ../

#Cleanup directory
echo "Cleaning up intermediate files"
rm *mpileup* *.bam *.bam.bai
cd ../

#Run methylpy
echo "run methylpy on v4 ref"
cd v4/
mkdir allc reports
python ../../scripts/run_methylpy.py Fvesca_v4 "../fastq/*.fastq" "../ref/v4/methylCseq/Fvesca_v4" "10" "4" "ChrC" > reports/Fvesca_v4_output.txt

#Format allc files
echo "Formatting allc files"
mv allc_* allc/
cd allc
mkdir tmp
head -1 allc_Fvesca_v4_ChrL.tsv > tmp/header
for i in allc_Fvesca_v4_*
do
sed '1d' "$i" > tmp/"$i"
done

tar -cjvf Fvesca_v4_allc.tar.bz2 allc_Fvesca_v4_*
rm allc_*
cd tmp
rm allc_Fvesca_v4_ChrL.tsv allc_Fvesca_v4_ChrC.tsv allc_Fvesca_v4_ChrM.tsv
cat header allc_* > ../Fvesca_v4_allc_total.tsv
cd ..
rm -R tmp
tar -cjvf Fvesca_v4_allc_total.tar.bz2 Fvesca_v4_allc_total.tsv
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
