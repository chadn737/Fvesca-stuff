#PBS -S /bin/bash
#PBS -q batch
#PBS -N methylC_analysis
#PBS -l nodes=1:ppn=2:rjslab
#PBS -l walltime=480:00:00
#PBS -l mem=100gb

#cd $PBS_O_WORKDIR
echo "Starting"
echo "Preparing gene lists"
cd ref
gunzip H4_1-1Syntelogs_V2-V4.txt.gz
cd v4/misc
cut -f2 ../../H4_1-1Syntelogs_V2-V4.txt > tmp
fgrep -f tmp Fvesca_v4.gff > Fvesca_v4_syntelogs.gff
rm tmp
cd ../../v2/misc
cut -f1 ../../H4_1-1Syntelogs_V2-V4.txt > tmp
fgrep -f tmp Fvesca_v2.gff > Fvesca_v2_syntelogs.gff
rm tmp
cd ../../
gzip H4_1-1Syntelogs_V2-V4.txt
cd ../

#echo "Analyzing v4 genome"
#cd v4
#mkdir results figures_tables
#echo "Running Python on v4"
#python ../../scripts/methylC_analysis.py v4
#echo "Running R on v2"
#Rscript ../../scripts/methylC_analysis.R v4
#cd ../

echo "Analyzing v2 genome"
cd v2
mkdir results figures_tables
echo "Running Python on v2"
python ../../scripts/methylC_analysis.py v2
echo "Running R on v2"
Rscript ../../scripts/methylC_analysis.R v2
cd ../
