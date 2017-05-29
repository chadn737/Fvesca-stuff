import sys
import os
import subprocess
import urllib.request

if not os.path.exists('data'):
    os.makedirs('data')

os.chdir('data')

if not os.path.exists('ref'):
    os.makedirs('ref')
if not os.path.exists('v2'):
    os.makedirs('v2')
if not os.path.exists('v4'):
    os.makedirs('v4')
if not os.path.exists('fastq'):
    os.makedirs('fastq')

os.chdir('ref')
if not os.path.exists('v2'):
    os.makedirs('v2')
if not os.path.exists('v4'):
    os.makedirs('v4')
os.chdir('v2')
if not os.path.exists('bowtie2'):
    os.makedirs('bowtie2'')
if not os.path.exists('misc'):
    os.makedirs('misc')
os.chdir('bowtie2')
urllib.request.urlretrieve('ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v2.0.a1/assembly/Fragaria_vesca_v2.0.a1_pseudomolecules.fasta.gz',
'Fragaria_vesca_v2.0.a1_pseudomolecules.fasta.gz')
os.chdir('../misc')
urllib.request.urlretrieve('ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v2.0.a1/genes/Fragaria_vesca_v2.0.a1.transcripts.gff3.gz',
'Fragaria_vesca_v2.0.a1.transcripts.gff3.gz')

os.chdir('../../fastq')
urllib.request.urlretrieve('ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX165/SRX1656919/SRR3286267/SRR3286267.sra',
                            'SRR3286267.sra')
subprocess.Popen('fastq-dump --gzip --split-3 SRR3286267.sra')
