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

os.chdir('fastq')
urllib.request.urlretrieve('ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX165/SRX1656919/SRR3286267/SRR3286267.sra',
                            'SRR3286267.sra')

subprocess.Popen('fastq-dump --gzip --split-3 SRR3286267.sra')
