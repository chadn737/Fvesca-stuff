#!/usr/local/apps/python/2.7.8/bin/python
from methylpy.call_mc import build_ref

#fasta file(s) of genome
#Multiple files like input_files=['chr1.fa','chr2.fa',...,'chrY.fa','chrL.fa'] should also work
input_files=['Fvesca_v4.fa']

#Prefix of output files
output='Fvesca_v4'

build_ref(input_files,output)
