#!/usr/bin/env python
from Bio import SeqIO
import sys
import os
import gzip

def filter_fastq(fastq, out):
    fq_id = []
    outfile = gzip.open(out, "wt")
    with gzip.open(fastq, "rt") as fq:
        fastq = SeqIO.parse(fq,"fastq")
        for fq in fastq:
            if fq.id not in fq_id:
                SeqIO.write(fq, outfile, "fastq")
            else:
                continue
    outfile.close()

    

if __name__ == '__main__':
    fastq = sys.argv[1]
    output = sys.argv[2]
    filter_fastq(fastq, output)
