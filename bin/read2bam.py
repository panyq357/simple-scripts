#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import uuid

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''\
Use fastp, bwa and samtools to generate BAM from FASTQ.
Unmapped and duplicated reads will be filtered out.
Suitable for ChIP-seq and BSA analysis.'''
)

parser.add_argument('-i', type=str, dest="r1_path", help='path to FASTQ file')
parser.add_argument('-I', type=str, dest="r2_path", help='path to the second FASTQ file of pair-end data')
parser.add_argument('-x', type=str, dest="index_prefix", help='prefix of bwa index')
parser.add_argument('-t', type=str, dest="thread_num", help='number of threads')
parser.add_argument('-o', type=str, dest="out_path", help='path to output BAM file')
parser.add_argument('-H', type=str, dest="html_report_path", help='path to output fastp html report file')

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

def get_common_prefix(s1, s2):
    common_prefix = ""
    for i, j in zip(s1, s2):
        if i != j:
            break
        common_prefix += i
    return(common_prefix)


common_prefix = get_common_prefix(args.r1_path, args.r2_path)
temp_fq = {
    "R1": os.path.join(tempfile.gettempdir(), str(uuid.uuid4()) + ".fq"),
    "R2": os.path.join(tempfile.gettempdir(), str(uuid.uuid4()) + ".fq"),
}

fastp = subprocess.run(
    f'''
        fastp \
            --in1 {args.r1_path} \
            --in2 {args.r2_path} \
            --out1 {temp_fq["R1"]} \
            --out2 {temp_fq["R2"]} \
            --html {args.html_report_path} \
            --json /dev/null
    ''',
    shell=True
)

bwa_samtools = subprocess.run(
    f'''
        bwa mem \
            -t {args.thread_num} \
            {args.index_prefix} \
            {temp_fq["R1"]} \
            {temp_fq["R2"]} | \
        samtools fixmate -u -m \
            - - | \
        samtools sort -u \
            -@{args.thread_num} \
            -T {args.temp_dir} | \
        samtools markdup -u \
            -@{args.thread_num} \
            - - |
        samtools view -b \
            -@{args.thread_num} \
            -F 1028 \
            -o {args.out_path}
    ''',
    shell=True
)

os.remove(temp_fq["R1"])
os.remove(temp_fq["R2"])

