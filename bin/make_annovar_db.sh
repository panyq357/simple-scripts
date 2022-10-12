#!/usr/bin/env bash

### Argument Parsing ###

show_help() {
cat << EOF
version: 0.1.0
author: panyq

usage: make_annovar_db.sh [-h] [-r REF_GENOME_FA] [-t GTF_FILEPATH] [-o OUTPUT_DIR]

Make a annotation database for ANNOVAR gene-based annotation using reference
genome FASTA file and GTF annotation file from EnsemblPlants.

options:
  -h                   show this help message and exit
  -f VCF_FILEPATH      path to VCF file
  -r REF_GENOME_FA     path to indexed reference genome FASTA file
  -o OUTPUT_FILEPATH   path to output table
  -V OUTPUT_DB_VERSION version of output database

dependencies:
  - gtfToGenePred
  - ANNOVAR's perl scripts (needs to be in PATH)
    - retrieve_seq_from_fasta.pl

example:
./make_annovar_db.sh
    -r Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz
    -t Oryza_sativa.IRGSP-1.0.54.chr.gtf.gz
    -o osdb
    -V OS

EOF
}

while getopts hr:t:o:V: opt
do
    case ${opt} in
        h) show_help ; exit ;;
        r) ref_genome_fa=${OPTARG} ;;
        t) gtf_filepath=${OPTARG} ;;
        o) output_dir=${OPTARG} ;;
        V) output_db_ver=${OPTARG} ;;
    esac
done

if [ $# -eq 0 ]
then
    show_help
    exit
fi

missing_argument=0
if [ -z ${ref_genome_fa} ]
then
    echo 'ERROR: missing argument for "-r REF_GENOME_FA"'
    missing_argument=1
fi
if [ -z ${gtf_filepath} ]
then
    echo 'ERROR: missing argument for "-t GTF_FILEPATH"'
    missing_argument=1
fi
if [ -z ${output_dir} ]
then
    echo 'ERROR: missing argument for "-o OUTPUT_DIR"'
    missing_argument=1
fi
if [ -z ${output_db_ver} ]
then
    echo 'ERROR: missing argument for "-V OUTPUT_DB_VER"'
    missing_argument=1
fi
if [ ${missing_argument} -eq 1 ]
then
    exit 1
fi

### Argument Parsing END ###

### Unzip input files if needed ###
if [ ${ref_genome_fa##*.} == "gz" ]
then
    old_ref_genome_fa=${ref_genome_fa}
    ref_genome_fa=$(mktemp make_annovar_db_temp.XXXXXX)
    echo "NOTICE: Decompressing input REF_GENOME_FA into tempfile: '${ref_genome_fa}'"
    gzip -dc ${old_ref_genome_fa} > ${ref_genome_fa}
    TEMPFILES="${TEMPFILES} ${ref_genome_fa}"
fi
if [ ${gtf_filepath##*.} == "gz" ]
then
    old_gtf_filepath=${gtf_filepath}
    gtf_filepath=$(mktemp make_annovar_db_temp.XXXXXX)
    echo "NOTICE: Decompressing input GTF_FILEPATH into tempfile '${gtf_filepath}'"
    gzip -dc ${old_gtf_filepath} > ${gtf_filepath}
    TEMPFILES="${TEMPFILES} ${gtf_filepath}"
fi
trap "rm -rf ${TEMPFILES}" EXIT
### Unzip END ###

### Make annotation database ###
mkdir ${output_dir}
refGene_filepath="${output_dir}/${output_db_ver}_refGene.txt"
gtfToGenePred -genePredExt \
    ${gtf_filepath} \
    ${refGene_filepath}
refGeneMrna_filepath="${output_dir}/${output_db_ver}_refGeneMrna.fa"
retrieve_seq_from_fasta.pl \
    --format refGene \
    --seqfile ${ref_genome_fa} \
    ${refGene_filepath} \
    --out ${refGeneMrna_filepath}
### Make annotation database END ###

echo "make_annovar_db.sh DONE"
exit
