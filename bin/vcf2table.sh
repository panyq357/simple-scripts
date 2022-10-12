#!/usr/bin/env bash

### Argument Parsing ###

# Default values
num_forks=4
filter_expression="QUAL<=100 | INFO/DP > 500 | INFO/DP < 10"
chromosome="1,2,3,4,5,6,7,8,9,10,11,12"
vep_species_name="oryza_sativa"
cache_version="54"
csq_fields="Consequence,Gene,EXON,INTRON,CDS_position,Protein_position,Amino_acids,Codons,DISTANCE"

show_help() {
cat << EOF
version: 0.1.0
author: panyq

usage: vcf2table.sh [-h] [-i VCF_FILEPATH] [-r REF_GENOME_FA] [-o OUTPUT_FILEPATH] ...

This script does the following.
  1. Use bcftools to filter VCF.
  3. Use ANNOVAR to annotate VCF.
  2. Use GATK's VariantsToTable to convert VCF to a table.
This script is suitable for subsequent PyBSASeq analysis.

Notes:
  1. Reference genome FASTA should be indexed using samtools faidx and
     GATK's CreateSequenceDictionary.
  2. VCF file should be indexed using tabix.

options:
  -h                   show this help message and exit
  -i VCF_FILEPATH      path to VCF file
  -r REF_GENOME_FA     path to indexed reference genome FASTA file
  -o OUTPUT_FILEPATH   path to output table
  -f NUM_FORKS         number of thread used in VEP annotation (default: ${num_forks})
  -e FILTER_EXPRESSION filter expression for bcftools (default: "${filter_expression}")
  -c CHROMOSOME        chromosome keeped (default: "${chromosome}")
  -s VEP_SPECIES_NAME  species name for VEP (default: "${vep_species_name}")
  -V CACHE_VERSION     VEP cache version (default: "${cache_version}")
  -F CSQ_FIELDS        VEP annotation fields in INFO/CSQ of VCF 
                         (default: "${csq_fields}")
  -H SUMMARY_HTML      path to VEP summary html file (default: "<OUTPUT_FILEPATH>_summary.html")

dependencies:
  - bcftools
  - gatk4
  - ensembl-vep
EOF
}

while getopts hi:r:o:f:e:c:s:V:F:H: opt
do
    case ${opt} in
        h) show_help ; exit ;;
        i) vcf_filepath=${OPTARG} ;;
        r) ref_genome_fa=${OPTARG} ;;
        o) output_filepath=${OPTARG} ;;
        f) num_forks=${OPTARG} ;;
        e) filter_expression=${OPTARG} ;;
        c) chromosome=${OPTARG} ;;
        s) vep_species_name=${OPTARG} ;;
        V) cache_version=${OPTARG} ;;
        F) csq_fields=${OPTARG} ;;
        H) summary_html=${OPTARG} ;;
    esac
done

if [ $# -eq 0 ]
then
    show_help
    exit
fi

if [ -z ${summary_html} ]
then
    summary_html="${output_filepath}_summary.html"
fi

variables=(vcf_filepath ref_genome_fa output_filepath)
descriptions=('-i VCF_FILEPATH' '-r REF_GENOME_FA' '-o OUTPUT_FILEPATH')

missing_argument=0
for i in $(seq 0 $(expr ${#variables[@]} - 1))
do
    if [ -z ${!variables[$i]} ]
    then
        echo "ERROR: missing argument for '${descriptions[$i]}'"
        missing_argument=1
    fi
done

if [ $missing_argument -eq 1 ]
then
    exit 1
fi
### Argument Parsing END ###

temp_vcf=$(mktemp vcf2table_temp.XXXXXX)
trap "rm -rf ${temp_vcf}" EXIT

bcftools view \
    -e "${filter_expression}" \
    -r "${chromosome}" \
    ${vcf_filepath} | \
vep --genomes --cache --offline --vcf --per_gene --force_overwrite \
    --format vcf \
    --fork ${num_forks} \
    --fields "${csq_fields}" \
    --species ${vep_species_name} \
    --cache_version ${cache_version} \
    --output_file ${temp_vcf} \
    --stats_file ${summary_html}

gatk VariantsToTable \
    -R ${ref_genome_fa} \
    -V ${temp_vcf} \
    -F CHROM -F POS -F REF -F ALT \
    -GF AD -GF DP -GF GQ -GF PL -F CSQ \
    -O ${output_filepath}

