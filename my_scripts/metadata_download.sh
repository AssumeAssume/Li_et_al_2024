#!/bin/bash
### Xiufeng Li 2019/10/6
### 

case "$1" in
    (-h|--help|?|'')
    echo -e "\nFunction:get the metadata from ebi\n"
    echo -e "Usage: bash metadata_download.sh PRJNA_ACCESSION\n\tor metadata_download.sh PRJNA_ACCESSION\n"
    echo -e "It will download the metadata in your current dir and name the file as PRJNAXXXXX.tsv\n"
    exit 0 
;;
esac


prj_accession=$1
url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${prj_accession}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,library_name,library_layout,library_strategy,run_alias,sample_alias,fastq_md5,fastq_ftp,submitted_ftp,sra_ftp,sample_title&format=tsv&download=true"

wget -c $url -O ${prj_accession}.tsv
