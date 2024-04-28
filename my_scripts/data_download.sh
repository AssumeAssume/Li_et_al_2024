#!/bin/bash
### Xiufeng Li 2019/10/6
### 
usage() {
	echo -e "\nFunction: download the data from ebi according to PRJ accession\n"
    echo "Usage:"
    echo -e "bash data_download.sh -a XXXXXX.tsv -g 'grep' \n"
    echo -e "or data_download.sh -a XXXXXX.tsv -g 'grep' \n "
    echo "-a ,the PRJ accession txt file ."
    echo "-g ,grep character to select which file to download."
    echo -e "\t pass '$' to grep if you want to download all file"
    echo -e "It will download the fastq data in your current dir \n"
    exit -1
}
no_args="true"
while getopts 'ha:g:' OPT; do
    case $OPT in
        a)
            input="$OPTARG";;
        g)
            grep="$OPTARG";;

        h) usage;;

        ?) usage

    esac
    no_args="FALSE"
done

  # 用于 准确定位 $1 
shift $((OPTIND - 1))

[[ "$no_args" == "true" ]] && { usage; exit 1; }


#example 
# bash data_download.sh XXXXXX.tsv $

# pass "$" to grep if you want to download all file

ftp_column=`head -n 1 $input|tr "\t" "\n"|cat -n|awk '$2=="fastq_ftp"{print $1}'`

 cat $input|grep -Pi "${grep}"|awk -v ftp="$ftp_column" -F "\t" '{print $ftp}' |tr ";" "\n"|sed  -e 's/ftp.*uk//g'  |grep -v "fastq_ftp" >ASCP_DOWNLOAD_FILE
ascp -QT -l 300m -P33001 -k 1 \
	-i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
	--mode recv \
	--host fasp.sra.ebi.ac.uk \
	--user era-fasp \
	--file-list ASCP_DOWNLOAD_FILE \
    -v \
	.
