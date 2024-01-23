#!/bin/bash

<< script_info
	Gets FTP file links for GEO datasets for Pseudomonas, Bacteroides, Staph and Strep 
	Arguments: text file with UIDs, output file for FTP links
	Usage: ./ftp_extract.sh <UID_file> <output_file>
	Script execution order: id_extract, ftp_extract, suppl_data_download
script_info

UIDFile=$1

output_file=$2

touch "$output_file"

while IFS= read -r line; do
	base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id="
	request="$base""$line"
	XmlOutput=$( exec curl "$request" --silent)
	ftpUrl=$( echo "$XmlOutput" | grep -o 'ftp:\/\/ftp\.ncbi\.nlm\.nih\.gov\/geo\/series\/GSE[0-9]*[n]*\/GSE[0-9]*\/' | sed 's/ftp:/https:/' )
	echo "$ftpUrl" >> "$output_file"
	sleep 1
done < "$UIDFile"