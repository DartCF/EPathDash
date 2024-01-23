#!/bin/bash

<< script_info
	Gets supplementary files stored at given FTP locations from GEO server 
	Arguments: Text file with FTP links to query, location to save supplmentary files
	Usage: ./suppl_data_download.sh <FTP_file> <output_directory>
	Script execution order: id_extract, ftp_extract, suppl_data_download
script_info

URLFile=$1
outdir=$2

while IFS= read -r line; do
	GEOAcc=$( echo "$line" | cut -d"/" -f7 )
	mkdir "$outdir"/"$GEOAcc"
	echo "$GEOAcc"
	request="$line"/suppl/
	XmlOutput=$( exec curl "$request" --silent )
	filename=$( echo "$XmlOutput" | grep 'href="GSE' | cut -d'"' -f2)
	for f in $filename; do
		curl "$request""$f" --output "$outdir"/"$GEOAcc"/"$f"
		sleep 1
	done
done < "$URLFile"