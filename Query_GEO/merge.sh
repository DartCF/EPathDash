#!/bin/bash
<< script_info
	Merges UIDs and GEO Accession IDs into single CSV file
	Arguments: Text file with UIDs, text file with FTP links, file location to save new table
		* Order of datasets must correspond across UID and FTP files
	Usage: ./merge.sh <UID_file> <FTP_file> <Output_file>
	Script execution order: id_extract, ftp_extract, merge.sh, get_metadata.sh
script_info

UIDfile=$1
FTPfile=$2
output=$3

Iter=1

echo Accession,UID >> "$output"

while IFS= read -r line; do
	GEOAcc=$( echo "$line" | cut -d"/" -f7 )	
	GEOUID=$( cat "$UIDfile" | head -n "$Iter" | tail -1 )
	echo "$GEOAcc","$GEOUID" >> "$output"
	((Iter++))
done < "$FTPfile"