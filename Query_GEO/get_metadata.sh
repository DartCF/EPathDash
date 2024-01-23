#!/bin/bash
<< script_info
	Gets metadata for given GEO datasets
	Arguments: CSV file with UIDs and GEO Accession numbers, directory location for retreived metadata
	Usage: ./get_metadata <CSV_file> <Output_directory>
	Script execution order: id_extract, ftp_extract, merge, get_metadata
script_info

TableFile=$1
outdir=$2

{
	read
	while IFS= read -r line; do
		# accessory variable declarations
		base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id="
		GEOAcc=$( echo "$line" | cut -d"," -f1 )
		echo GeoAccession:"$GEOAcc"
		GEOUID=$( echo "$line" | cut -d"," -f2 )
		echo GEOUID:"$GEOUID"
		outfile="$outdir"/"$GEOAcc"_metadata.csv
		request="$base""$GEOUID"

		# xml request and parse metadata
		XmlOutput=$( exec curl "$request" --silent)
		date=$( echo "$XmlOutput" | grep -Eo '[0-9]{4}\/[0-9]{2}\/[0-9]{2}' )
		title=$( echo "$XmlOutput" | grep title | cut -d">" -f2 | cut -d"<" -f1 )
		summary=$( echo "$XmlOutput" | grep 'Name="summary"' | cut -d">" -f2 | cut -d"<" -f1)
		link=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="$GEOAcc"
		echo Date,GEO Accession,Title,Description,Strain,Medium,Treatment,Genotype,Link >> "$outfile"
		echo "$date","$GEOAcc",\""$title"\",\""$summary"\",NA,NA,NA,NA,"$link" >> "$outfile"
		sleep 1
	done
} < "$TableFile"
