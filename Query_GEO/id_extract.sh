#!/bin/bash

<< script_info
	Gets UIDs of GEO datasets for Pseudomonas, Bacteroides, Staph and Strep
	Arguments: Output directory for generated text files of UIDs returned from GEO
	Usage: ./id_extract.sh <output directory>
	Change "reldate" argument in URLs to update how many days "back" to query
	Script execution order: id_extract, ftp_extract, suppl_data_download
script_info

outdir=$1

# get XML files from query
PseudomonasHits=$( exec curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=%22pseudomonas%20aeruginosa%22+AND+%22Pseudomonas%20aeruginosa%22%5Borganism%5D+AND+%22Expression%20profiling%20by%20high%20throughput%20sequencing%22%5BFilter%5D&reldate=885&datetype=pdat&retmax=1000&usehistory=y" --silent)
BacteroidesHits=$( exec curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=%22bacteroides%20thetaiotaomicron%22+AND+%22Bacteroides%20thetaiotaomicron%22%5Borganism%5D+AND+%22Expression%20profiling%20by%20high%20throughput%20sequencing%22%5BFilter%5D&reldate=885&datetype=pdat&retmax=1000&usehistory=y" --silent)
StaphylococcusHits=$( exec curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=%22staphylococcus%20aureus%22+AND+%22Staphylococcus%20aureus%22%5Borganism%5D+AND+%22Expression%20profiling%20by%20high%20throughput%20sequencing%22%5BFilter%5D&reldate=885&datetype=pdat&retmax=1000&usehistory=y" --silent)
StreptococcusHits=$( exec curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=%22streptococcus%20sanguinis%22+AND+%22Streptococcus%20sanguinis%22%5Borganism%5D+AND+%22Expression%20profiling%20by%20high%20throughput%20sequencing%22%5BFilter%5D&reldate=885&datetype=pdat&retmax=1000&usehistory=y" --silent)

HitList=( "$PseudomonasHits" "$BacteroidesHits" "$StaphylococcusHits" "$StreptococcusHits" )
IdVars=( "PseudomonasIds" "BacteroidesIds" "StaphylococcusIds" "StreptococcusIds" )
Species=("Pseudomonas" "Bacteroides" "Staphylococcus" "Streptococcus")

# parse for UIDs
Iter=0
for hit in "${HitList[@]}"; do
	declare "${IdVars[Iter]}=$( echo "$hit" | grep -E 'Id>[0-9]+<\/Id>' | sed 's/<Id>\(.*\)<\/Id>/\1/g')"
	((Iter++))
done

# Output the extracted Id values
IdList=( "$PseudomonasIds" "$BacteroidesIds" "$StaphylococcusIds" "$StreptococcusIds" )
Iter=0
for Ids in "${IdList[@]}"; do
	filename="${Species[Iter]}"_UIDs.txt
	echo "$Ids" > "$outdir"/"$filename"
	((Iter++))
done
