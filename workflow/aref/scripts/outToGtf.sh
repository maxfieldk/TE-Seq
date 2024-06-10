#!/usr/bin/env bash
#this is the repeatmasker.out file
RMOUT=$1
cat ${RMOUT} | tail -n +4 | \
awk -v OFS="\t" '{ if ($12 ~ /)/) print $5, "RepeatMasker", "exon", $6, $7, $1, $9, ".", "gene_id \""$11"/"$10"/"$5"_"$15"\"; transcript_id \""$11"/"$10"/"$5"_"$15"\"; Target \""$11"/"$10" "$14" "$13" "$12"\"; pctdiv "$2"; pctdel "$3"; "; \
else print $5, "RepeatMasker", "exon", $6, $7, $1, $9, ".", "gene_id \""$11"/"$10"/"$5"_"$15"\"; transcript_id \""$11"/"$10"/"$5"_"$15"\"; Target \""$11"/"$10" "$12" "$13" "$14"\"; pctdiv "$2"; pctdel "$3"; " }' | \
awk -v FS="\t" -v OFS="\t" '{ if ($7 == "C") print $1, $2, $3, $4, $5, $6, "-", $8, $9; else print $0 }' | \
sort -k1,1 -k4,4n -k5,5n > $2