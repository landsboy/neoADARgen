#!/bin/bash
chr="$1"
start_seq="$2"
end_seq="$3"
k="$4"
genome="$5"
# path to create temp file of BED format for query
out_temp="sup/TEMP/TEMP${k}.bed"
# create the temp file with the BED entry
printf "%s\t%s\t%s\n" "$chr" "$start_seq" "$end_seq" > "$out_temp"
# extract the sequence and print it
total_seq=$(bedtools getfasta -fi "$genome" -bed "$out_temp" | grep -A 1 "^>" | tail -n 1 | tr -d "\n" | awk '{print toupper($0)}')
echo "$total_seq"
