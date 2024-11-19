#!/bin/bash

# Define output file
output_file="reps.fa"

# Start with an empty output file
> $output_file

# Fetch a list of human families in summary format
curl -s "https://dfam.org/api/families?format=summary&clade=Homo%20sapiens" | jq -r '.results[].accession' | while read accession; do
  # Fetch family details including the consensus sequence
  response=$(curl -s "https://dfam.org/api/families/$accession")

  # Extract family name and sequence
  family_name=$(echo "$response" | jq -r '.name')
  sequence=$(echo "$response" | jq -r '.consensus_sequence')

  # Append to output file in FASTA format
  echo -e ">$family_name\n$sequence" >> $output_file
done

echo "Consensus sequences saved to $output_file"