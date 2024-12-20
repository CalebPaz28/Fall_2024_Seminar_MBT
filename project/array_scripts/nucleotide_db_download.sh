## This is for downloading a nucleotide database from NCBI
update_blastdb.pl --decompress nt_viruses

## Moving the sequence files to a desired location
mv nt_viruses.* /scratch/biol726310/BIOL7263_Genomics/array_project/databases/nucleotide_db

## Looping through the decompress all of the files
for file in nt_viruses.*.tar.gz; do
    tar -xzvf "$file"
done

## Making the database (not required with the above code)
#makeblastdb -in /scratch/biol726310/BIOL7263_Genomics/array_project/databases/nucleotide_db/nt_viruses -parse_seqids -blastdb_version 5 -title "NT Virus Database" -dbtype nucl -out /scratch/biol726310/BIOL7263_Genomics/array_project/databases/nt_virus_db