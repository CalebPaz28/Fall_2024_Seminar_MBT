## This is for downloading the entire virus protein database
wget "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%2810239%29&sort=organism_name+asc" -O /scratch/biol726310/BIOL7263_Genomics/array_project/databases/protein_db/virus_aa_seq.fasta


## This is for downloading a nucleotide database from NCBI
## update_blastdb.pl --decompress nt_viruses