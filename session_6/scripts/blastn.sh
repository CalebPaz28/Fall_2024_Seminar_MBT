## blastn search from the command line
blastn -db /home/mbtoomey/BIOL7263_Genomics/Example_data/blastdb/human_rna_db -query transcripts.fasta -outfmt "6 qseqid sseqid stitle" -num_threads 20 -num_alignments 1 > TTC_rna_blast.tsv
