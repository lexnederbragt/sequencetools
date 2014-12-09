sequencetools
=============

Set of scripts etc. to work with DNA sequence and related files:

####assemblathon_stats.pl
Basic statistics on contigs and scaffolds. Modifed after https://github.com/ucdavis-bioinformatics/assemblathon2-analysis

####scaffoldgap2bed.py
Generates a bed-file for the gaps for a given fasta file
Splits the sequences on gaps of (default): 20 bases

####trf2gff.py
Converts the output from the TRF program to a gff file
Requires the `-ngs` flag set (available from TRF version 4.07b)
Really rough code, use at your own risk
Usage:

```
python trf2gff.py trf_outfile.dat >trf_outfile.gff
```

