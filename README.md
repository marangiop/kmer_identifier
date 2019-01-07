# kmer_identifier
This is the repository for the Bioinformatics Programming and System Management coursework as part of the MSc in Synthetic Biology and Biotechnology.

Kmer_Identifier is a Python-based programme facilitates the identification of short nucleotide subsequences (i.e. kmers) with unexpected distributions in biological data generated from Next generation sequencing platforms. It is specifically designed for analysing FASTQ-formatted files, one of the most common and most accepted file format for storing the output of high-throughput next generation sequencing platforms. 

Inputs

Kmer length = integer that specifies the length of the subsequence in number of nucleotide bases that is used to check the sequencing data against (i.e. kmer length of 2 = 2 bases long kmers ; 5 = 5 bases long kmers).

Reporting threshold = integer that allows to specify how many times does an observed count for a given kmer have to be greater than the expected count for that kmer such that the kmer is reported in a given run of Kmer_Identifier. For example, a reporting threshold of 2 means that Kmer_Identifier will only report kmer sequences that have an observed count at least 2 times greater than it would be expected if that kmer was at that position just by chance.
