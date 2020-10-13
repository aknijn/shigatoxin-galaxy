# shigatoxin-galaxy
This tool implemented for use in a Galaxy instance (https://galaxyproject.org/) performs Shiga toxin sub-typing of Escherichia coli

requirements for the tool:  
python 3.7 (https://www.python.org/)  
perl 5.26.2 (https://www.perl.org/)  
perl-bioperl 1.7 (https://bioperl.org/)  
blast 2.9 (https://doi.org/10.1186/1471-2105-10-421)  
trimmomatic 0.39 (https://doi.org/10.1093/bioinformatics/btu170)  
spades 3.14 (https://doi.org/10.1089/cmb.2012.0021)  
skesa 2.3 (https://doi.org/10.1186/s13059-018-1540-z)  
fastqc 0.11.9 (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
muscle 3.8 (https://doi.org/10.1093/nar/gkh340)  
duk (https://www.osti.gov/servlets/purl/1016000/)  
fastq-pair (https://github.com/linsalrob/EdwardsLab/)  

The input can be of two types: raw reads (FASTQ, single end or paired end) or contigs (FASTA). In the case of contigs, the tool simply performs a blastn search against 
the Shiga toxin subtype database (STSTDB) from the Statens Serum Institut SSI and Technical University of Denmark DTU (https://doi.org/10.1128/JCM.00008-15, 
https://bitbucket.org/%7Bec84c234-a1e2-4442-8d73-bc3bdc479f29%7D/).
Instead, in the case of raw reads, the tool performs a quality assessment applying FastQC and then trimming by Trimmomatic followed by several operations aimed to construct stx 
consensus sequences on which the final blastn search will be performed against the forementioned database STSTDB (see flow chart below).  
The tool performs four different assemblies: 1. SPAdes on all contigs; 2. SKESA on all contigs; 3. SPAdes on contigs filtered using duk/fastq_pair against STSTDB; 
4. SKESA on contigs filtered using duk/fastq_pair against STSTDB. On each assembly a blastn search is performed against STSTDB, extracting the best matching contig with an 
e-value < 0.001 and an identity > 95%. The sequences of all four results are divided between stx1 and stx2 and each four type files put together, thus obtaining two multifasta 
files: stx1.fasta and stx2.fasta. To both files their corrisponding reference sequences from STSTDB are added and the two resulting files are aligned by MUSCLE. From these 
alignments the reference sequences are filtered out before all possible consensus sequences are reconstructed. These stx1 and stx2 consensus sequences are combined into a 
single multifasta file on which a blastn query is performed against the STSTDB, extracting the best matching sequence with an e-value < 0.001 and an identity > 95%. 
Furthermore, all sequences shorter than 1200 base pairs are filtered out.  
The tool outputs a report with a table of all matching references with the corrisponding values for the pident, length and positive parameters.  
A summary is given of the Shiga toxin subtypes that matched between 95% < identity <= 100% with the identity value indicated in parentheses in case of partial matches.  
Also a link is given to the FastQC web page (two in case of paired end reads).


![flow chart of the tool](https://github.com/aknijn/shigatoxin-galaxy/blob/master/stx.png?raw=true)
