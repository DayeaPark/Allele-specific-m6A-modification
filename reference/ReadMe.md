### Generating transcriptome-specific VCF and masked transcriptome for ASM analysis

ASM analysis necessitates a modified transcriptome and a VCF file that has undergone location conversion from the genome to the transcriptome. Our custom Python script identifies heterozygous SNPs to produce a transcriptome-specific VCF and compiles a masked transcriptome in the FASTA format. Essential for file generation are the transcriptome reference (FASTA), transcript annotations (GTF), and genetic VCF files.
In this research, we utilized hybrid mouse embryonic stem cells (mESC from C57BL/6J and CAST/EiJ strains) and the human NA12878 cell line for ASM analysis. The references used in our study were created using scripts located in each directory.

Each data were downloaded from following: <br> <br>
Mouse: 
* Transcriptome reference and GTF (appris) https://appris.bioinfo.cnio.es/#/downloads
* Genomic VCF https://www.mousegenomes.org/

NA12878 : 
* Transcriptome reference and GTF (appris) https://appris.bioinfo.cnio.es/#/downloads
* Genomic VCF https://hgdownload.soe.ucsc.edu/gbdb/hg38/platinumGenomes/
