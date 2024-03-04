
# Allele-specific m6A RNA modification analysis 

This tool facilitates allele-specific m6A RNA modification (ASM) analysis by identifying differential modification levels across alleles with Oxford Nanopore direct RNA sequencing. Aligned long-reads are allocated to their respective alleles using our Python script, while m6A modification levels are quantified with [m6Anet](https://m6anet.readthedocs.io/en/latest/), a recent supervised machine learning model from Goke lab designed for precise m6A detection. Subsequently, the data from each allele read group are analyzed to assess ASM using an R script designed for this purpose. This integrated approach offers a comprehensive solution for ASM analysis, available in this repository.

![github_scheme](https://github.com/DayeaPark/Allele-specific-m6A-modification/assets/99752377/aa39b0f8-e6bc-45d4-8dca-c160e8aa6196)

## Contents 
* Installation 
* Procedure 
* Reference file

## Installation 
The environment requires these tools. 
* python3
* guppy3 
* m6Anet 
* minimap2
* nanopolish
* samtools

### conda option
Install  Conda.
All other dependencies can be installed using the environment file, ASM_env.yaml, in this repository.

```
git clone https://github.com/DayeaPark/Allele-specific-m6A-modification.git  
conda env create -f env/ASM_env.yaml
```

## Procedure 

### 1. Basecalling
Fastq files are generated from basecalling the raw fast5 files using guppy basecaller (quality score cutoff = 7). The passed reads were concatenated in a fasta file and used for further steps.  

```
guppy_basecaller --input_path fast5_directory/ --save_path output_directory/ --flowcell FLO-MIN106 --kit SQK-RNA002 --calib_detect --num_callers 16 ----gpu_runners_per_device 80 -x "cuda:all" --compress_fastq --reverse_sequence --u_substitution

# The passed fastq files are concatenated in a fasta
cat output_directory/pass/*.fastq.gz > passed_all.fastq.gz
zless passed_all.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' > passed_all.fasta
```

### 2. Alignment 
Reads were aligned to the transcriptome using minimap2. For precise modification detection, we utilized an N-masked transcriptome reference where SNPs were obscured with the letter ‘N’. The process for generating this reference is detailed in the Reference File section. The mapping quality was filtered with 10. 

```
minimap2 -ax map-ont transcriptome.fa passed_all.fasta -o output.sam
samtools view -F0x900 -q 10 -h -Sb output.sam  | samtools sort > output_sort.bam
samtools index output_sort.bam
```

### 3. Allele assignment of individual long-read 
In the aligned BAM file, allelic origins of individual reads are determined using a developed Python script (allele_assignment.py). This script locates SNP positions within each read, accounting for any insertions or deletions relative to the reference transcriptome. For this operation, it necessitates an aligned BAM file and a VCF file as inputs. Users can modify the allelic bias ratio threshold (default=0.7) and set a minimum SNP count (default=2). The script generates three output files: alternative.bam (bias ratio > 0.7 in default setting), reference.bam (bias ratio < 0.3 in default setting), and undefined.bam (bias ratio between 0.3 and 0.7). Additionally, a log.txt file is produced, detailing the SNP counts for each single read analyzed (transcript name, alternative SNP count, and reference SNP count). Each bam file needs an index for next steps. 

```
python allele_assignment.py --threshold_alt 0.7  –threshold_snp 2 --bam output_sort.bam --vcf vcf_file > log.txt

samtools index reference.bam
samtools index alternative.bam 
samtools index undefined.bam 

```

### 4. Allelic m6A modification detection using m6Anet 

To detect m6A modifications, Nanopolish requires an index created from the basecalled FASTA file. Each BAM file is then processed individually using a specified command. Detailed instructions and further information can be found in the [m6Anet documentation](https://m6anet.readthedocs.io/en/latest/). 

```
nanopolish index -d passed_all.fasta
nanopolish eventalign -r passed_all.fasta -b [bam file] -g transcriptome.fa --scale-events --signal-index -t 50 > eventalign.txt

## dataprep
m6anet-dataprep --eventalign eventalign.txt --out_dir dataprep/ --n_processes 10 --readcount_min 2 --readcount_max 50000

## run inference 
m6anet inference --input_dir dataprep/ --out_dir m6Anet_output/ --n_processes 10 --num_iterations 1000

```
### 5. Identification of allele-specific m6A modifications (ASM)
ASM were assessed through bootstrapping strategy based on R. 

## Reference files 
In this analysis, we utilized modified transcriptome references and transcriptomic VCF. Heterogeneous SNPs, as identified in the VCF file, were either masked with 'N' (resulting in an N-masked transcriptome) or substituted with the corresponding SNP from either the alternative or reference allele (resulting in masked_alt_transcriptome and masked_ref_transcriptome, respectively). These modifications, including creation of transcriptomic VCF, were facilitated by custom Python scripts, which were developed based on the methodology outlined in the published paper by [Ozadam et al. (2023)](https://www.nature.com/articles/s41586-023-06228-9). These files are available in the reference directory of this repository.
