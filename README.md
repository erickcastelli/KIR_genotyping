# UNDER DEVELOPMENT

# KIR Genotyping, haplotyping, and allele calling from next-generation sequencing data
Tutorial for genotyping, haplotyping, and allele calling for KIR genes

Version 1.0 (May 25th, 2021)

Author: Erick C. Castelli (erick.castelli@unesp.br)

The advantage of this method is that it calls variants, haplotypes, and then KIR alleles from the phased VCF. Thus, you have many different levels of information.

This tutorial is suitable for whole-genome sequencing (WGS), whole-exome sequencing (WES), and amplicon sequencing. It has been tested with Illumina sequencing data. Please note that coverage is essential. We recommend a coverage of at least 30x for WGS, and 100x for WES and amplicons. Read size is also important, and you will get better results when dealing with a read size larger than 75 nucleotides. **Copy number evaluation is not suitable for amplicon sequencing.**

This procedure needs sample size. The minimum sample size we have tested is 150 samples. For a single-sample allele call, please use another method (PING). However, if you have a large sample size, this procedure is highly accurate and will provide you with SNPs (in the hg38 reference genome), haplotypes, allele calls, copy numbers, and allow you to detect new variants.


## How to cite this pipeline
This pipeline is described in:

KIR2DL4 genetic diversity in a Brazilian population sample: implications for transcription regulation and protein diversity in samples with different ancestry backgrounds. Immunogenetics volume 73, pages 227–241 (2021). doi 10.1007/s00251-021-01206-9

Immunogenetics of resistance to SARS-CoV-2 infection in discordant couples. MedRxiv 2021, doi 10.1101/2021.04.21.21255872

## Packages and software needed
- hla-mapper 4 (www.castelli-lab.net/apps/hla-mapper) 
- GATK 4 (https://gatk.broadinstitute.org/hc/en-us)
- GATK 3.8 or WhatsHap (https://whatshap.readthedocs.io/en/latest/)
- vcfx 2 (www.castelli-lab.net/apps/vcfx)
- phasex 0.8.2 (https://github.com/erickcastelli/phasex)
- samtools 1.12 (http://samtools.sourceforge.net)
- BWA 0.7.17 (https://sourceforge.net/projects/bio-bwa/files/)
- bcftools 1.12 (http://samtools.github.io/bcftools/)
- vcftools (http://vcftools.sourceforge.net)
- IGV (https://software.broadinstitute.org/software/igv/)
- Emboss (transeq)

## STEP 1: Using hla-mapper to get unbiased read mapping for KIR genes
This step is essential. You won't retrieve correct genotypes in KIR genes unless using a method tailored for these genes. We recommend the use of hla-mapper.

hla-mapper supports many genes in the LRC region, incluing KIR genes, LILRB1, LILRB2, LAIR1, and LAIR2.

Please check its website for instructions (www.castelli-lab.net/apps/hla-mapper)

There are two possible inputs for hla-mapper, a BAM file (step 1a), or FASTQ files (step 1b).
 
### STEP 1A (USING A BAM FILE): 
- Download a copy of the human reference genome (hg38) and prepare it for BWA.
- Prepare it for BWA:
> bwa index reference_genome
- Using BWA MEM, map your reads against the reference genome, and sort it:
> bwa mem reference_genome R1.fastq R2.fastq > sample.sam
> 
> samtools sort sample.sam > sample.bam
> 
> samtools index sample.bam

If you downloaded a BAM file from another source, be sure that it refers to the hg38 reference genome and skip the procedure above.

Now, you can optimize read mapping in the KIR region
> hla-mapper dna bam=sample.bam db=hla_mapper_database sample=Sample_Name output=output_folder
- You need to repeat this for each sample.
- You need to indicate a different output folder for each sample


### STEP 1B (USING FASTQ FILES): 
First, you can perform the analysis as step 1A, or you can proceed with this step directly.

> hla-mapper dna r1=R1.fastq.gz r2=R2.fastq.gz db=hla_mapper_database sample=Sample_Name output=output_folder
- You need to repeat this for each sample.
- You need to indicate a different output folder for each sample

## STEP 2 - Check some of the hla-mapper BAM files using IGV 
Using IGV, please check some of the hla-mapper outputted BAM files (Sample_Name.adjusted.bam). Make sure everything is OK. You can compare the original BAM (using BWA MEM) with the new one (using hla-mapper).

***You should perform the following steps for each KIR gene separatedly.***

## STEP 3 - Copy number evaluation
After the hla-mapper optimization of all samples, proceed with the calculation of copy number for each loci.

To do that, use the script provided in /copy_number, and follow the instruction in folder /copy_number or https://github.com/erickcastelli/KIR_genotyping/blob/main/copy_number/README.md

For the next step (variant calling)
