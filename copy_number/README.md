# Estimation of Copy Numbers for KIR genes using hla-mapper
Version: 1.0, May 26th, 2021
Written by: Erick C. Castelli

This Perl script calculates the ratio between the read depths observed for a target (the KIR gene) and a reference (the TNF/LTA/LTB cluster on chromosome 6).

The intervals for the regions that we consider is at file **intervals.txt**. You may edit this file to adapt for your own needs, but use word REF to indicate the references.

To calculate copy number, firt you need to optimize the alignments using hla-mapper (www.castelli-lab.net/apps/hla-mapper).

After the alignment optimization, copy all the ".adjusted.bam" files to a single location.

Use the Perl script **compare_kir_coverage.pl** to calculate read depths and ratios, as follows KIR2DL1 as an example:

> perl compare_kir_coverage.pl -i intervals.txt -b the_directory_with_all_bams -t KIR2DL1 -o output_TXT_file -i

The output_TXT_file is a file containing the read depths observed in each interval, the averages for the target and reference, and the ratio between them.


