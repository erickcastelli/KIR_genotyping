#!/usr/bin/perl
use strict;
use File::Basename;
use Getopt::Std;

our ($opt_t, $opt_b, $opt_o, $opt_i);
getopts('t:b:o:i:');

my $target = $opt_t;
my $bams = $opt_b;
my $out = $opt_o;
my $intervals = $opt_i;

if ($out eq "") {help();} 
if ($bams eq "") {help();} 
if ($target eq "") {help();} 
if ($intervals eq "") {help();} 

if (! -e $bams) {help();} 
if (! -e dirname $out) {help();} 
if (! -e $intervals) {help();} 

open (IN, $intervals);
my %interval = ();
while (<IN>)
{
	chomp $_;
	if ($_ eq "") {next;}
	my @data = split("\t",$_);
	$interval{$data[0]}{$data[1]} = $data[2];
}
close (IN);

if (! exists $interval{$target}) {print "\nThis target is not on the list!\n"; help();}

my @files = <$bams/*.bam>;
my @gene_list = sort keys $interval{$target};
my @references = sort keys $interval{"REF"};


my %db = ();
my $count = 0;
foreach (@files)
{
	my $sample = basename $_;
	my $file = $_;
	$sample =~ s/.adjusted.bam//;
	my $progress = ($count / scalar(@files)) * 100;
	$progress = sprintf("%.1f", $progress);
	print "\rProcessing sample $sample ... $progress % done          ";
	$| = 1;
	
	foreach (@gene_list)
	{
		my $cmd = `samtools coverage '$file' -r $interval{$target}{$_} --verbosity 0 --ff SECONDARY,UNMAP,DUP -q 1`;
		my @data = split("\n",$cmd);
		my @fields = split(" ",$data[1]);
		my $depth = $fields[6];
		$db{$sample}{$_} = $depth;
	}
	foreach (@references)
	{
		my $cmd = `samtools coverage '$file' -r $interval{"REF"}{$_} --verbosity 0 --ff SECONDARY,UNMAP,DUP -q 1`;
		my @data = split("\n",$cmd);
		my @fields = split(" ",$data[1]);
		my $depth = $fields[6];
		$db{$sample}{$_} = $depth;
	}	
	$count++;
}
print "\rProcessing samples ... 100 % done             ";
$| = 1;

print "\nAll samples done!";
print "\n";
$| = 1;

open (OUT, ">$out");

print OUT "SAMPLE";
foreach (@gene_list)
{
	print OUT "\t$_";
}
foreach (@references)
{
	print OUT "\t$_";
}
print OUT "\t$target\_average\treferences_average\tratio";
print OUT "\n";


my @samples = sort keys %db;
foreach (@samples)
{
	my $target = 0;
	my $ref = 0;
	
	print OUT $_;
	my $sample = $_;
	foreach (@gene_list)
	{
		print OUT "\t$db{$sample}{$_}";
		$target = $target + $db{$sample}{$_};
	}
	$target = $target / scalar(@gene_list);
	$target = sprintf("%.2f", $target);
	
	foreach (@references)
	{
		print OUT "\t" . $db{$sample}{$_};
		$ref = $ref + $db{$sample}{$_};
	}
	$ref = $ref / scalar(@references);
	$ref = sprintf("%.2f", $ref);
	
	my $ratio = $target / $ref;
	$ratio = sprintf("%.2f", $ratio);
	
	print OUT "\t$target\t$ref\t$ratio";
	print OUT "\n";
}
close (OUT);
print "Please check file $out for results.\n";
exit;

sub help
{
	print "\nPerl script to detect and compare coverage for KIR genes and references";
	print "\nVersion 1.0, written by Erick C. Castelli\n";
	print "\nHow to use it:";
	print "\nperl $0 -b bam_folder -o output_file -t target -i intervals.txt\n";
	print "\ntargets: just one per run. E.g.,KIR2DL4, KIR3DL2";
	print "\nbam_folder: A directory in which all hla-mapper BAMS are placed";
	print "\nintervals.txt: This file is provided with the script\n\n";
	
	exit;
}