#!/usr/bin/perl -w
use strict;

########
# This script filters mapped reads for single cells using a mapping quality threshold.
# It also removes reads with high with high levels of mCH, indicating they may have not been bisulfite converted.
# Script takes a list of all samples in a file named sample_list.txt
# BAM files with mapped reads should be named <sample>.bam and be in the current directory


my $command; my $return; my @return;
my $suffix="_q30";
my $read_mch_threshold=0.7;
my $ch_threshold=3;
my $threshold=500000;


$command="ls -d \*bismark/ > sample_list.txt";
print "$command\n"; $return=system($command);

open sample_list, "sample_list.txt" or die $!;
open script_out, ">filter.sh" or die $!;
open log_out, ">bam_filter.log" or die $!;
print script_out "mkdir recycle\n";
my $sample; my @name; my $name; my @sample; my $um; my $m; my $high_mc_count;
my @data; my $tmp; my $sam_name;

## # Iterate through each sample (cell):
while (<sample_list>)
{
	chop $_; chop $_; $sample=$_;
	
	print "processing $sample\n";
	chdir $sample;


	### RUN QUALITY FILTERING ON BAM FILES
	@sample=split(/_bismark/,$sample);
	if (-e "$sample[0]\.bam")
	{
		$command="samtools view -bhq30 $sample[0]\.bam > ${sample[0]}$suffix.sam";
		print "$command\n"; $return=system($command);
		$sam_name="${sample[0]}$suffix";
	}
	else
	{
		$command="ls *_rmdup.bam";
		print "$command\n"; $return=qx/$command/; @return=split(/\s/,$return);
		@name=split(/\.bam/,$return[0]); 
		$command="samtools view -hq30 $return[0] > ${name[0]}$suffix.sam";
		print "$command\n"; $return=system($command);
		$sam_name="${name[0]}$suffix";
	}
	

	### FILTER OUT THE READS WITH POOR BISULFITE CONVERSION 
	$high_mc_count=0;
	open sam_in, "$sam_name\.sam" or die $!;
	open sam_out, ">$sam_name\_filtered.sam" or die $!;
	while (<sam_in>)
	{
		chop $_;
		if (substr($_,0,1) eq '@') { print sam_out "$_\n"; }
		else
		{
			$um=0; $m=0;
			@data=split(/\t/,$_);
			$tmp=$data[14]; $um=($tmp=~tr/hx//);
			$tmp=$data[14]; $m=($tmp=~tr/HX//);	
			
			if (($um+$m>=$ch_threshold)and(($m/($um+$m))>$read_mch_threshold)) { $high_mc_count++; } else { print sam_out "$_\n"; }
		}
	}
	close sam_out;
	close sam_in;


	### PRINT STATISTICS FOR SAMPLE
	print log_out "$sam_name\tremoved $high_mc_count reads with high mCH level\n";
	$command="samtools view -bh ${name[0]}$suffix\_filtered.sam > ${name[0]}$suffix\_filtered.bam";
	print "$command\n"; $return=system($command);
	$command="rm $sam_name\.sam ${name[0]}$suffix\_filtered.sam";
	print "$command\n"; $return=system($command);	
	
	$command="samtools view -c ${name[0]}$suffix\_filtered.bam";
	print "$command\n"; $return=qx/$command/; @return=split(/\s/,$return);	
	if ($return[0]<$threshold)
	{
		print script_out "mv $sample recycle\n";
	}
	chdir "../";
}	
close sample_list;
close script_out;
close log_out;

$command="bash filter.sh";
print "$command\n"; $return=system($command);
