#!/usr/bin/perl

# Addition of the -rl (read length) parameter

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

# variables
my $output = 'output_dir';
my $version = "1.1";
my $last_update = "2021-01-15";
my $help;
my $annotation_files;
my ($query_name, $E_value, $score);
my $chave = 0;
my (%qname, %evalue, %scores);
my (%hnumseq, %hfamily);
my $hmm_db;
my ($usr_score,$usr_evalue);
my ($numseq, $family);
my @virfamilies;
my $hmmsearch_tabular_file;
my @score;
my $transeq_file;
my $fasta_qual;
my $cpu;
my $input_file;
my $file_log = "file.log";
my $version_op;
my $replace = "yes";
my $read_length;
my $help_print = "##### HMM-Prospector - version $version ($last_update) - L. Oliveira & A. Gruber #####

HMM-Prospector is a script that uses a single or multiple profile HMMs as a query in 
similarity searches against a FASTQ/FASTA dataset using hmmsearch program. HMM-Prospector then processes 
the results and generated tabular files with qualitative and quantitative results. Previously run 
hmmsearch result files (short tabular format) can also be used as datasets. 


Usage:
hmm-prospector.pl -d <file> -s|-e <decimal>  <optional parameters>

Example:
hmm-prospector.pl -d all_reads.fastq -i model.hmm -s 10
hmm-prospector.pl -d all_reads.fastq -i vFam-A_2014.hmm -s 10 -a /usr/local/genome/databases/vfam/annotationFiles_2014/
hmm-prospector.pl -d transeq_result.fasta -i vFam-A_2014.hmm -s 1 -a  /usr/local/genome/databases/vfam/annotationFiles_2014/
hmm-prospector.pl -d hmmsearch_tabular.txt -e 0.001 -a /usr/local/genome/databases/vfam/annotationFiles_2014
hmm-prospector.pl -d hmmsearch_tabular.txt -e 0.001 

Mandatory parameters:
-d  <file>        : Dataset (FASTQ, FASTA or hmmsearch's tabular output file)

OPTIONAL PARAMETERS:
-a                : Directory containing profile HMM annotations (valid only when using vFam models as input).
-cpu              : Number of threads to be used by hmmsearch.
-e | -s <decimal> : E-value (-e) or score (-s) threshold value. Report hmmsearch hits that present values equal to 
		    or lower than E-value or equal to or larger than score. One parameter and the respective value 
		    must be provided. If an hmmsearch tabular result file is used as input, then parameters -e or -s 
		    become mandatory.
-h|help           : Show this help message.
-i                : Input file (single or multiple profile HMMs) - mandatory when using a FASTA or FASTQ dataset.
-o             	  : Output directory (default = output_dir).
-r 		  : Ignore cutoff scores in the profile HMMs and use a custom value defined by parameters -e or -s 
		    for all input models (default = yes). If -r no is used, HMM-Prospector will use the cutoff scores 
		    specified in the respective CUTOFF SCORE tag of each profile HMM. For models not containing cutoff 
		    values, HMM-Prospector will use the cutoff value specified by the parameter -e ou -s. If none of these 
		    parameters has been specified, the program will then use hmmsearch's default cutoff value (-E 10). 
-rl		  : Length of the reads contained in the input file.
-v|version        : Version.
 \n";

my $optret = GetOptions ("d=s"			=> \$input_file,
                         "i=s"  	        => \$hmm_db,
                         "s=f"                  => \$usr_score,
                         "e=f"                  => \$usr_evalue,
                         "o=s"		        => \$output,
                         "h|help"               => \$help,
			 "p=s"			=> \$transeq_file,
			 "a=s" 			=> \$annotation_files,
			 "cpu=i"		=> \$cpu,
			 "r=s"			=> \$replace,
			 "rl=i"			=> \$read_length,
			 "v|version"            => \$version_op);

if($help){
    print $help_print;
    die "\n";    
}

if($version_op){
    die "Version $version.\nLast update: $last_update\n";
}

if(!$input_file){
    die "ERROR: Missing mandatory argument -d.\n$help_print\n";
}

if(($usr_score) and ($usr_evalue)) {
    die "ERROR: Please use -s or -e, not both.\n$help_print\n";
}

if(!$cpu){
    $cpu = `grep CPU -c /proc/cpuinfo`;
    $cpu = $cpu/2;
}

my $taxon_tag = 0;
my $protein_tag = 0;
my $range_tag = 0;
my $type_tag = 0;

# This routines counts the number of profile HMMs.
if(defined $hmm_db){
    $taxon_tag = `grep -c TAXON $hmm_db`;
    $protein_tag = `grep -c PROTEIN $hmm_db`;
    $range_tag = `grep -c RANGE $hmm_db`;
    $type_tag = `grep -c TYPE $hmm_db`;
} 
my %hash_att = ();

# If more than one profile HMM is informed, it is necessary to store the individual 
# characteristics of each one
if($taxon_tag > 0){
     $/ = "//\n";      
     open(DATA, "<$hmm_db");
     my $name;
     foreach my $unique_hmm (<DATA>){
      	my @hmm_lines = split("\n", $unique_hmm);
	my $taxon = "";
	my $protein = "";
	my $range = "";
	my $type = "";
	my $name;
        foreach my $line (@hmm_lines){
	    if ($line =~ /NAME/g) {
                $name = $line;
                $name =~ s/NAME//;
                $name =~ s/^\s+|\s+$//g;
            }
            if ($line =~ /TAXON/g) {
		$taxon = $line;
		$taxon =~ s/TAXON//;
		$taxon =~ s/^\s+|\s+$//g;
            }
	    if ($line =~ /PROTEIN/g) {
                $protein = $line;
                $protein =~ s/PROTEIN//;
                $protein =~ s/^\s+|\s+$//g;
            }
	    if ($line =~ /RANGE/g) {
                $range = $line;
                $range =~ s/RANGE//;
                $range =~ s/^\s+|\s+$//g;
            }
	    if ($line =~ /TYPE/g) {
                $type = $line;
                $type =~ s/TYPE//;
                $type =~ s/^\s+|\s+$//g;
            }
        }
     	if($taxon eq ""){
	    $taxon = "-";
	}
 	if($protein eq ""){
            $protein = "-";
        }
	if($range eq ""){
            $range = "-";
        }
	if($type eq ""){
            $type = "-";
        }
	$hash_att{$name} = "$taxon\t$protein\t$range\t$type";
    }
    close (DATA);
}

$/ = "\n";

# Creating the output directory
$output = output_dir_name($output);
system "mkdir $output";
$file_log = $output."/".$file_log;

my $aux_score;
if(defined $usr_score){
    $aux_score = $usr_score;
}

# Creating the log file
my $log = "$output/error.log";
open(my $log_file_handle, ">$log") or die "ERROR: Could not create log file $log : $!\n";

my $transeq_name;
my $len = length($input_file);
my $dir = undef;
my $pos = rindex($input_file, "/");
if($pos != -1){
    $dir = substr $input_file, 0, $pos;
}
my $aux_name = substr $input_file, (rindex($input_file, "/") + 1), $len;
open(my $fl, ">$file_log") or die "ERROR: Could not create log file $file_log : $!\n";
print $fl "Dataset: $input_file\n";
print STDERR "Dataset: $input_file\n";

my $file_name = $aux_name;
my $aux;
my %hmms = ();
my %hmms_score = ();
my $resp;

#Checks if the input file exists.
if(-e $input_file){
    if($hmm_db){
    	print $fl "Input file: $hmm_db\n";
    	print STDERR "Input file: $hmm_db\n";
    }

    #Checks if an annotation file has been entered.
    if(defined $annotation_files){
    	my $caracter = substr $annotation_files, -1;
    	if($caracter eq '/'){
            my $l = length($annotation_files);
            $annotation_files = substr $annotation_files, 0, ($l-1);
    	}
    	print STDERR "Annotation data: $annotation_files\n";
    	print $fl "Annotation data: $annotation_files\n";
    }
    else{
    	print STDERR "Annotation data not specified. Generating quantitative reports without annotation data\n";
    	print $fl "Annotation data not specified. Generating quantitative reports without annotation data\n";
    }
    print $fl "File: $input_file\n";

    #Identifies the type of input file (fastq, fasta or tabular file).
    my $ext = verify_file_type($input_file);
    if($ext == 1){ #The input file is a fasta file  
	if(!defined $hmm_db){
	    print $fl "ERROR: Missing argument -i\n";
	    close($fl);
	    system "mv $file_log $output";
	    die "ERROR: Missing argument -i.\n$help_print\n";	
   	}
	else{
	    # Checks if the user entered a score or e-value to be used as a cutoff parameter.
            $resp = verifyScoreEvalue($hmm_db);
        }

	# Verifies the composition of the fasta file (nucleotide or amino acid)
        
	my $type = verifyFastaFileComposition($input_file);
	if($type == 1){ # Nucleotide
	    print $fl "Type: Nucleotide Fasta\n";
	    
	    # Runs transeq program to translate the nucleotide sequences into proteins	    
	    $transeq_name = runTranseq($file_name, $input_file, $dir);	
     	}
    	else{# Amino acid
	    print $fl "Type: Protein Fasta\n";
	    $transeq_name = $input_file;	    
    	}	

	#Runs hmmsearch program.
	$aux = $file_name."_hmmsearch.tab";
	$hmmsearch_tabular_file = runHmmsearch($file_name, $transeq_name, $resp);
    }
    elsif($ext == 2){ # The input file is a fastq file
	print $fl "Type: Fastq\n";
	print $fl "Step: Converting fastq file to fasta file.\n";
	
	if(!defined $hmm_db){
	    print $fl "ERROR: Missing argument -i\n";
            close($fl);
	    system "mv $file_log $output";
            die "ERROR: Missing argument -i.\nTry -h for more details.\n";
        }
	else{
	    # Checks if the user entered a score or e-value to be used as a cutoff parameter.
            $resp = verifyScoreEvalue($hmm_db);
        }
	my $aux1;
 	my $aux2;
	my $aux3;
	my $aux4;
  	my $do;

	# Generates a fasta file from fastq.
	if(defined $dir){
            $aux1 = $dir."/".$file_name.".fastq.fasta";
	    $aux2 = $dir."/".$file_name.".fasta";
	    $aux3 = $file_name.".fastq.fasta";
            $aux4 = $file_name.".fasta";	
	    if((-e $aux1) or (-e $aux2) or (-e $aux3) or (-e $aux4)){
            	print $fl "FASTA file found. Skipping FASTQ to FASTA conversion.\n";
            	print "FASTA file found. Skipping FASTQ to FASTA conversion.\n";
            }
	    else{
		$do = 1;
	    }
    	}
    	else{
	    $aux1 = $file_name.".fastq.fasta";
            $aux2 = $file_name.".fasta";
    	    if((-e $aux1) or (-e $aux2)){
	    	print $fl "FASTA file found. Skipping FASTQ to FASTA conversion.\n";
 	    	print "FASTA file found. Skipping FASTQ to FASTA conversion.\n";
    	    }
	    else{
		$do = 1;
	    }
	}
    	if(defined $do){	   
	    my $aux = $file_name.".fasta"; 
	    print $fl "Runnig fastq_to_fasta program (Parameters: -i $input_file -n -Q 33 -o $aux)\n";
	    print STDERR "Runnig fastq_to_fasta program (Parameters: -i $input_file -n -Q 33 -o $aux)\n";
	    system "fastq_to_fasta -i $input_file -n -Q 33 -o $aux 2>> $log";	    
	    print $fl "Done.\n";
    	}

	# Runs transeq program to translate the nucleotide sequences into proteins
	$transeq_name = runTranseq($file_name, $input_file, $dir);
  	
	# Runs hmmsearch program.
	$hmmsearch_tabular_file = runHmmsearch($file_name, $transeq_name, $resp);
	
    }
    else{ #The input file is a tabular file

	#Checks if the tabular file is a hmmsearch output tabular file.
	my $type = verifyTabularFile($input_file);	
	print $fl "Type: tabular file\n";
	$resp = 0;
	if($type == 1){ #The file is a hmmsearch output tabular file. 
	    $hmmsearch_tabular_file = $input_file;	    
	    if((!($usr_score)) and (!($usr_evalue))) {
              	system "rm -rf $output";
               		die "ERROR: Missing mandatory argument -s or -e.\n$help_print\n";
            }
	}
	else{ # The file is not a hmmsearch output tabular file. Execution was aborted.
	    print $fl "The $input_file file is not compatible with hmmsearch tabular file!\n";
	    die "The $input_file file is not compatible with hmmsearch tabular file!\n";	    
	}
    }
}
else{ # The input file does not exist. Execution was aborted.
    close($fl);
    system "rm -rf $output";
    die "File $input_file not found! Directory $output removed!\n";
}

print $fl "Creating table with results and saving it in $output\n"; 

my %vfams;
# If annotation files are entered, the program stores information of the families related 
# to each vFam model.
if(defined $annotation_files){
    #Associates vFAMs to their families
    my $filestring = `ls $annotation_files`;
    my @filenames = split(/\n/, $filestring);
    foreach my $name (@filenames){	
    	my $file = $annotation_files."/".$name;
    	my $num_seq;
    	open(my $fn, "$file");
    	while(<$fn>){
	    chomp($_);
	    if($_ =~ /NUM_SEQ/){
	    	$num_seq = $_;
	    	$num_seq =~ s/NUM_SEQ//g;
	    	$num_seq =~ s/'//g;
            	$num_seq =~ s/^\s+//;
	    }
            if($_ =~ /FAMILIES/){
	    	my $families = $_;
            	$families =~ s/FAMILIES//g;
            	$families =~ s/'//g;
            	$families =~ s/^\s+//;
	    	my $vfam_name = $name;
	    	$vfam_name =~ s/_annotations.txt//g; 
            	if($families eq "{}"){
		    $vfams{$vfam_name} = $num_seq."/"."{}";
	    	}
	    	else{
		    $families =~ s/{//g;
                    $families =~ s/}//g;
                    $families =~ s/^\s+//;
		    my @aux1 = split(",", $families);
		    my $families_name = "";
		    foreach my $fam(@aux1){
		    	my @aux2 = split(":", $fam);
		    	$families_name .= $aux2[0]."_";
		    }
		    $vfams{$vfam_name} = $num_seq."/".$families_name;
	    	}	
	    }
    	}
    	close($file);
    }
}

# The program analyzes the results of hmmsearch and generates the output files.

open (my $filehadlefromtabular, $hmmsearch_tabular_file) or die "ERROR: Could not read file $hmmsearch_tabular_file : $!\n";
my @lines = <$filehadlefromtabular>;

my $str_table1_vfam = "";
my $str_table2_vfam = "";
my $str_table3_vfam = "";

my $str_table1 = "";
my $str_table2 = "";

if(defined $annotation_files){ # If annotation files are entered (for vFAMs or viralOGs models)
    my %families;
    my %hmms = ();
    foreach my $line (@lines){  	
    	if ($line =~ m/^#/){
     	    next;
        }
     	else {	    
   	    my @aux = split(" ", $line);
	    if(defined $vfams{$aux[2]}){
		if($resp == 0 || $resp == 1 || $resp == 2){
                    if(defined $usr_score){
			if($aux[5] >= $usr_score) {
			    my @num_fam = split("/",$vfams{$aux[2]});
                            my @split = split("_", $num_fam[1]);
                            for(my $i = 0; $i < scalar(@split); ++$i){
                            	if($split[$i] eq "{}"){
                                    $split[$i] = "Empty";
                            	}
                            }
                            my $join = join(",", @split);
			    my $desc = "";
                            for(my $k = 18; $k < scalar(@aux); ++$k){
                            	$desc .= $aux[$k]." ";
                            }
                            $str_table1_vfam .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$num_fam[0]\t$join\t$desc\n";
                            foreach my $fam (@split){
                            	if(defined $families{$fam}){
                                    if(defined $families{$fam}->{$aux[2]}){
                                    	$families{$fam}->{$aux[2]}++;
                                    }
                                    else{
                                    	$families{$fam}->{$aux[2]} = 1;
                                    }
                            	}
                            	else{
                                    my %auxVfam;
                                    $auxVfam{$aux[2]} = 1;
                                    $families{$fam} = \%auxVfam;
                            	}
                            }   
                    	}
		    }
                    elsif(defined $usr_evalue) {
			my $calc_str = sprintf("%.10f", $aux[4]);			
                    	if( $calc_str <= $usr_evalue) {
                            my @num_fam = split("/",$vfams{$aux[2]});
                            my @split = split("_", $num_fam[1]);
                            for(my $i = 0; $i < scalar(@split); ++$i){
                            	if($split[$i] eq "{}"){
                                    $split[$i] = "Empty";
                            	}
                            }
                            my $join = join(",", @split);
			    my $desc = "";
                            for(my $k = 18; $k < scalar(@aux); ++$k){
                                $desc .= $aux[$k]." ";
                            }
                            $str_table1_vfam .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$num_fam[0]\t$join\t$desc\n";
                            foreach my $fam (@split){
                            	if(defined $families{$fam}){
                                    if(defined $families{$fam}->{$aux[2]}){
                                    	$families{$fam}->{$aux[2]}++;
                                    }
                                    else{
                                    	$families{$fam}->{$aux[2]} = 1;
                                    }
                            	}
                            	else{
                                    my %auxVfam;
                                    $auxVfam{$aux[2]} = 1;
                                    $families{$fam} = \%auxVfam;
                            	}
                            }
                    	}	
                    }
                    else{
			my @num_fam = split("/",$vfams{$aux[2]});
                    	my @split = split("_", $num_fam[1]);
                    	for(my $i = 0; $i < scalar(@split); ++$i){
                            if($split[$i] eq "{}"){
                            	$split[$i] = "Empty";
                            }
                    	}
                    	my $join = join(",", @split);
                    	my $calc_str = sprintf("%.6f", $aux[4]);
			my $desc = "";
                        for(my $k = 18; $k < scalar(@aux); ++$k){
                            $desc .= $aux[$k]." ";
                        }
                    	$str_table1_vfam .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$num_fam[0]\t$join\t$desc\n";
                    	foreach my $fam (@split){
                            if(defined $families{$fam}){
                            	if(defined $families{$fam}->{$aux[2]}){
                                    $families{$fam}->{$aux[2]}++;
                            	}
                            	else{
                                    $families{$fam}->{$aux[2]} = 1;
                            	}
                            }
                            else{
                            	my %auxVfam;
                            	$auxVfam{$aux[2]} = 1;
                            	$families{$fam} = \%auxVfam;
                            }
                    	}
                   }
		}
	    }
	    else{
		my @aux = split(" ", $line);
            	if($resp == 0 || $resp == 1){
                    if(defined $usr_score){
                    	if($aux[5] >= $usr_score) {
                            if(defined $hmms{$aux[2]}){
                            	++$hmms{$aux[2]};
                            }
                            else{
                            	$hmms{$aux[2]} = 1;
                            }
			    my @aux_tag = split("\t", $hash_att{$aux[2]});
			    my $desc = "";
                            for(my $k = 18; $k < scalar(@aux); ++$k){
                                $desc .= $aux[$k]." ";
                            }
                            $str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
                            if($taxon_tag > 0){
                                $str_table1 .= "\t$aux_tag[0]";
                            }
                            if($protein_tag > 0){
                                $str_table1 .= "\t$aux_tag[1]";
                            }
                            if($range_tag > 0){
                                $str_table1 .= "\t$aux_tag[2]";
                            }
                            if($type_tag > 0){
                                $str_table1 .= "\t$aux_tag[3]";
                            }
                            $str_table1 .= "\n";
                    	}
                    }
                    elsif(defined $usr_evalue) {
                    	my $calc_str = sprintf("%.10f", $aux[4]);
                    	if( $calc_str <= $usr_evalue) {
                            if(defined $hmms{$aux[2]}){
                            	++$hmms{$aux[2]};
                            }
                            else{
                            	$hmms{$aux[2]} = 1;
                            }
			    my @aux_tag = split("\t", $hash_att{$aux[2]});
			    my $desc = "";
                            for(my $k = 18; $k < scalar(@aux); ++$k){
                                $desc .= $aux[$k]." ";
                            }
                            $str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
                            if($taxon_tag > 0){
                                $str_table1 .= "\t$aux_tag[0]";
                            }
                            if($protein_tag > 0){
                                $str_table1 .= "\t$aux_tag[1]";
                            }
                            if($range_tag > 0){
                                $str_table1 .= "\t$aux_tag[2]";
                            }
                            if($type_tag > 0){
                                $str_table1 .= "\t$aux_tag[3]";
                            }
                            $str_table1 .= "\n";
                    	}
                    }
                    else{
                    	if(defined $hmms{$aux[2]}){
                            ++$hmms{$aux[2]};
                    	}
                    	else{
                           $hmms{$aux[2]} = 1;
                    	}
			my @aux_tag = split("\t", $hash_att{$aux[2]});
			my $desc = "";
                        for(my $k = 18; $k < scalar(@aux); ++$k){
                            $desc .= $aux[$k]." ";
                        }
                        $str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
                        if($taxon_tag > 0){
                            $str_table1 .= "\t$aux_tag[0]";
                        }
                        if($protein_tag > 0){
                            $str_table1 .= "\t$aux_tag[1]";
                        }
                        if($range_tag > 0){
                            $str_table1 .= "\t$aux_tag[2]";
                        }
                        if($type_tag > 0){
                            $str_table1 .= "\t$aux_tag[3]";
                        }
                        $str_table1 .= "\n";
                    }
            	}
		elsif($resp == 2){
                    my $name = lc($aux[2]);
                    $name =~ s/\s+//g;
                    if(defined $hmms_score{$name}){
                    	my $aux_sc = $hmms_score{$name};
                    	if($aux[5] >= $aux_sc) {
                            if(defined $hmms{$aux[2]}){
                            	++$hmms{$aux[2]};
                            }
                            else{
                            	$hmms{$aux[2]} = 1;
                            }
			    my @aux_tag = split("\t", $hash_att{$aux[2]});
			    my $desc = "";
                            for(my $k = 18; $k < scalar(@aux); ++$k){
                                $desc .= $aux[$k]." ";
                            }
                            $str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
                            if($taxon_tag > 0){
                                $str_table1 .= "\t$aux_tag[0]";
                            }
                            if($protein_tag > 0){
                                $str_table1 .= "\t$aux_tag[1]";
                            }
                            if($range_tag > 0){
                                $str_table1 .= "\t$aux_tag[2]";
                            }
                            if($type_tag > 0){
                                $str_table1 .= "\t$aux_tag[3]";
                            }
                            $str_table1 .= "\n";
                    	}
                    }
                    elsif(defined $usr_score){
                    	if($aux[5] >= $usr_score) {
                            if(defined $hmms{$aux[2]}){
                            	++$hmms{$aux[2]};
                            }
                            else{
                            	$hmms{$aux[2]} = 1;
                            }
			    my @aux_tag = split("\t", $hash_att{$aux[2]});
			    my $desc = "";
                            for(my $k = 18; $k < scalar(@aux); ++$k){
                                $desc .= $aux[$k]." ";
                            }
                            $str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
                            if($taxon_tag > 0){
                            	$str_table1 .= "\t$aux_tag[0]";
                            }
                            if($protein_tag > 0){
                            	$str_table1 .= "\t$aux_tag[1]";
                            }
                            if($range_tag > 0){
                            	$str_table1 .= "\t$aux_tag[2]";
                            }
                            if($type_tag > 0){
                            	$str_table1 .= "\t$aux_tag[3]";
                            }
                            $str_table1 .= "\n";
                    	}
                    }
                    elsif(defined $usr_evalue){
                    	my $calc_str = sprintf("%.10f", $aux[4]);
                    	if( $calc_str <= $usr_evalue) {
                            if(defined $hmms{$aux[2]}){
                            	++$hmms{$aux[2]};
                            }
                            else{
                            	$hmms{$aux[2]} = 1;
                            }
			    my @aux_tag = split("\t", $hash_att{$aux[2]});
			    my $desc = "";
                            for(my $k = 18; $k < scalar(@aux); ++$k){
                                $desc .= $aux[$k]." ";
                            }
                            $str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
                            if($taxon_tag > 0){
                            	$str_table1 .= "\t$aux_tag[0]";
                            }
                            if($protein_tag > 0){
                            	$str_table1 .= "\t$aux_tag[1]";
                            }
                            if($range_tag > 0){
                            	$str_table1 .= "\t$aux_tag[2]";
                            }
                            if($type_tag > 0){
                            	$str_table1 .= "\t$aux_tag[3]";
                            }
                            $str_table1 .= "\n";
                    	}
                    }
		    else{
                    	if(defined $hmms{$aux[2]}){
                            ++$hmms{$aux[2]};
                    	}
                    	else{
                            $hmms{$aux[2]} = 1;
                    	}
			my @aux_tag = split("\t", $hash_att{$aux[2]});
			my $desc = "";
                        for(my $k = 18; $k < scalar(@aux); ++$k){
                            $desc .= $aux[$k]." ";
                        }
                        $str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
                        if($taxon_tag > 0){
                            $str_table1 .= "\t$aux_tag[0]";
                        }
                        if($protein_tag > 0){
                            $str_table1 .= "\t$aux_tag[1]";
                        }
                        if($range_tag > 0){
                            $str_table1 .= "\t$aux_tag[2]";
                        }
                        if($type_tag > 0){
                            $str_table1 .= "\t$aux_tag[3]";
                        }
                        $str_table1 .= "\n";
                    }
		}
	    }
    	}
    }
    foreach my $fam (sort keys %families){
    	my $count = 0;
    	my $string = "";
    	my %vFams = %{$families{$fam}};
    	my $nVfam = 0;
    	foreach my $vf (sort keys %vFams){
	     $count += $vFams{$vf};
	     $string .= $vf."(".$vFams{$vf}.") ";
	     ++$nVfam;
	     $str_table3_vfam .= "$vf\t$fam\t$vFams{$vf}\n";
    	}
        $fam =~ s/ //i;
        $fam =~ s/family//i;
	$str_table2_vfam .= "$fam\t$count\t$nVfam\t$string\n";
    }
    foreach my $key (sort keys %hmms){
	$str_table2 .= "$key\t$hmms{$key}";
	my @aux_tag = split("\t", $hash_att{$key});
        if($taxon_tag > 0){
            $str_table2 .= "\t$aux_tag[0]";
        }
        if($protein_tag > 0){
            $str_table2 .= "\t$aux_tag[1]";
        }
        if($range_tag > 0){
            $str_table2 .= "\t$aux_tag[2]";
        }
        if($type_tag > 0){
            $str_table2 .= "\t$aux_tag[3]";
        }
        $str_table2 .= "\n";
	print $fl "Profile HMM $key does not found in vFam database. Generating quantitative reports without annotation data.\n";
    }
}
else{ # If no annotation files are entered.
    print STDERR "Annotation data not specified. Generating quantitative reports without annotation data\n";
    print $fl "Annotation data not specified. Generating quantitative reports without annotation data\n";
    foreach my $line (@lines){
    	if ($line =~ m/^#/){
            next;
     	}
     	else {
            my @aux = split(" ", $line);
	    if($resp == 0 || $resp == 1){
                if(defined $usr_score){
		    if($aux[5] >= $usr_score) {
                    	if(defined $hmms{$aux[2]}){
                            ++$hmms{$aux[2]};
                    	}
                    	else{
                            $hmms{$aux[2]} = 1;
                    	}
			my $desc = "";
			for(my $k = 18; $k < scalar(@aux); ++$k){
			    $desc .= $aux[$k]." ";
			}
			$str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
			if(defined $hash_att{$aux[2]}){
			    my @aux_tag = split("\t", $hash_att{$aux[2]});
			    if($taxon_tag > 0){
			    	$str_table1 .= "\t$aux_tag[0]";
			    }
			    if($protein_tag > 0){
                            	$str_table1 .= "\t$aux_tag[1]";
                            }
			    if($range_tag > 0){
                           	$str_table1 .= "\t$aux_tag[2]";
                            }
			    if($type_tag > 0){
                            	$str_table1 .= "\t$aux_tag[3]";
                            }
			}
		        $str_table1 .= "\n";
                    }
                }
		elsif(defined $usr_evalue) {
                    my $calc_str = sprintf("%.10f", $aux[4]);
                    if( $calc_str <= $usr_evalue) {
                    	if(defined $hmms{$aux[2]}){
                            ++$hmms{$aux[2]};
                    	}
                    	else{
                            $hmms{$aux[2]} = 1;
                    	}
			my $desc = "";
                        for(my $k = 18; $k < scalar(@aux); ++$k){
                            $desc .= $aux[$k]." ";
                        }
                        $str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
			if(defined $hash_att{$aux[2]}){
                            my @aux_tag = split("\t", $hash_att{$aux[2]});
                            if($taxon_tag > 0){
                            	$str_table1 .= "\t$aux_tag[0]";
                            }
                            if($protein_tag > 0){
                            	$str_table1 .= "\t$aux_tag[1]";
                            }
                            if($range_tag > 0){
                         	$str_table1 .= "\t$aux_tag[2]";
                            }
                            if($type_tag > 0){
                            	$str_table1 .= "\t$aux_tag[3]";
                            }
			}
                        $str_table1 .= "\n";
                    }
            	}
            	else{
                    if(defined $hmms{$aux[2]}){
                    	++$hmms{$aux[2]};
                    }
                    else{
                    	$hmms{$aux[2]} = 1;
                    }
		    my $desc = "";
                    for(my $k = 18; $k < scalar(@aux); ++$k){
                        $desc .= $aux[$k]." ";
                    }
                    $str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
		    if(defined $hash_att{$aux[2]}){
                        my @aux_tag = split("\t", $hash_att{$aux[2]});
                        if($taxon_tag > 0){
                    	    $str_table1 .= "\t$aux_tag[0]";
                    	}
                    	if($protein_tag > 0){
                            $str_table1 .= "\t$aux_tag[1]";
                    	}
                    	if($range_tag > 0){
                            $str_table1 .= "\t$aux_tag[2]";
                    	}
                    	if($type_tag > 0){
                            $str_table1 .= "\t$aux_tag[3]";
                    	}
		    }
                    $str_table1 .= "\n";
            	}
	    }
	    elsif($resp == 2){
		my $name = lc($aux[2]);
                $name =~ s/\s+//g;
                if(defined $hmms_score{$name}){
		    my $aux_sc = $hmms_score{$name};
		    if($aux[5] >= $aux_sc) {
                    	if(defined $hmms{$aux[2]}){
                            ++$hmms{$aux[2]};
                    	}
                    	else{
                            $hmms{$aux[2]} = 1;
                    	}
			my $desc = "";
                        for(my $k = 18; $k < scalar(@aux); ++$k){
                            $desc .= $aux[$k]." ";
                        }
			$str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
			if(defined $hash_att{$aux[2]}){
			    my @aux_tag = split("\t", $hash_att{$aux[2]});			
                            if($taxon_tag > 0){
                            	$str_table1 .= "\t$aux_tag[0]";
                            }
                            if($protein_tag > 0){
                            	$str_table1 .= "\t$aux_tag[1]";
                            }
                            if($range_tag > 0){
                            	$str_table1 .= "\t$aux_tag[2]";
                            }
                            if($type_tag > 0){
                            	$str_table1 .= "\t$aux_tag[3]";
                            }
			}
                        $str_table1 .= "\n";
                    }
		}
                elsif(defined $usr_score){
                    if($aux[5] >= $usr_score) {
                        if(defined $hmms{$aux[2]}){
                            ++$hmms{$aux[2]};
                        }
                        else{
                            $hmms{$aux[2]} = 1;
                        }
			my $desc = "";
                        for(my $k = 18; $k < scalar(@aux); ++$k){
                            $desc .= $aux[$k]." ";
                        }
                        $str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
			if(defined $hash_att{$aux[2]}){
                            my @aux_tag = split("\t", $hash_att{$aux[2]});
                            if($taxon_tag > 0){
                            	$str_table1 .= "\t$aux_tag[0]";
                            }
                            if($protein_tag > 0){
                            	$str_table1 .= "\t$aux_tag[1]";
                            }
                            if($range_tag > 0){
                            	$str_table1 .= "\t$aux_tag[2]";
                            }
                            if($type_tag > 0){
                           	 $str_table1 .= "\t$aux_tag[3]";
                            }
			}
                        $str_table1 .= "\n";
                    }
		}
                elsif(defined $usr_evalue){
	 	    my $calc_str = sprintf("%.10f", $aux[4]);
                    if( $calc_str <= $usr_evalue) {
                    	if(defined $hmms{$aux[2]}){
                            ++$hmms{$aux[2]};
                    	}
                    	else{
                            $hmms{$aux[2]} = 1;
                    	}
			my $desc = "";
                        for(my $k = 18; $k < scalar(@aux); ++$k){
                            $desc .= $aux[$k]." ";
                        }
                        $str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
			if(defined $hash_att{$aux[2]}){
                            my @aux_tag = split("\t", $hash_att{$aux[2]});
                            if($taxon_tag > 0){
                            	$str_table1 .= "\t$aux_tag[0]";
                            }
                            if($protein_tag > 0){
                            	$str_table1 .= "\t$aux_tag[1]";
                            }
                            if($range_tag > 0){
                            	$str_table1 .= "\t$aux_tag[2]";
                            }
                            if($type_tag > 0){
                            	$str_table1 .= "\t$aux_tag[3]";
                            }
			}
                        $str_table1 .= "\n";
                    }
                }
		else{
		    if(defined $hmms{$aux[2]}){
                        ++$hmms{$aux[2]};
                    }
                    else{
                        $hmms{$aux[2]} = 1;
                    }
		    my $desc = "";
                    for(my $k = 18; $k < scalar(@aux); ++$k){
                        $desc .= $aux[$k]." ";
                    }
                    $str_table1 .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\t$desc";
		    if(defined $hash_att{$aux[2]}){
                    	my @aux_tag = split("\t", $hash_att{$aux[2]});
                    	if($taxon_tag > 0){
                    	    $str_table1 .= "\t$aux_tag[0]";
                    	}
                    	if($protein_tag > 0){
                           $str_table1 .= "\t$aux_tag[1]";
                    	}
                    	if($range_tag > 0){
                            $str_table1 .= "\t$aux_tag[2]";
                    	}
                    	if($type_tag > 0){
                            $str_table1 .= "\t$aux_tag[3]";
                    	}
		    }
                    $str_table1 .= "\n";
		}
	    }
	}
    }
    foreach my $key (sort keys %hmms){
	$str_table2 .= "$key\t$hmms{$key}";
	if(defined $hash_att{$key}){
            my @aux_tag = split("\t", $hash_att{$key});
            if($taxon_tag > 0){
            	$str_table2 .= "\t$aux_tag[0]";
            }
            if($protein_tag > 0){
            	$str_table2 .= "\t$aux_tag[1]";
            }
            if($range_tag > 0){
            	$str_table2 .= "\t$aux_tag[2]";
            }
            if($type_tag > 0){
            	$str_table2 .= "\t$aux_tag[3]";
            }
	}
        $str_table2 .= "\n";
    }    
}

# Printing the results in output files.
if(!$str_table1_vfam eq ""){
    my $table_1_vfam = "$output/vFam_table1.csv";
    open(TBL1, ">$table_1_vfam") or die "ERROR: Could not create file $table_1_vfam!\n";
    print TBL1 "Target Name\tquery_pHMM\tE-value\tScore\t# of seqs\tFamily\tDescription of target\n";
    print TBL1 $str_table1_vfam;
    close(TBL1);
}

if(!$str_table2_vfam eq ""){
    my $table_2_vfam = "$output/vFam_table2.csv";
    open(TBL2, ">$table_2_vfam") or die "ERROR: Could not create file $table_2_vfam!\n";
    print TBL2 "Family\t# of seqs\t# of vFams\tvFams\n";
    print TBL2 $str_table2_vfam;
    close(TBL2);    
}

if(!$str_table3_vfam eq ""){
    my $table_3_vfam = "$output/vFam_table3.csv";
    open(TBL3, ">$table_3_vfam") or die "ERROR: Could not create file $table_3_vfam!\n";
    print TBL3 "vFam\tFamily\t# of seqs\n";
    print TBL3 $str_table3_vfam;
    close(TBL3);
}

if(!$str_table1 eq ""){
    my $table_1 = "$output/table1.csv";
    open(TBL1, ">$table_1") or die "ERROR: Could not create file $table_1!\n";
    print TBL1 "Target Name\tquery_pHMM\tE-value\tScore\tDescription of target";
    if($taxon_tag > 0){
	print TBL1 "\tTaxon";
    }
    if($protein_tag > 0){
	print TBL1 "\tProtein";
    }
    if($range_tag > 0){
	print TBL1 "\tRange";
    }
    if($type_tag > 0){
	print TBL1 "\tType";
    }
    print TBL1 "\n";
    print TBL1 $str_table1;
    close(TBL1);
}

if(!$str_table2 eq ""){
    my $table_2 = "$output/table2.csv";
    open(TBL2, ">$table_2") or die "ERROR: Could not create file $table_2!\n";
    print TBL2 "HMM\t# of seqs";
    if($taxon_tag > 0){
        print TBL2 "\tTaxon";
    }
    if($protein_tag > 0){
        print TBL2 "\tProtein";
    }
    if($range_tag > 0){
        print TBL2 "\tRange";
    }
    if($type_tag > 0){
        print TBL2 "\tType";
    }
    print TBL2 "\n";
    print TBL2 $str_table2;
    close(TBL2);
}

close ($filehadlefromtabular);

print "Done.\n";
print $fl "Done.\n";
close($fl);

exit;

################################################################
###                       Subroutines                        ###
################################################################

# This routine checks if profile HMMs contain cutoff scores values to be used in the 
# analysis.
sub verifyCutoff{
    my $hmm = shift;
    my $value = `grep SCORE $hmm`;
    $value =~ s/CUTOFF SCORE//;
    $value =~ s/\s//g;
    $value =~ s/\://g;
    if((defined $value) and !($value eq "")){
        return $value;
    }
    else{
        return undef;
    }
}

# This routine checks if a cutoff score or e-value was entered by the user.
sub verifyScoreEvalue{
    my $hmm = shift;
    my $aux_score = verifyCutoff($hmm);
    my $resp = 0;
    if(defined $aux_score){
	$replace =~ s/\s//g;
	if(lc($replace) eq "no"){	    
	    my $number = `grep -c NAME $hmm`;
	    if($number == 1){
		$usr_score = $aux_score;
            	$usr_evalue = undef;
		if(defined $read_length){
                    my $model_length = `grep LENG $hmm`;
                    $model_length =~ s/LENG//;
                    $model_length =~ s/\s//g;
                    $model_length =~ s/\://g;
                    my $read_length_prot = ($read_length/3);
                    if($read_length_prot < $model_length){
                        print $fl "Model: $hmm\n";
                        print $fl "Cutoff score: $usr_score\n";
                        print $fl "Read Length: $read_length\n";
                        my $new_value = ($read_length_prot/$model_length)*$usr_score;
                        $usr_score = $new_value;
                        print $fl "Adjusted cutoff score: $usr_score\n";
                    }
                }
		else{
            	    print $fl "Cutoff score: $usr_score\n";
		}
		return 1;
	    } 
	    else{
		my $dir = $output."/hmms";
		system "mkdir $dir";
		$/ = "//\n";
    		open(DATA, "<$hmm");
    		foreach my $unique_hmm (<DATA>){
        	    my @hmm_lines = split("\n", $unique_hmm);
		    my $aux_n;
		    my $achei = 0;
        	    foreach my $line (@hmm_lines){
			if($line =~ m/(NAME\d*)/g){
			    $aux_n = $line;
			    $aux_n =~ s/NAME//;
			    $aux_n =~ s/\s+//g;
			}
			elsif($line =~ m/(CUTOFF SCORE\d*)/g){
			    my $value = $line;
                            $value =~ s/CUTOFF SCORE//;
                            $value =~ s/\s//g;
                            $value =~ s/\://g;
			    my $file = $dir."/".$aux_n.".hmm";
			    open(HMMFILE, ">$file") or die "ERROR: Couldn't create file $output/$file : $!\n";
                            print HMMFILE $unique_hmm;
                            close (HMMFILE);
			    if(defined $read_length){
                                my $model_length = `grep LENG $file`;
				$model_length =~ s/LENG//;
				$model_length =~ s/\s//g;
				$model_length =~ s/\://g;
				my $read_length_prot = ($read_length/3);
            			if($read_length_prot < $model_length){
				    print $fl "Model: $file\n";
                                    print $fl "Cutoff score: $value\n";
                                    print $fl "Read Length: $read_length\n";
                		    my $new_value = ($read_length_prot/$model_length)*$value;
                		    $value = $new_value;
				    print $fl "Adjusted cutoff score: $value\n";
            			}
                            }
                            $hmms_score{lc($aux_n)} = $value;
			    $achei = 1;			    
                            next;
			}
        	   }
		   if($achei == 0){
			my $file = $dir."/".$aux_n.".hmm";
			open(HMMFILE, ">$file") or die "ERROR: Couldn't create file $output/$file : $!\n";
                        print HMMFILE $unique_hmm;
                        close (HMMFILE);
		   }
              }
    	    close (DATA);
	    $/ = "\n";
	    return 2;
           }
	}
	elsif(lc($replace) eq "yes"){
	    if(defined $usr_evalue){
                print $fl "Cutoff evalue: $usr_evalue\n";
            }
            elsif(defined $usr_score){
             	print $fl "Cutoff score: $usr_score\n";
            }
	    return 0;
	}	
    }
    else{
	if(($usr_score) and ($usr_evalue)) {
	    system "rm -rf $output";
	    die "ERROR: Please use -s or -e, not both.\n$help_print\n";
        }
	if(defined $usr_evalue){
            print $fl "Cutoff evalue: $usr_evalue\n";
        }
        elsif(defined $usr_score){
            print $fl "Cutoff score: $usr_score\n";
        }
	return 0;
    }
}

# This routine indetifies input file type.
sub verifyDatabaseType{
    my $file = shift;
    unless( open(DATA, "<$file") ) {
    	die "ERROR: Could not open database file $file: $!\n\n";
    }

    my $fasta_qual;
    my $ext = 'fasta';

    while (my $line = <DATA>) {
    	# check if DB is FASTA
    	if ($line =~ /^>/) {
            $ext = "fasta";
            $fasta_qual = $file . ".qual" if (-e "$file\.qual");
            last;
    	}
    	# check if DB is FASTQ
    	elsif($line =~ /^@/) {
            $ext = "fastq";
            last;
    	}
    	else {
            last;
    	}
    }
    return $ext;
}

# This routine checks if a directory with the user-specified output name already exists. 
# If so, a numeric suffix is added to the name.
sub output_dir_name {
    my $output_dir_name = shift;
    my $count = 2;
    my $flag = 0;
    print STDERR "Creating output directory $output_dir_name\n";
    if (-d $output_dir_name) {
        $flag = 1;
        while (-d "$output_dir_name\_$count") {
            $count++;
        }
        $output_dir_name = "$output_dir_name\_$count";
    }
    print STDERR "\nOutput directory already exists, saving results to $output_dir_name instead.\n\n" if $flag;
    return ($output_dir_name);
}

# This routine identifies if the fasta file contains nucleotide or protein sequences.
sub verifyFastaFileComposition{
    my $file = shift;
    my $count = 0;
    my $num_lines = 0;
    my $type = 1; # 1 - nucleotide; 2 - protein
    open(my $vf, "$file");
    while(<$vf>){
    	chomp($_);
	if($_ =~ /^>/){}
        else{
	    if ((index($_, 'L') != -1) or (index($_, 'l') != -1)) {
	    	++$count;
	    }
	    if ((index($_, 'I') != -1) or (index($_, 'i') != -1)) {
            	++$count;
            }
	    if ((index($_, 'P') != -1) or (index($_, 'p') != -1)) {
            	++$count;
            }
	    if ((index($_, 'F') != -1) or (index($_, 'f') != -1)) {
            	++$count;
            }
	    if ((index($_, 'Q') != -1) or (index($_, 'q') != -1)) {
            	++$count;
            }
	    if ((index($_, 'E') != -1) or (index($_, 'e') != -1)) {
            	++$count;
            }
	    if($count > 0){
	    	$type = 2;
	    	last;
	    }
	    else{
		++$num_lines;
		if($num_lines > 20){
		    last;
		}
	    }
        }
    }
       
    close($file);
    return $type;
}

# This routine checks if the user-specified tabular input file is a hmmsearch tabular 
# output file.
sub verifyTabularFile{
    my $file = shift;
    open(my $tab, "$file");
    while(<$tab>){
        chomp($_);
	if($_ =~ /target name/){
	    my @line = split(" ", $_);
	    if(($line[3] eq "accession") and ($line[4] eq "query") and ($line[6] eq "accession") and ($line[7] eq "E-value")){
		close($file);
		return 1;
	    }
	}
    }
    close($file);
    return 0;
}

# This routine runs transeq program.
sub runTranseq{
    my $name = shift;
    my $file = shift;
    my $dir = shift;
    my $string;
    my $aux;    
    if(defined $dir){
    	$string = $dir."/".$name."*_transeq*";
	$aux = `ls $string 2> /dev/null`;
	if($aux eq ""){
	    $string = $name."*_transeq*";
	    $aux = `ls $string 2> /dev/null`;
	}
    }
    else{
	$string = $name."*_transeq*";
	$aux = `ls $string 2> /dev/null`;	
    }

    if($aux eq ""){	
	$aux = undef;
    }
    else{
	my @aux_t = split("\n", $aux);
	$aux = $aux_t[0];
    }
    my $transeq;
    my $auxiliar;
    print $fl "Step: translation of nucleic acid sequences\n";
    if(defined $aux){ #Verify if transeq file exists
	print $fl "Transeq protein file found. Skipping translation of nucleic acid sequences\n";
	print "Transeq protein file found. Skipping translation of nucleic acid sequences\n";
	$transeq = $aux;
    }
    else{ #run transeq program:
     	$transeq = $name."_transeq.fasta";
	print $fl "Performing conceptual translation of $file and creating $transeq (Parameters: frame = 6) \n";
        print "Performing conceptual translation of $file and creating $transeq... \n";
        system "transeq -frame 6 $file $transeq";
        print "Done.\n";
	print $fl "Done\n";
    }
    return $transeq;
}

# This routine runs hmmsearch program.
sub runHmmsearch{
    my $name = shift;
    my $transeq = shift;
    my $option = shift;
    my $file1 = "$output/".$name."_hmmsearch.tab";
    my $full = "$output/".$name."_hmmsearch.txt";
    print $fl "Step: similarity search with hmmsearch.\n";
    if(-e $file1){
	my $type = verifyTabularFile($file1);
        if($type == 1){
	    print $fl "Similarity search file found. Skipping similarity search with hmmsearch.\n";
            print "Similarity search file found. Skipping similarity search with hmmsearch.\n";
            system "mv $file1 $output";
	}
	else{
	    print $fl "Step: similarity search with hmmsearch.\n";
	    die "Similarity search file found, but it is not in a tabular format.\n";
	}
	
    }
    else{
    	if(!$hmm_db){
	    print $fl "ERROR: Missing mandatory argument -i.\nTry -h for more details.\n";
            print "ERROR: Missing mandatory argument -i.\nTry -h for more details.\n";
            die "\n";
        }
        else{
	    if($option == 0 || $option == 1){
		if(!defined $usr_evalue and !defined $usr_score){
		    print $fl "Performing similarity search with hmmsearch (Parameters: -E 10 --tblout $output/$file1 -o $output/$full --cpu $cpu $hmm_db $transeq)\n";
                    print "Performing similarity search with hmmsearch... \n";
                    !system "hmmsearch -E 10 --tblout $file1 -o $full --cpu $cpu $hmm_db $transeq 2>> $log" or
                        die "ERROR: Could not run hmmsearch : $!\nCommand: hmmsearch --tblout $output/$file1 -o $output/$full $hmm_db $transeq  2>> $log\n";
		}
		else{
	    	    print $fl "Performing similarity search with hmmsearch (Parameters: --tblout $output/$file1 -o $output/$full --cpu $cpu $hmm_db $transeq)\n";
            	    print "Performing similarity search with hmmsearch... \n";
	    	    !system "hmmsearch -T 1 --tblout $file1 -o $full --cpu $cpu $hmm_db $transeq 2>> $log" or
                    	die "ERROR: Could not run hmmsearch : $!\nCommand: hmmsearch --tblout $output/$file1 -o $output/$full $hmm_db $transeq  2>> $log\n";
		}
	    	print "Done\n";
	    	print $fl "Done\n";
	    }
	    else{
		my $full = "$output/".$name."_hmmsearch.txt";
		my $aux_full = "$output/full.txt";
		my $aux_tab = "$output/tab.txt";
		open(FULL, ">$aux_full");
		open(TAB, ">$aux_tab");
		my $dir = $output."/hmms";		
		opendir(DIR, "$dir");
        	my @hmms = readdir(DIR);
        	closedir(DIR);
        	my $n = 0;
		my $count = 0;
		print $fl "Performing similarity search with hmmsearch (Parameters: --tblout $output/$file1 -o $output/$full --cpu $cpu $hmm_db $transeq)\n";
                print "Performing similarity search with hmmsearch... \n";
        	foreach my $hmm (@hmms){
            	    if($hmm eq "." or $hmm eq ".."){}
            	    else{
			my $file_full = $output."/".$hmm."_full";
		  	my $file_tab = $output."/".$hmm."_tab";
			my $name = $hmm;
			$name =~ s/\.hmm//ig;
			if(defined $hmms_score{$name}){
			    !system "hmmsearch -T 10 --tblout $file_tab -o $file_full --cpu $cpu $dir/$hmm $transeq 2>> $log" or
                                die "ERROR: Could not run hmmsearch : $!\nCommand: hmmsearch --tblout $file_tab -o $file_full $hmm $transeq  2>> $log\n";
			}
			else{
        	            !system "hmmsearch -E 10 --tblout $file_tab -o $file_full --cpu $cpu $dir/$hmm $transeq 2>> $log" or
                	    	die "ERROR: Could not run hmmsearch : $!\nCommand: hmmsearch --tblout $file_tab -o $file_full $hmm $transeq  2>> $log\n";				      }
			open(FILE, $file_full);
			while(<FILE>){
			    chomp($_);
			    if($_ =~ /#/){
				if($count == 0){
				    print FULL "$_\n";
				}
			    }
			    else{
				print FULL "$_\n";
			    }
			}
			close(FILE);
			system "rm $file_full";
			open(FILE, $file_tab);
                        while(<FILE>){
                            chomp($_);
                            if($_ =~ /#/){
				if($count == 0 and ($_ =~ /full sequence/ or $_ =~ /target name/)){
				    print TAB "$_\n";
				}
			    }
                            else{
                                print TAB "$_\n";
                            }
                        }
                        close(FILE);
			system "rm $file_tab";
			++$count;
            	    }
        	}
		close(FULL);
		close(TAB);
		system "mv $aux_tab $file1";
            	system "mv $aux_full $full";
		system "rm -rf $dir";
	    }
        }
    }
    return $file1;
}

# This routine checks the input file type (fastq, fasta or tabular file).
sub verify_file_type{ # 1 - fasta file; 2 - fastq file; 3 - tabular file
    my $file = shift;
    open(FILE, $file);    
    my $found_fastq = 0;
    my $found_fasta = 0;
    my $count = 0;
    my $count_fasta = 0;
    my $count_h = 0;
    my $seq = "";
    while(<FILE>){
        chomp($_);
        if($_ =~ /^>/){
	    if(($found_fasta == 1) and ($seq ne "")){
		++$count_fasta
	    }
	    $found_fasta = 1;
	    if($count_fasta > 0){
            	close(FILE);
            	return 1;
	    }
        }
	elsif($_ =~ /^@/){
	    $found_fastq = 1;
	}
	elsif(($_ =~ /^\+/) and ($found_fastq == 1)){
	    close(FILE);
            return 2;
	}
	else{
	    if($found_fasta == 1){
		$seq .= $_;
	    }
	}
	++$count;
	if($count > 100){
	    last;
	}
    }
    if(($found_fasta == 1) and ($seq ne "")){
   	close(FILE);
        return 1;
    }
    close(FILE);
    return 3;
}
