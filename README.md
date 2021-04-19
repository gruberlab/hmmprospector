# HMM-prospector: A tool for surveying genomic and metagenomic data with profile HMMs

HMM-Prospector is a Perl script that uses a single or multiple profile HMM file as a query in similarity searches against a FASTQ/FASTA dataset using hmmsearch program. HMM-Prospector processes the results and generates tabular files with qualitative and quantitative results. Previously run hmmsearch result files (short tabular format) can also be used as datasets.

## Instalation

HMM-Prospector also does not need to be installed. The user should only download the hmm-prospector.pl file.

## Requirements

- hmmsearch: HMMER3package-http://hmmer.org/. 

- transeq: EMBOSS package - http://emboss.sourceforge.net/ .

- fastq_to_fasta: FASTX-Toolkit–http://hannonlab.cshl.edu/fastx_toolkit/.

P.S: All programs must be located in a directory listed in the PATH of the operating system.

## Usage

```
perl hmm-prospector.pl -d <file> -s|-e <decimal>  <optional parameters>
```  

### Mandatory parameters:

```
-d  <file>        : Dataset (FASTQ, FASTA or hmmsearch's tabular output file)
```

### Optional parameters:

```
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

-v|version        : Version.
```
## Tutorial
Follow the instructions in the HMM-Prospector Manual file to learn how to use HMM-Prospector program and interpret the results.

## Reference
If you use this program for your publication, please cite:

Oliveira, L.S. and Gruber, A. (2021) Rational design of profile HMMs for viral classification and discovery. In Nakaya, H. (Ed.), Bioinformatics. Exon Publications, Brisbane, Australia, pp. 151-170.

## Contact

To report bugs, to ask for help and to give any feedback, please contact Arthur Gruber (argruber@usp.br) or Liliane S. Oliveira (liliane.sntn@gmail.com).

