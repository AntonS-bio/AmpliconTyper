# AmpliconTyper
Tool for genotyping Oxford Nanopore amplicon sequencing data

## Purpose
This tool is intended to be used with a partner [AmpliSeqDesigner](https://github.com/AntonS-bio/AmpliSeqDesigner) which helps to design primers for environmental surveillance. The partner tools designs primers for specific genotypes (based on SNPs unique to those genotypes) while minimising risk of amplyfing other DNA that might be present in the environment.

This tool, AmpliconTyper, has two modes: 
    1. Train a classification model using existing Nanopore data - this can be your own data, data from genomic repositories (ENA, NCBI, DDBJ) or a mix of these
    2. Classify ONT sequencing data using the model trained in (1)



## Setup
The easiest way to setup the tool is to use conda or mamba. If you are using macOS or Linux, you can install via command line. If you are using Windows 10+, the best option is to setup Windows Subsystem Linux (WSL) first which will give you access to Linux functionality from you Windows system. After that, you can use conda or mamba.

The best practice is to install packages in dedicated environment to avoid software conflicts. To create new environment and install AmpliconTyper into it use:
```
conda create --name  amplicontyperENV -c bioconda -c conda-forge amplicontyper
```
Once installed, use
```
conda activate amplicontyperENV
```
to activate the environment. Now you are ready to use the AmpliconTyper. You would need to run activation command (but not create command) every time you start a new terminal. In both commands above "amplicontyperENV" can be replaced with whatever you want to call the environment. 




## Model training
If you have been given a model pickle file - you don't have to do this, so go to Read Classification section. If you haven't been given one, ask if anyone on your project already has it - chances are bioinformatician does. There are some models in this repository in "models" directory. These are versioned and maintained. You need to be sure you are using the right model for you set of amplicons.

The idea of training is that to distinguish target and non-target organism, AmpliconTyper needs to see a lot of examples of amplicon sequencing data from both. The examples of target organisms should be fairly easy to find - your own sequencing data should work. For non-target organism, you probably won't be able to generate sufficient data yourself, so you have to rely on public data. SRA database at NCBI is a good source. I recommend finding taxa for your target organism/serovar/strain and taking all available ONT sequencing data within two-three taxonomic layers above. For example, if you are looking at Salmonella, you may decide to take all ONT data from all taxa within Enterobacterales. Once raw ONT data is mapped to FASTQs these will be your negative BAMs. If in your test sequencing runs you often see an off-target amplification of same organism, add it to a negative set.

Once you have mapped the the sequencing data you've amassed to reference amplicons, you can train your model. The training process is quick, but generating BAM files for it can take several days. This is not related to AmpliconTyper, but to volume of data you need to process. However, training only needs to be done once and the "train" function itself below should take fewer than 10 minutes to run. It will generate a model file that you can share with colleagues to avoid them having to train the model.

```
usage: train -t  -r  -p  -v  -o  [-m] [-h] [-n] 

Classify reads in BAM file using existing model or train a model from bam files. You will get an error if you use "train.py" instead of "train"

options:
  -h, --help            show this help message and exit
  -t , --target_regions 
                        Bed file that specifies reference intervals to which reads where mapped
  -r , --reference      FASTA file with the reference to which reads where mapped
  -p , --positive_bams 
                        Directory with or list of NON-target BAM files and corresponding BAM index files (.bai)
  -n , --negative_bams 
                        Directory with or list of TARGET BAM files and corresponding BAM index files (.bai)
  -v , --vcf            VCF file containing positions that will be excluded from training as well as genotype defining SNPs (also excluded)
  -s,  --special_cases  Tab delimited file specifying amplicon for which presence/absense should be reported
  -o , --output_dir     Directory for output files
  -m , --bams_matrix    Matrix with precalculated BAM matrices

```

## Reads classification (must have a trained model, see Model training above)
The purpose of classifier is to take all reads provided to it and identify which came from your target and which from non-target. The classifier will use extra information you provide to generate a nice HTML report file with results - you can see sample here in test_data/report.html. If you download this file it can be opened in any browser. 

**IMPORTANT all BAMs must be indexed using samtools index command.**

To run classifier you have to main options. First, if you have already mapped FASTQs to amplicons reference sequence, you can provide the report with either directory with BAM files, a single BAM file or a plain text file with list of BAM files (full address including directory) to classify. Examples below use supplied test data from test_data.


### I already mapped FASTQs to references and generated and indexed BAM files

This will classify all BAMs in directory ./bams/ 
```
classify -b ./bams/ -m model.pkl -d metadata.tsv -c Description -o report.html
```
This will only classify sample_1.bam
```
classify -b ./bams/sample_1.bam -m model.pkl -d metadata.tsv -c Description -o report.html
```




### I only have FASTQ files, but no BAMs

If you don't have bams, you can use "-f" option to map FASTQs using minimap2 to the amplicon reference. The command only differs from above by this extra parameter. Not that "-b" is still present - it will be a directory to which the mapped BAMs will be saved. 
```
classify -b ./bams/ -f ./fastqs/ -m model.pkl -d metadata.tsv -c Description -o report.html
```
If you look at structure of test_data/fastqs you will see it's a bit odd. It's very difficult to capture all possible ways the FASTQs will be provided, but the most likely one if a result of single run from ONT devices. When working with multiplex sequencing run (most likely use case), an ONT device will create a directory for each barcode and will place multiple FASTQ files for this barcode into that directory. It will likely look like this:
```
Run_20
      ↳ barcode1
               ↳ FAX83461_pass_barcode1_503f1897_ffa628d1_1.fastq.gz
               ↳ FAX83461_pass_barcode1_503f1897_ffa628d1_2.fastq.gz
               ↳ FAX83461_pass_barcode1_503f1897_ffa628d1_3.fastq.gz
    
      ↳ barcode2
               ↳ FAX83461_pass_barcode2_503f1897_ffa628d1_1.fastq.gz
               ↳ FAX83461_pass_barcode2_503f1897_ffa628d1_2.fastq.gz
               ↳ FAX83461_pass_barcode2_503f1897_ffa628d1_3.fastq.gz
    
      ↳ barcode3
               ↳ FAX83461_pass_barcode3_503f1897_ffa628d1_1.fastq.gz
               ↳ FAX83461_pass_barcode3_503f1897_ffa628d1_2.fastq.gz
               ↳ FAX83461_pass_barcode3_503f1897_ffa628d1_3.fastq.gz
```
Normally, you'd need to merge the the files from each barcode before mapping, but AmpliconTyper will do it for you.

**IMPORTANT if AmpliconTyper is doing the mapping, it will call each output BAM by the name with corresponding barcode (ex. barcode2.bam), this may overwrite some old BAMs**

Neither the metadata file (-d), nor genotypes hierarchy (-g) are required, but you probably already have them and they substantially enrich the classification report, so it's worth using them. The metadata file is simply a delimited file that lists the samples in the first column, the genotypes hierarchy file was likely used when your PCR primers were generated. 

```
usage: classify -b  -m  -o  [-d] [-c] [--column_separator] [-g] [--cpu] [-h]

Classify reads in BAM file using existing model or train a model from bam files. You will get an error if you use "classify.py" instead of "classify"

options:
  -h, --help            show this help message and exit
  -b , --bam_files      Directory with, list of, or individual BAM and corresponding BAM index files (.bai)
  -m , --model          Pickle (.pkl) file containing pretrained model. Model must be trained on same reference
  -o , --output_file    Name of HTML file to store classification results
  -d , --sample_descriptions 
                        File with sample descriptions (tab delimited), first column must be the BAM file name without .bam
  -c , --description_column 
                        Column in sample descriptions file to use to augment samples descriptions
  --column_separator    The column separator for the sample descriptions file, default=tab
  --cpus                Number of CPUs to use, default=1
  --target_reads_bams   Directory to which to write reads classified as target organism
  -s, --max_soft_clip   Specifies maximum permitted soft-clip (length of read before first and last mapped bases) of mapped read.
  -l, --max_read_len_diff
                        Specifies maximum nucleotide difference between mapped portion of read and target amplicon.

```
Here is illustration of options -s and -l where "-" below is a nuclotide and " " is a gap in alignment.

![image](https://github.com/user-attachments/assets/bfb3fcc0-99b4-4ff9-bc64-9deb9614512c)

By default, reads with length different from reference (i.e. target amplicon sequence) are discarded. 
Option -s allows the read to overhang the reference on either side. With "-s 5" the Read_2 would be kept as it overhangs the reference by 5 nucleotides. Option "-s" would not have any effect on Read_3 and Read_4, they would be rejected. 
Options -l allows the aligned portion of read to be shorter or longer than reference lenth. This allows accomodation of deletion/insertion sequencing errors. With "-l 5" Read_3 and Read_4 will both be kept. Read_2 would also be kept because it's 5nt longer than reference.
**IMPORTANT While option -s is quite forgiving, -l can lead to problems as it does will not distiguish between sequencing errors and true InDel that may be a sign of different species**



### Understanding results

![GitHubReadMeImage](https://github.com/user-attachments/assets/9f86e6f1-ddad-4a86-8561-1d92c4aa9be3)






The summary of results consists of the following fields:

A) the name of the model file used and the date that model was trained.

B) A diagram of how many amplicons have a certain number of SNPs in consensus sequence. The position of the number indicates the number of SNPs and the number itself - the number of amplicons. In first sample, 3 amplicons have 0 SNPs (0th place), and 1 amplicon has 2 SNPs (third place). In the second sample, 12 amplicons have 0 SNPs (0th place) and 2 amplicons has 1 SNP (second place). This simultaneously shows the presence of products and identified possible off-target amplification.

C) In S. Typhi, some genotypes are labels are not strictly hierarchical. For example, genotypes 2, 3 and 4 are a subset of genotype 1, genotypes 3 and 4 are subset of genotype 2 and genotype 4 is subset of genotype 3. This means that identification of these genotypes cannot be done on basis of a single SNP, because SNP that separates genotypes 0-1 from genotype 2 will also be present in genotypes 3 and 4. For these high level genotypes special logic is used to determine which of these genotypes (0, 1, 2, 3 and 4) are possible with within a sample. The possible genotypes are reported in this field. When "Any" is reported, amplicons provide no information of high level genotypes.

D) List of genotype alleles identified in the amplicons. The number in brackets indicates percentage of reads (read depth) that support the genotype assignment. Sometimes, a reference position that is linked to specific genotype has unexpected nucleotide. This is reported as "Unknown allele at known position".

E) List of AMR linked alles idenfitied in the amplicons. The number in brackets indicates percentage of reads (read depth) that support the genotype assignment. In cases where genotype is based on presence/absence of amplicon (eg. S. Typhi PST6 plasmid) the amplicon name will be reported here if detected..

F) When classification was done from FASTQ files, this show the total number of reads in each sample.

G) This shows the total number of reads that mapped to some amplicon and also were of correct length taking options "-s" and "-l" into account. 

H) This is number of reads of wrong length. Note that normally, not all reads are mapped, so H + I <> G

The rest of report contains details on each sample as well as each sample's amplicon and these can be jumped to via hyperlinks in the report file.

In addtion to the HTML report, the classify command also outputs a tab delimited file with the same name as the HTML file. This delimited file contains all the fields from the HTML report for simpler aggregation of results and downstream analysis.


## Getting data out of model file

You may have noticed that the "classify" does not require FASTA as input. To avoid excessive inputs the reference sequence is stored in the model file. To access it, you can use genotyper_utilities function. 

For example, to print the reference sequences in file model.pkl you should run 

```
genotyper_utilities --reference -m model.pkl
```
to get the SNPs stored in the model.pkl you should run
```
genotyper_utilities --snps -m model.pkl
```
