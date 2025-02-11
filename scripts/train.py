from os.path import expanduser
import traceback
from sys import exc_info
import argparse
from typing import List
from inputs_validation import ValidateFiles
from data_classes import Amplicon
from model_manager import ModelManager
from read_classifier import ReadsMatrix
from input_processing import InputProcessing

def main():

    parser = argparse.ArgumentParser(description='Classify reads in BAM file using existing model or train a model from bam files')
    parser.add_argument('-t','--target_regions',  type=str,
                        help='Bed  file that specifies reference intervals to which reads where mapped', required=True)
    parser.add_argument('-r','--reference', type=str,
                        help='FASTA file with the reference to which reads where mapped', required=True)
    parser.add_argument('-p','--positive_bams', type=str,
                        help='Directory with or list of NON-target BAM files and corresponding BAM index files (.bai)', required=False)
    parser.add_argument('-n','--negative_bams', type=str,
                        help='Directory with or list of TARGET BAM files and corresponding BAM index files (.bai)', required=False)
    parser.add_argument('-v','--vcf', type=str,
                        help='VCF file containing positions that will be excluded from training as well as genotype defining SNPs (also excluded)', required=True)
    parser.add_argument('-s','--special_cases', type=str,
                        help='Tab delimited file specifying amplicon for which presence/absense should be reported', required=False)
    parser.add_argument('-o','--output_dir', type=str,
                        help='Directory for output files', required=True)
    parser.add_argument('--cpus', type=int,
                        help='Directory for output files', required=False, default=1)
    parser.add_argument('-m', '--bams_matrix',  type=str, default=None,
                        help='Matrix with precalculated BAM matrices ()', required=False)

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        exit(0)

    cpu_to_use=int(args.cpus)

    input_processing=InputProcessing()
    if args.bams_matrix is None:
        positive_bams=input_processing.get_bam_files( args.positive_bams )
        print(positive_bams)
        negative_bams=input_processing.get_bam_files( args.negative_bams )
        if len(positive_bams)==0 or len(negative_bams)==0:
            exit()

    bed_file = expanduser(args.target_regions)
    if not input_processing.file_exists(bed_file):
        exit(0)

    fasta_file = expanduser(args.reference)
    if not input_processing.file_exists(fasta_file):
        exit(0)

    vcf_file=expanduser(args.vcf)
    if not input_processing.file_exists(vcf_file):
        exit(0)

    output_dir=args.output_dir
    input_processing.make_dir(output_dir)
    model_file=f'{output_dir}/model.pkl'
    model_quality_file=f'{output_dir}/quality.tsv'

    use_existing_bams_matrix=False
    if not args.bams_matrix is None:
        use_existing_bams_matrix=True
        matrix_file=expanduser(args.bams_matrix)
        if not input_processing.file_exists(matrix_file):
            exit(0)
    else:
        matrix_file=f'{output_dir}/bams_matrix.pkl'


    file_validator=ValidateFiles()
    file_validator.validate_bed(bed_file)
    file_validator.validate_fasta(fasta_file)
    file_validator.contigs_in_fasta(bed_file, fasta_file )

    target_regions: List[Amplicon] = []
    with open(bed_file) as input_file:
        for line in input_file:
            target_regions.append(Amplicon.from_bed_line(line, fasta_file))


    model_manager=ModelManager(model_file, cpu_to_use)
    model_manager.model_evaluation_file=model_quality_file
    model_manager.load_genotype_snps(vcf_file)
    ReadsMatrix.permitted_read_soft_clip=2000
    ReadsMatrix.permitted_mapped_sequence_len_mismatch=5
    if use_existing_bams_matrix:
        model_manager.train_from_existing_matrices(target_regions=target_regions, matrices_file=matrix_file)
    else:
        model_manager.train_from_bams(target_regions, positive_bams, negative_bams, matrix_file)



if __name__=="__main__":
    try:
        main()
    except Exception as e:
        exc_type, exc_obj, exc_tb = exc_info()
        lf = traceback.extract_tb(exc_tb)[-1]
        print("Program terminated with an error:")
        print("File: ", lf.filename,", Line: ", lf.lineno)
        print(e)
    
