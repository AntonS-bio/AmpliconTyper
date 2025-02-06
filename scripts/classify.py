import argparse
import subprocess
from typing import List, Dict
import pickle
from os.path import  expanduser, join, exists, basename
from os import mkdir
from shutil import which, rmtree
import uuid
import asyncio
from check_for_update import UpdateChecker
from data_classes import Amplicon
from input_processing import InputProcessing
from model_manager import ModelManager
from classifier_report import ClasssifierReport
from read_classifier import Classifier, ReadsMatrix, ModelsData, ClassificationResult
import warnings
warnings.filterwarnings("ignore")
from map import ReadMapper
import pysam as ps

VERSION="0.1.28"

def generate_amplicons(model_file: str, fasta_file:str) -> List[Amplicon]:
    target_regions: List[Amplicon] = []
    with open(model_file, "rb") as input_model:
        models_data: ModelsData = pickle.load(input_model)
        model_manager: Dict[str, Classifier] = models_data.classifiers
    for name, model in model_manager.items():
        bed_line=f'{model.name}\t0\t{len(model.nucletoide_seq)}\t{model.name}'
        target_regions.append(Amplicon.from_bed_line(bed_line, fasta_file))
    return target_regions

def generate_ref_fasta(model_file: str, ref_file: str) -> None:
    """Output reference FASTA from model file. This FASTA will be used for 
    mapping of reads and will be removed afterwards

    :param models_file: file name with saved model
    :type models_file: str
    :param ref_file: file name to which to save FASTA sequences
    :type models_file: str
    """
    with open(model_file, "rb") as input_model:
        models_data: ModelsData = pickle.load(input_model)
        model_manager: Dict[str, Classifier] = models_data.classifiers
    with open(ref_file, "w" ) as ref_output:
        for name, model in model_manager.items():
            ref_output.write(">"+str(model.name)+"\n")
            ref_output.write(str(model.nucletoide_seq)+"\n")

def get_postive_reads(results: List[ClassificationResult], target_reads_bams: str):
    """Output bam files (one per sample) with reads classified as target organism
    
    :param results: List of classification results 
    :type results: List[ClassificationResult]
    :param target_reads_bams: Directory to which bam files will be saved
    :type target_reads_bams: str
    """
    postive_read_ids: Dict[str, List[str]]={}
    for result in results:
        if result.sample not in postive_read_ids:
            postive_read_ids[result.sample]=[]
        for read_id in result.positive_read_ids:
            postive_read_ids[result.sample].append(read_id)

    for sample, ids in postive_read_ids.items():
        bam = ps.AlignmentFile(sample, "rb",check_sq=False)
        output_file=join(target_reads_bams,basename(sample))
        target_reads_bam=ps.AlignmentFile( output_file, "wb", template=bam)
        for read in bam.fetch():
            if read.query_name in ids:
                target_reads_bam.write(read)
        target_reads_bam.close()
        subprocess.run(f'samtools index {output_file} {output_file+".bai"}', shell=True, executable="/bin/bash", stdout=subprocess.PIPE,  stderr=subprocess.PIPE)

def get_arguments():
    parser = argparse.ArgumentParser(description='Classify reads in BAM file using existing model or train a model from bam files')
    parser.add_argument('-f','--fastqs', type=str,
                        help='Directory with ONT run results or individual FASTQ file', required=False)
    parser.add_argument('-b','--bams',  type=str,
                        help='Directory with, list of, or individual BAM and corresponding BAM index files (.bai)', required=True)
    parser.add_argument('-m','--model', type=str,
                        help='Pickle (.pkl) file containing pretrained model. Model must be trained on same reference', required=True)
    parser.add_argument('-d','--sample_descriptions',  type=str,
                        help='File with sample descriptions (tab delimited), first column must be the BAM file name without .bam', required=False)
    parser.add_argument('-c','--description_column',  type=str,
                        help='Column in sample description file to use to augmnet samples descriptions', required=False)
    parser.add_argument('--column_separator',  type=str,
                        help='The column separator for the sample descriptions file', required=False, default="\t")
    parser.add_argument('--cpus', type=int, help='Number of CPU cores to use', required=False, default=1)
    parser.add_argument('-o','--output_file', type=str,
                        help='File to store classification results', required=True)
    parser.add_argument('--target_reads_bams', type=str,
                        help='Directory to which to write reads classified as target organism. Using this will substantially increase run time.', required=False)
    parser.add_argument('-s', '--max_soft_clip', type=int,
                        help='Specifies maximum permitted soft-clip (length of read before first and last mapped bases) of mapped read', required=False, default=25)
    parser.add_argument('-l', '--max_read_len_diff', type=int,
                        help='Specifies maximum nucleotide difference between mapped portion of read and target amplicon', required=False, default=10)
    parser.add_argument('--organism_presence_cutoff', type=float,
                        help='Currently not in use', required=False, default=20.0)
   
    parser.add_argument('-v', '--version', action='version', version='%(prog)s '+VERSION)

    try:
        args = parser.parse_args()
    except Exception as e:
        parser.print_help()
        return

    return args

async def classify(temp_dir, args):

    mkdir(temp_dir)

    if args.organism_presence_cutoff<1 and args.organism_presence_cutoff!=0:
        print(f'Option "organism_presence_cutoff" has value {args.organism_presence_cutoff} which is below 1, but not 0. Valid values are 0 to 100.')
        return

    ReadsMatrix.permitted_read_soft_clip=int(args.max_soft_clip)
    ReadsMatrix.permitted_mapped_sequence_len_mismatch=int(args.max_read_len_diff)

    #check minimap2 is installed
    if which("minimap2") is None:
        print(f'Missing minimap2 program. It is available via Bioconda. Exiting.')
        return

    cpu_to_use=int(args.cpus)

    input_processing=InputProcessing()

    model_file = args.model
    if not exists(model_file):
        print(f'Could not located model file {model_file}. Exiting.')
        return

    fasta_file = join(temp_dir,"reference.fasta")
    generate_ref_fasta(model_file, fasta_file)
    generate_amplicons(model_file, fasta_file)

    if not args.sample_descriptions is None: #this has to be checked early in case file does not exist
        metadata_file=expanduser(args.sample_descriptions)
        if not input_processing.file_exists(metadata_file):
            return
        metadata_column, metadata_sep = input_processing.check_metadata(args, metadata_file)
        if not input_processing.check_column_in_descriptions(metadata_file,metadata_column, metadata_sep):
            return

    if not args.target_reads_bams is None:
        if not exists(args.target_reads_bams):
            mkdir(args.target_reads_bams)
        target_reads_bams=expanduser(args.target_reads_bams)


    #Map the raw reads or collect bams to classify
    bams_dir=expanduser(args.bams)
    if not args.fastqs is None:
        args.fastqs=expanduser(args.fastqs)
        if not input_processing.file_exists(args.fastqs):
            return
        print("Starting read mapping to reference")
        mapper=ReadMapper(temp_dir,
                        args.fastqs,
                        model_file, fasta_file, cpu_to_use)
        
        #check that sample labels are present in metadatafile
        sample_labels: Dict[str, str] = {}
        if not  args.sample_descriptions is None:
            files_to_check=mapper.output_names(bams_dir)
            if not input_processing.get_sample_labels(metadata_file,metadata_sep, metadata_column, files_to_check, sample_labels):
                return #terminate if some sample is not present in metadata.

        if not mapper.map(bams_dir):
            return
        file_to_classify=[f.bam_file for f in mapper.results]
    else: #use the bams provided
        file_to_classify=input_processing.get_bam_files( args.bams )
        sample_labels: Dict[str, str] = {}
        if not  args.sample_descriptions is None:
            mapper=ReadMapper(temp_dir,
                        bams_dir,
                        model_file, fasta_file, cpu_to_use)
            if not input_processing.get_sample_labels(metadata_file,metadata_sep, metadata_column, file_to_classify, sample_labels):
                return #terminate if some sample is not present in metadata.

    if len(file_to_classify)==0:
        print("No files to classify. Please check that you have specified correct FASTQs directory if using -f or existing BAMs if not using -f")
        return
    

    output_file=input_processing.check_output_dir(args.output_file)

    target_regions: List[Amplicon] = generate_amplicons(model_file, fasta_file)

    model_manager=ModelManager(model_file, cpu_to_use)
    results=model_manager.classify_new_data(target_regions, file_to_classify)

    if not args.target_reads_bams is None:
        get_postive_reads(results, target_reads_bams)
        generate_ref_fasta(model_file, join(target_reads_bams,"reference.fasta"))

    print("Generating report")
    length_parameters={"-s": ReadsMatrix.permitted_read_soft_clip, "-l":ReadsMatrix.permitted_mapped_sequence_len_mismatch}
    if not args.fastqs is None:
        report=ClasssifierReport(output_file, model_file, args.organism_presence_cutoff, length_parameters, sample_labels, mapping_results=mapper.results)
    else:
        report=ClasssifierReport(output_file, model_file, args.organism_presence_cutoff, length_parameters, sample_labels)
    report.create_report(results)
    rmtree(temp_dir)
    print("Done")
    
async def run():
    temp_dir=expanduser( join("./",str(uuid.uuid4())) )

    try:
        input_arguments=get_arguments()
        checker=UpdateChecker(input_arguments.model)
        await asyncio.gather(checker.get_result(VERSION), classify(temp_dir, input_arguments))
        print("\n"+checker.result)

    except Exception as e:
       print(e)
    if exists(temp_dir): #clean up
        rmtree(temp_dir)

def main():
    pass

if __name__=="classify" or __name__=="__main__":
    loop = asyncio.get_event_loop()
    loop.run_until_complete(run())