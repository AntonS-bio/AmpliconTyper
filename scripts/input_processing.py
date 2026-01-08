from typing import List, Tuple, Dict
from os.path import expanduser, isdir, exists, dirname, join, split
from os import listdir, mkdir
from pathlib import Path

class InputProcessing:

    def __init__(self) -> None:
        pass

    def get_bam_files(self, bam_source: str) -> List[str]:
        bam_source=expanduser(bam_source)
        file_to_classify=[]
        if isdir(bam_source):  #Source is directory with BAM files
            file_to_classify=[join(bam_source,f) for f in listdir(bam_source)  if
                            f.split(".")[-1]=="bam"]
            if len(file_to_classify)==0:
                print(f'No files ending with .bam found in {bam_source}')
        elif bam_source.split(".")[-1]=="bam": #Source is single BAM file
            if exists(bam_source):
                file_to_classify.append(bam_source)
            else:
                print(f'BAM file {bam_source} does not exist. Maybe the location is incorrect?')
        elif not exists(bam_source):
            mkdir(bam_source)
        else:
            #check that list of bams contains valid bams
            with open(bam_source) as bam_list:
                for line in bam_list:
                    tentative_bam=expanduser(line.strip())
                    if exists(tentative_bam):
                        file_to_classify.append(tentative_bam)
                    else:
                        print(f'BAM file {tentative_bam} does not exist. Maybe the location is incorrect?')
        return file_to_classify
    
    def check_column_in_descriptions(self, description_file: str, column: str, separator: str) -> bool:
        with open(description_file) as input_file:
            first_line=input_file.readline().strip().split(separator)
            if column in first_line:
                return True
            else:
                print(f'Column {column} not found in file {description_file} using column separator "{separator}"')
                return False

    def get_sample_labels(self, description_file: str, separator: str, label_column: str, bam_files: List[str], sample_labels: Dict[str, str]) -> bool:
        first_column_values:Dict[str, str]={}
        with open(description_file) as input_file:
            first_line = input_file.readline().strip().split(separator)
            column_index = first_line.index(label_column)
            for line in input_file:
                if len(line.strip())==0: #The line is empty, likely last line of file.
                    continue
                first_column_values[line.strip().split(separator)[0]]=line.strip().split(separator)[column_index]

        outcome=True
        for file in bam_files:
            if Path(file).stem not in first_column_values:
                print(f'BAM derived sample name "{Path(file).stem}" not found in first column of metadata file.')
                outcome=False
            else:
                sample_labels[Path(file).stem]=first_column_values[Path(file).stem]
        return outcome

    def get_ref_bed_model(self, args) -> Tuple[str, str, str]:
        bed_file=expanduser(args.target_regions)
        fasta_file=expanduser(args.reference)
        model_file=expanduser(args.model)
        for file, name in zip( [bed_file, fasta_file, model_file], ["BED", "FASTA", "Model"]):
            if not exists(file):
                print(f'{name} file {file} does not exist.')
                return ("", "", "")
        return (bed_file, fasta_file, model_file)

    def file_exists(self, file_name: str) -> bool:
        if not exists(file_name):
            if isdir(file_name):
                print(f'Directory {file_name} does not exist. Check the address spelling.')
                return False
            else:
                print(f'File {file_name} does not exist. Check the address spelling.')
                return False
        return True

    def make_dir(self, dir_name: str) -> bool:
        if not exists(dir_name):
            mkdir(dir_name)
        else:
            #print(f'Directory {dir_name} already exists.')
            return False
        return True

    
    def check_output_dir(self, address:str) -> str:
        full_address=expanduser(address)
        dir, file = split(full_address)
        if dir=="":
            #only file is supplied, output the results in current directory
            return expanduser(join(".",address))
        else:
            if not exists(dir):
                print(f'File {address} requires directory {dir} which does not exist. Please create and try again.')
                raise ValueError()
            else:
                return full_address
            

    def check_metadata(self, args, metadata_file: str) -> Tuple[str, str]:
        if args.description_column is None and not args.sample_descriptions is None:
            print("Description file was provided (-d), but the column to use was not (-c)")        
            raise ValueError()
        if not args.description_column is None and args.sample_descriptions is None:
            print("Description column was provided (-c), but the file was not (-d)")
            raise ValueError()
        elif not args.description_column is None and not args.sample_descriptions is None:
            metadata_column = args.description_column
            metadata_sep = args.column_separator

            return (metadata_column, args.column_separator)