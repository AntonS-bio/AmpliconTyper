import subprocess
import gzip
from typing import  Dict, List
from os import  walk,  remove, mkdir
from os.path import isfile, join,  exists, isdir, expanduser, split
import uuid
import pysam as ps
# from read_classifier import Classifier
# import pickle

class MappingResult:
    def __init__(self, bam_file:str) -> None:
        self._bam_file=bam_file
        self._total_reads=-1
        self._mapped_reads=-1
        self._mapped=False

    @property
    def description(self) -> str:
        if self._mapped:
            return f'Total reads: {self._total_reads} of which {self.mapped_reads} mapped. {self.bam_file}'
        else:
            return f'File {self._bam_file} was unmapped, something may have crashed'

    @property
    def sample_name(self) -> str:
        head, tail = split(self.bam_file)
        return tail.split(".")[0]


    @property
    def mapped(self) -> bool:
        return self._mapped

    @mapped.setter
    def mapped(self, value: bool):
        self._mapped = value

    @property
    def bam_file(self) -> str:
        return self._bam_file

    @bam_file.setter
    def bam_file(self, value: str):
        self._bam_file = value

    @property
    def total_reads(self) -> int:
        return self._total_reads

    @total_reads.setter
    def total_reads(self, value: int):
        self._total_reads = value

    @property
    def mapped_reads(self) -> int:
        return self._mapped_reads

    @mapped_reads.setter
    def mapped_reads(self, value: int):
        self._mapped_reads = value

    @property
    def unmapped_reads(self) -> int:
        return self._total_reads-self._mapped_reads

class ReadMapper:

    VALID_EXTENSION=["fq", "fastq"]

    def __init__(self, temp_dir:str, source: str, model_file: str, ref_fasta: str, cpus:int) -> None:
        self._source=expanduser(source)
        self.temp_dir=expanduser(temp_dir)
        if not exists(self.temp_dir):
            mkdir(self.temp_dir)
        self.model_file=expanduser(model_file)
        self.reference=ref_fasta
        if not exists(self.model_file):
            raise IOError(f'Model file {self.model_file} not found. Please check spelling')
        self._barcode_files: Dict[str, List[str]] = {}
        self.results: List[MappingResult] = []
        self._cpus=cpus

    def source_to_samples(self) -> None:
        self._barcode_files: Dict[str, List[str]]={} #key=path to barcode, value=list of FASTQ files
        if not exists(self._source):
            raise IOError(f'When mapping reads, FASTQ source file/dir {self._source} was not found. Please check spelling.')
        
        if isfile(self._source):
            self._barcode_files[split(self._source)[0]]=[split(self._source)[1]] #this preserves type continuity
        elif isdir(self._source):
            for path, subdirs, dir_files in walk(self._source):
                for name in dir_files:
                    name_elements=name.split(".")
                    if isfile(join(path, name)) and \
                        ( ( len(name_elements)>3 and name_elements[1] in ReadMapper.VALID_EXTENSION and name_elements[2]=="gz" ) or \
                          (  len(name_elements)==3  and name_elements[1] in ReadMapper.VALID_EXTENSION ) ):
                        if path not in self._barcode_files:
                            self._barcode_files[path] = []
                        self._barcode_files[path].append(name)

    def _merge_samples(self, source_dir: str) -> str:
        """Merges all files ending in .fq/.fastq with or without .gz into single file
        Nanopore usually outputs multiple files for same barcode

        :param source_dir: Source files directory, should match key from self._barcode_files
        :type source_dir: str
       
        :return GUID string value for the merged FASTQ
        :rtype: str
        """
        merged_file=join(self.temp_dir, str(uuid.uuid4())+".fq.gz" )
        buffer_size = 16000
        with open(merged_file, 'wb') as dest_file:
            for filename in self._barcode_files[source_dir]:
                with open( join(source_dir,filename), 'rb') as source_file:
                    chunk = True
                    while chunk:
                        chunk = source_file.read(buffer_size)
                        dest_file.write(chunk)
        return merged_file

    def _map_helper(self, fastq: str,  output_name:str) -> MappingResult:
        """Maps supplied FASTQ to reference and returns total and mapped read count
        
        :param fastq: either file ending in .fq/fq.gz or folder with FASTQ files
        :type fastq: str
        
        :param reference: Reference FASTA sequence for mapping
        :type reference: str
        
        :param output_name: Full name to output BAM file
        :type output_name: str

        :return Tuple with total and mapped reads count
        :rtype: Tuple(int,int)
        """
        result=MappingResult(output_name)
        result.total_reads = sum(1 for f in gzip.open(fastq, 'rb')) / 4
        subprocess.run(f'minimap2 -ax map-ont {self.reference} {fastq} | samtools view -@ {self._cpus} -F 2308 -bS  | samtools sort -@ {self._cpus} -o {output_name}', \
                shell=True, executable="/bin/bash", stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
        # 2308=256(not primary) + 4(not aligned)+2064(not supplementary alignment)
        subprocess.run(f'samtools index {output_name}', shell=True, executable="/bin/bash", stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
        bam = ps.AlignmentFile(output_name, "rb",check_sq=False)
        result.mapped_reads = bam.mapped
        bam.close()
        result.mapped = True
        return result

    def map(self, bams_dir: str) -> bool:
        if not exists(bams_dir):
            mkdir(bams_dir)
        if not isdir(bams_dir):
            print("You specified --fastqs, so --bams must be a directory into which mapped files will be placed.")
            return False

        self.source_to_samples()
        self.results: List[MappingResult]=[]
        for key in self._barcode_files.keys():
            merged_file=self._merge_samples(key)
            self.results.append( self._map_helper(merged_file,  join(bams_dir, key.split("/")[-1]+".bam") ) )
            print(self.results[-1].description)
            remove(merged_file)
        return True

    def output_names(self, bams_dir: str) -> bool:
        if not exists(bams_dir):
            mkdir(bams_dir)
        if not isdir(bams_dir):
            print("You specified --fastqs, so --bams must be a directory into which mapped files will be placed.")
            return False

        result=[]
        self.source_to_samples()
        self.results: List[MappingResult]=[]
        for key in self._barcode_files.keys():
            result.append(join(bams_dir, key.split("/")[-1]+".bam") )
        return result

