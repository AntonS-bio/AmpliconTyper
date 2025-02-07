import pysam as ps
import numpy as np
from typing import List, Dict
from os.path import expanduser
import uuid
from Bio import SeqIO

base_dic={"A":1,"C":2,"G":3,"T":4,"N":5,"-":0}

train_mode=True
permitted_read_soft_clip=2000
permitted_mapped_sequence_len_mismatch=5

log=open("./debug.log","w")
def write_log(value: str):
    print(value)
    log.write(value+"\n")

class ReferenceSequence:
    def __init__(self, contig_id, seq_start, seq_end, sequence) -> None:
        self._refseq_id=contig_id
        self._ref_start=seq_start
        self._ref_end=seq_end
        self._sequence=sequence

    @classmethod
    def from_bed_line(cls, bed_line:str, ref_fasta_file: str):
        """Constructor using bedfile lines. 
        :param bed_line: String from bedfile, if the line has fourth column, this will be included in amplicon name
        :type bed_file: str

        :param ref_fasta_file: Path to fasta file on which the amplicon is based
        :type ref_fasta_file: str
        """
        bed_line_values = bed_line.strip().split("\t")
        contig_id, seq_start, seq_end=bed_line_values[0:3]
        seq_start=int(seq_start)
        seq_end=int(seq_end)
        with open(ref_fasta_file) as fasta_input:
            for record in SeqIO.parse(fasta_input,"fasta"):
                if record.id==contig_id:
                    if seq_end>len(str(record.seq)):
                        raise ValueError(f'Contig {contig_id} is shorter, {len(str(record.seq))}nt, than end position in the bedfile: {seq_end}')
                    new_refseq=cls(contig_id, seq_start, seq_end, str(record.seq[seq_start:seq_end]))
                    return new_refseq
        raise ValueError(f'Contig {contig_id} is not found in fasta file: {ref_fasta_file}')
    
    @property
    def refseq_id(self) -> str:
        return self._refseq_id

    @property
    def ref_start(self) -> int:
        return self._ref_start

    @property
    def ref_end(self) -> int:
        return self._ref_end
    
    @property
    def sequence(self) -> str:
        return self._sequence
    
class Amplicon:
    def __init__(self, name: str, seq: str) -> None:
        self._name: str=name
        self._seq: str=seq
        self._left_flanking_id=""
        self._right_flanking_id=""
        self._has_homologues=False
        self._uuid=str(uuid.uuid4())
        self._has_reference=False #indicates if the amplicon has associated reference sequence

    @classmethod
    def from_bed_line(cls, bed_line:str, ref_fasta_file: str):
        """Constructor using bedfile lines. 
        :param bed_line: String from bedfile, if the line has fourth column, this will be included in amplicon name
        :type bed_file: str

        :param ref_fasta_file: Path to fasta file on which the amplicon is based
        :type ref_fasta_file: str

        """
        bed_line_values = bed_line.strip().split("\t")
        if len(bed_line_values)>=4:
            name=bed_line_values[3]
        else:
            name='_'.join( [str(f) for f in bed_line_values[0:3] ] )

        new_refseq=ReferenceSequence.from_bed_line(bed_line,ref_fasta_file)
        new_amplicon=cls(name, new_refseq.sequence)
        new_amplicon.ref_seq=new_refseq
        return new_amplicon

    @property
    def ref_seq(self) -> ReferenceSequence:
        if not self._has_reference:
            raise ValueError(f'Amplicon {self._name} has no reference')
        return self._ref_seq

    @ref_seq.setter
    def ref_seq(self, value: ReferenceSequence):
        self._has_reference=True
        self._ref_seq=value

    @property
    def has_reference(self) -> bool:
        return self._has_reference

    @property
    def ref_contig(self) -> str:
        if not self._has_reference:
            raise ValueError(f'Amplicon {self._name} has no reference')
        return self._ref_seq.refseq_id

    @property
    def seq(self) -> str:
        if self._has_reference:
            return self.ref_seq.sequence
        else:
            return self._seq

    @property
    def name(self) -> str:
        return self._name

    @property
    def id(self) -> str:
        return self._uuid

    @property
    def len(self) -> int:
        return len(self.seq)

    @property
    def has_homologues(self) -> bool:
        return self._has_homologues

    @has_homologues.setter
    def has_homologues(self, value: bool):
        self._has_homologues = value

    def __hash__(self):
        return hash(self.id)

def read_bam(target_regions: List[Amplicon]) -> None:
    _read_matrices:  Dict[str, np.array]={}
    _wrong_len_reads:  Dict[str, int]={}
    bam = ps.AlignmentFile(bam_file, "rb",check_sq=False)
    for target_amplicon in target_regions:#amplicon_coordinates.index:
        write_log("processing amplicon "+target_amplicon.name )
        target_start, target_end = [target_amplicon.ref_seq.ref_start, target_amplicon.ref_seq.ref_end]
        target_contig=target_amplicon.ref_seq.refseq_id
        alignment_matrix=np.zeros( (bam.count(contig=target_contig, start=target_start, end=target_end), target_end-target_start), dtype=np.int8)
        
        used_rows=[False]* alignment_matrix.shape[0] #this will allow removal of empty rows when only full length reads are required
        read_ids=[""]* alignment_matrix.shape[0]
        orientation=[]
        _wrong_len_reads[target_amplicon.name]=0
        for i,read in enumerate(bam.fetch(contig=target_contig, start=target_start, end=target_end)):
            if not read.is_unmapped:
                if read.query_sequence is None:
                    write_log(f'Bam file {bam_file} has an empty read aligning to {target_amplicon.name} this should not happen')
                    continue
                #allows for soft clippping on both sides (first line), reads shorter than target (second line), reads longer than target (third line)
                if train_mode and not (read.query_alignment_start < permitted_read_soft_clip \
                                                and (read.query_length-read.query_alignment_end) < permitted_read_soft_clip \
                                                and (read.query_alignment_end - read.query_alignment_start) < (target_end - target_start) + permitted_mapped_sequence_len_mismatch \
                                                and (read.query_alignment_end - read.query_alignment_start) > (target_end - target_start) - permitted_mapped_sequence_len_mismatch):
                    #read too short/too long
                    _wrong_len_reads[target_amplicon.name]+=1
                    continue
                orientation.append(read.is_forward)
                query_nt = [ read_ref_pair[0] for read_ref_pair in read.get_aligned_pairs() if not read_ref_pair[1] is None and
                                read_ref_pair[1]<target_end and read_ref_pair[1]>=target_start and not read_ref_pair[0] is None]
                ref_nt = [ read_ref_pair[1]-target_start for read_ref_pair in read.get_aligned_pairs() if not read_ref_pair[1] is None and
                            read_ref_pair[1]<target_end and read_ref_pair[1]>=target_start and not read_ref_pair[0] is None]
                bases=set([ read.query_sequence[f] for f in query_nt])
                if True in [f not in base_dic for f in bases]:
                    write_log(f'Ambigous bases not supported: {bam_file} target {target_amplicon.name}. Target will be ignored')
                    break

                alignment_matrix[i, ref_nt]=[ base_dic[read.query_sequence[f]] for f in query_nt]
                used_rows[i] = True
                read_ids[i] = read.query_name
        write_log("pysam indicates "+str(alignment_matrix.shape[0])+" reads in this "+bam_file+ " for amplicon "+target_amplicon.name+\
        " of these  "+str(_wrong_len_reads[target_amplicon.name])+" were rejected due to lengths")
        _read_matrices[target_amplicon.name]=alignment_matrix[used_rows]
    return _read_matrices

write_log("using pysam version " + ps.__version__)
bed_file=expanduser("./amplicons.bed")
fasta_file=expanduser("./amplicons.fna")
#Load the BED file into list of AmpliconObjects
target_regions: List[Amplicon] = []
with open(bed_file) as input_file:
    for line in input_file:
        target_regions.append(Amplicon.from_bed_line(line, fasta_file))



bams=['./positive_bams/ERR13868196.bam', './positive_bams/ERR13868190.bam', './positive_bams/ERR13868195.bam']
for bam in bams:
    write_log("reading bam "+bam)
    bam_file=expanduser(bam)
    results=read_bam(target_regions)
    for item, matrix in results.items():
        write_log("read_bam returned "+str(matrix.shape)+" matrix for amplicon "+item )
log.close()