from collections import Counter
from os.path import exists, basename, dirname, isdir
from multiprocessing import Pool
from typing import Dict, List, Set
import pickle
from tqdm import tqdm
import numpy as np
from read_classifier import ReadsMatrix, Classifier, ClassificationResult, number_dic, GenotypeSNP, ModelsData
from data_classes import Amplicon
import hashlib

class ModelManager:
    """Set of function to train a model from either BAM files or, in case
     BAMs have previously been converted to Numpy arrays, from file containing
     Numpy arrays

    :param models_file: file in/from which to save/load models
    :type models_file: List[Amplicon]     
    """
    def __init__(self, models_file: str, cpus: int) -> None:
        self.cpus_to_use=cpus
        self._model_output_file=models_file
        self._matrices_file=""
        self._target_regions: List[Amplicon]=[]
        self._excluded_samples:Set[str]=set()
        self._trained_models: Dict[str, Classifier] ={}
        self._genotype_snps:List[GenotypeSNP]=[]
        self._model_evaluation_file=""

    def train_from_bams(self, target_regions: List[Amplicon], positive_bams: List[str], negative_bams: List[str], matrices_output_file:str) -> None:
        """Creates nucleotides matrix from bams
        and trains model to differentiate reads from positive and negative bams

        :param target_regions: List of amplicons that to extract from BAMs
        :type target_regions: List[Amplicon]
        :param positive_bams: List of BAM files that correspond to positive cases
        :type positive_bams: List[str]
        :param negative_bams: List of BAM files that correspond to negative cases
        :type negative_bams: List[str]
        :param matrices_output_file: File where matrices would be pickled
        :type matrices_output_file: str
        """
        if not isdir(dirname(matrices_output_file)):
            raise ValueError(f'Terminating model training, directory {dirname(matrices_output_file)} does not exist')
        self._target_regions.clear()
        self._target_regions=[f for f in target_regions]
        self._matrices_file=matrices_output_file
        ReadsMatrix.train_mode=True
        postive_matrix=self._bams_to_matrix(positive_bams)
        negative_matrix=self._bams_to_matrix(negative_bams)
        with open(matrices_output_file, "wb") as output:
            result={"positive": postive_matrix, "negative": negative_matrix}
            pickle.dump(result, output)
            del postive_matrix
            del negative_matrix
        self._train_classifier(matrices_output_file)

    def train_from_existing_matrices(self, target_regions: List[Amplicon], matrices_file: str) -> None:
        """Creates nucleotides matrix from existing nucleotide matrices
        and trains model to differentiate reads from positive and negative matrices

        :param target_regions: List of amplicons that to extract from BAMs
        :type target_regions: List[Amplicon]
        :param matrices_file: Pickle file containing pickled matrices
        :type matrices_file: str
        """
        self._target_regions.clear()
        self._target_regions=[f for f in target_regions]
        self._train_classifier(matrices_file)

#region Properties
    @property
    def excluded_samples(self) -> Set[str]:
        """Samples to be excluded from model training
        Convinience function to exclude a few bams from training

        returns Set[str]
        """
        return self._excluded_samples

    def exclude_sample(self, value: str):
        """Add sample to list
        Convinience function to exclude a few bams from training

        :param value: Sample name (sampleXYZ.bam)
        :type value: str
        """
        self._excluded_samples.add(value)

    @property
    def model_evaluation_file(self) -> str:
        """File name which will contain model testing result

        returns str
        """
        return self._model_evaluation_file

    @model_evaluation_file.setter
    def model_evaluation_file(self, value: str):
        """File name which will contain model testing result

        :param value: File name 
        :type value: str
        """
        self._model_evaluation_file=value

    @property
    def trained_models(self) -> Dict[str, Classifier]:
        """Dictionary of trained models. Key is amplicon name.

        returns Dict[str, Classifier]
        """
        return self._trained_models
    
    @property
    def models_fingerprint(self) -> str:
        """Hash of the uuids of all models in this ModelManager object
        """
        hash256=hashlib.sha256()
        sorted_uuids=sorted([f.uuid for f in  self._trained_models ])
        uuids_concatenated="".join( [f for f in sorted_uuids] )
        hash256.update(uuids_concatenated.encode())
        return f'{hash256.hexdigest()[0:3]}-{hash256.hexdigest()[3:6]}-{hash256.hexdigest()[6:9]}-{hash256.hexdigest()[9:12]}' #This works because there are few models
        #and they don't need to be cryptographically secure.
       
    @trained_models.setter
    def trained_models(self, value: Dict[str, Classifier]):
        """Dictionary of trained models. Key is amplicon name.

        :param value: Dictionary containing trained model, key=amplicon name
        :type value: Dict[str, Classifier]
        """
        self._trained_models=value

#endregion
    
    def add_genotype_snp(self, vcf_line:str) -> None:
        """Converts a VCF line into GenotypeSNP object
        param vcf_line: SNP line from VCF file
        type vcf_line: str
        """
        vcf_values=vcf_line.strip().split("\t")
        vcf_values[1]=int(vcf_values[1])-1 #correct for beds and BAMs being 0 indexed and VCF 1 indexed
        for genotype_snp in self._genotype_snps:
            if genotype_snp.coordinates_match(vcf_values[0], vcf_values[1]):
                genotype_snp.add_from_vcf(vcf_line)
                return None #existing GenotypeSNP is updated, no need to continue
        genotype_snp: GenotypeSNP = GenotypeSNP(vcf_values[0], vcf_values[1], vcf_values[3])
        genotype_snp.add_from_vcf(vcf_line)
        self._genotype_snps.append(genotype_snp)

    def load_genotype_snps(self, vcf_file: str) -> None:
        if not exists(vcf_file):
            raise ValueError(f'VCF file {vcf_file} does not exist, did you specify just file name instead of full address?')
        with open(vcf_file) as input_file:
            for line in input_file:
                if line[0]=="#":
                    continue
                values=line.strip().split("\t")
                if len(values)>9:
                    self.add_genotype_snp(line)

    def _train_classifier(self, matrice_file: str) -> Classifier:
        """Trains ensemble classifier to distinguish 
        reads from positive and negative sets of BAMs
        
        :param matrice_file: Pickle file containing pickled matrices
        :type matrice_file: str
        """
        self.trained_models.clear()
        self.trained_models={}
        with open(matrice_file, "rb") as input_matrices:
            raw_data=pickle.load(input_matrices)
        
        for i, target in enumerate(self._target_regions):
            print(f'Starting classifier training {target.name} target {str(i+1)}/{str(len(self._target_regions))}')
            target_genotype_snps=[f for f in self._genotype_snps if f.contig_id==target.ref_contig]
            classifier=Classifier(target.name, target.seq, target_genotype_snps)
            #remove data which is typhi, by not classified as such by taxa
            for file in self.excluded_samples:
                if file in raw_data["negative"]:
                    raw_data["negative"].pop(file)
            negative_data=np.vstack( [ values[target.name] for key, values in raw_data["negative"].items() ] )
            negative_data_bams=np.vstack( [ np.asarray([key]*values[target.name].shape[0]).reshape(-1,1) for key, values in raw_data["negative"].items() ] )

            positive_data=np.vstack( [ values[target.name] for key, values in raw_data["positive"].items()] )

            #check that the nucleotides in the VCF with SNP and nucleotides in the positive matrix match

            print(f'Positive data shape: {positive_data.shape[0]} rows and {positive_data.shape[1]} columns')
            print(f'Negative data shape: {negative_data.shape[0]} rows and {negative_data.shape[1]} columns')

            columns_to_exclude=self._check_bam_vcf_consistency(positive_data, target.ref_contig)

            classifier.train_classifier(positive_data, negative_data, columns_to_exclude, negative_data_bams)
            del positive_data
            del negative_data
            self.trained_models[target.name]=classifier
        #Training is finished, clear the data
        del raw_data

        with open(self.model_evaluation_file,"w") as temp:
            for name, model in self.trained_models.items():
                temp.write(name+"\t"+str(model.sensitivity)+"\t"+str(model.specificity)+"\n")

        output = ModelsData()
        output.classifiers=self.trained_models
        output.metadata={}
        with open(self._model_output_file, "wb") as model_file:
            pickle.dump(output, model_file)

    def _check_bam_vcf_consistency(self, positive_data, ref_contig) -> List[int]:
        """Checks if the nucletides in VCF of variable positions 
        is consistent with the dominant nucleotides in BAM.
        Mismatch is not necessarily a problem, but lots of mismatches
        suggest a problem
        
        :param positive_data: Matrix of positive reads
        :type positive_data: np.array

        :return List of matrix positions to exclude
        :rtype: List[int]
        """
        columns_to_exclude=[f.position for f in self._genotype_snps if f.contig_id==ref_contig]
        if len(columns_to_exclude)>0:
            postive_matrix_nucleotides=np.apply_along_axis(lambda x: Counter(x).most_common(1)[0], 0, positive_data[:,columns_to_exclude])[0]
            postive_matrix_nucleotides=[number_dic[f] for f in postive_matrix_nucleotides]
            vcf_nucleotides=[f.reference_nucl for f in self._genotype_snps if f.contig_id==ref_contig]
            mismatches=[a==f for a, f in zip(vcf_nucleotides,postive_matrix_nucleotides)]
            if False in mismatches:
                print(f'VCF nucleotides do not match dominant nucleotide in positive BAMs, contig: {ref_contig}: {postive_matrix_nucleotides} vs {vcf_nucleotides}')
        return columns_to_exclude

    def _process_bams(self, bam_file: str,) -> ReadsMatrix:
        """Function for multiprocessesor processing of BAMs 
        Creates ReadMatrix object for specified 
        
        :param bam_file: Pickle file containing pickled matrices
        :type bam_file: str

        :return ReadsMatrix object with data from BAM file
        :rtype: ReadsMatrix
        """
        reads_matrix=ReadsMatrix(bam_file=bam_file, is_positive=False)
        reads_matrix.read_bam(self._target_regions)
        return reads_matrix

    def _bams_to_matrix(self, bam_files: List[str]) -> Dict[str, Dict[str,np.array]]:
        """Function for multiprocessesor processing of BAMs 
        Creates ReadMatrix object for specified 
        
        :param bam_files: Bam files to convert into numpy array
        :type bam_files: List[str]

        :return Dictionary with two numpy arrays corresponding to "positive" and "negative"
        BAMs
        :rtype: Dict[str, Dict[str,np.array]]
        """

        if __name__ == 'model_manager':
            print(self.cpus_to_use)
            with Pool(self.cpus_to_use) as p:
                msa_results = list(tqdm( p.imap(func=self._process_bams, iterable=bam_files), total=len(bam_files) ))
                bam_matrices:Dict[str, Dict[str,np.array]]=dict( (basename(f.bam_file), f.read_matrices) for f in msa_results  )
                return bam_matrices
            

    def _classify_new_data_helper(self, bam_file:str) -> ClassificationResult:
        '''Function called by multithreading pool to allow faster processing
        of data. The multithreading has to be based on simultaneous processing
        of amplicons, not BAMs due large size of the data inputs which prevents
        spawning classifiers for each BAM

        :param input_data: Tuple containing reads matrix, target amplicon and name of sample
        :type: Tuple[np.array, Amplicon, str]
        '''

        with open(self._model_output_file, "rb") as model_file:
            self.models_data: ModelsData = pickle.load(model_file)
            self.trained_models=self.models_data.classifiers

        reads_data: ReadsMatrix=self._process_bams(bam_file)

        results: List[ClassificationResult] = []
        for amplicon in self._target_regions:
            reads=reads_data.read_matrices[amplicon.name]

            if amplicon.name not in self.trained_models:
                raise ValueError(f'No models available for amplicon {amplicon.name}')
            
            relevant_model = self.trained_models[amplicon.name]
            classification: ClassificationResult=ClassificationResult(amplicon, bam_file)
            classification.read_ids = reads_data.read_ids[amplicon.name]
            classification.predicted_classes=relevant_model.classify_new(reads)
            classification.model_fingerprint=relevant_model.uuid()
            classification.model_timestamp=relevant_model.training_timestamp()
            classification._wrong_len_reads=reads_data.wrong_len_reads
            classification.get_consensus(reads[classification.predicted_classes==1])
            results.append(classification)
        del reads_data
        return results

    def classify_new_data(self, target_regions: List[Amplicon], bam_files: List[str]) -> List[ClassificationResult]:
        """Creates nucleotides matrix from existing nucleotide matrices
        and trains model to differentiate reads from positive and negative matrices

        :param target_regions: List of amplicons that to extract from BAMs
        :type target_regions: List[Amplicon]
        :param bam_file: BAM file with reads to classify
        :type bam_file: str
        """

        self._target_regions=target_regions
        print("Processing data")
        results: List[ClassificationResult] = []
        if __name__ == 'model_manager':
            with Pool(self.cpus_to_use) as p:
                predictions: List[ClassificationResult] = list(tqdm( p.imap(func=self._classify_new_data_helper, iterable=bam_files), total=len(bam_files) ))
                results = [k for f in predictions for k in f]
        return results

