from collections import Counter
from os.path import exists
import warnings
from typing import Dict, List
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.mixture import GaussianMixture
# from sklearn.neighbors import KNeighborsClassifier
# from sklearn.naive_bayes import GaussianNB, CategoricalNB
# from sklearn.neural_network import MLPClassifier
# from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.linear_model import SGDClassifier
import numpy as np
import pysam as ps
from data_classes import Amplicon




base_dic={"A":1,"C":2,"G":3,"T":4,"N":5,"-":0}
number_dic=dict([ (value, key) for key,value in base_dic.items() ])

class GenotypeSNP:
    UNINFORMATIVE_ALLLE="uninformative"

    """Contains information on genotypes implied by alleles

    :param contig_id: Name of contig that contans allele
    :type contig_id: str
    :param position: Position of allele on contig
    :type position: int
    :param reference_nucl: Reference nucleotide value
    :type reference_nucl: str    
    """
    def __init__(self, contig_id: str, position: int, reference_nucl: str) -> None:
        self._genotypes={} #key = allele: str, value = genotype :str
        self._contig_id = contig_id
        self._ref_nucl = reference_nucl
        self._position = int(position)

    @property
    def contig_id(self) -> str:
        return self._contig_id

    @property
    def position(self) -> int:
        return self._position

    @property
    def genotypes(self):
        return self._genotypes

    @property
    def reference_nucl(self) -> str:
        return self._ref_nucl

    def get_genotype(self, allele: str) -> str:
        """Return genotype for give allele value

        :param allele: Nucleotide sequence at target site
        :type allele: str

        :return Either genotype or, "uninformative" if allele
        is not informative
        :rtype: str
        """
        if allele in self._genotypes:
            return self._genotypes[allele]
        else:
            return GenotypeSNP.UNINFORMATIVE_ALLLE

    def add_genotype(self, allele: str, genotype: str=None) -> None:
        """Contains information on genotypes implied by alleles
        
        :param allele: nucleotide at target position
        :type allele: str
        :param genotype: genotype implied by the allele
        :type genotype: str
        """
        self._genotypes[allele]=genotype if not genotype is None else GenotypeSNP.UNINFORMATIVE_ALLLE

    def add_from_vcf(self, vcf_line: str) -> None:
        '''Add genotype from VCF line

        param vcf_line: SNP line from VCF file
        type vcf_line: str
        '''
        vcf_values=vcf_line.strip().split("\t")
        vcf_values[1]=int(vcf_values[1])-1 #correct for beds and BAMs being 0 indexed and VCF 1 indexed
        if not self.coordinates_match(vcf_values[0], vcf_values[1]):
            raise ValueError("Incorrect VCF line was given to GenotypeSNP object, either contig or position mismatch")
        if vcf_values[2].find(":")==-1:
            self.add_genotype( vcf_values[4] )
        else:
            gt_values=vcf_values[2].split(":")
            self.add_genotype(vcf_values[4], gt_values[0])

    def coordinates_match(self, contig: str, position: int) -> bool:
        return contig==self.contig_id and int(position)==self.position

#class Report

class ReadsMatrix:
    def __init__(self, bam_file: str, is_positive: bool) -> None:
        bam_file=bam_file
        if not exists(bam_file):
            raise ValueError(f'Bam file {bam_file} does not exist, did you specify just file name instead of full address?')
        self._bam_file: str=bam_file
        self._is_positive: bool=is_positive
        self._read_matrices: Dict[str, np.array]=None

    @property
    def bam_file(self) -> str:
        return self._bam_file

    @property
    def is_positive(self) -> bool:
        return self._is_positive
    
    @property
    def read_matrices(self) -> np.array:
        return self._read_matrices
    
    def get_read_matrix(self, amplicon_name: str) -> np.array:
        if amplicon_name in self._read_matrices:
            return self._read_matrices[amplicon_name]
        else:
            ValueError(f'Amplicon {amplicon_name} not found in file: {self._bam_file}')

    def read_bam(self, target_regions: List[Amplicon], **kwargs) -> None:
        full_length_only = kwargs.get('full_length_only', True)
        self._read_matrices:  Dict[str, np.array]={}
        bam = ps.AlignmentFile(self.bam_file, "rb",check_sq=False)
        for target_amplicon in target_regions:#amplicon_coordinates.index:
            target_start, target_end = [target_amplicon.ref_seq.ref_start, target_amplicon.ref_seq.ref_end]
            target_contig=target_amplicon.ref_seq.refseq_id
            alignment_matrix=np.zeros( (bam.count(contig=target_contig, start=target_start, end=target_end), target_end-target_start), dtype=np.int8)
            used_rows=[False]* alignment_matrix.shape[0] #this will allow removal of empty rows when only full length reads are required
            orientation=[]
            for i,read in enumerate(bam.fetch(contig=target_contig, start=target_start, end=target_end)):
                if not read.is_unmapped:
                    if read.query_sequence==None:
                        warnings.warn(f'Bam file {self.bam_file} has an empty read aligning to {target_amplicon.name} this should not happen')
                        continue
                    if full_length_only and (read.reference_start>target_start or read.reference_end<target_end):
                        continue #This will skip non-full length queries
                    orientation.append(read.is_forward)
                    query_nt = [ read_ref_pair[0] for read_ref_pair in read.get_aligned_pairs() if not read_ref_pair[1] is None and
                                    read_ref_pair[1]<target_end and read_ref_pair[1]>=target_start and not read_ref_pair[0] is None]
                    ref_nt = [ read_ref_pair[1]-target_start for read_ref_pair in read.get_aligned_pairs() if not read_ref_pair[1] is None and
                                read_ref_pair[1]<target_end and read_ref_pair[1]>=target_start and not read_ref_pair[0] is None]
                    bases=set([ read.query_sequence[f] for f in query_nt])
                    if True in [f not in base_dic for f in bases]:
                        warnings.warn(f'Ambigous bases not supported: {self.bam_file} target {target_amplicon.name}. Target will be ignored')
                        break

                    alignment_matrix[i, ref_nt]=[ base_dic[read.query_sequence[f]] for f in query_nt]
                    used_rows[i]=True

            self._read_matrices[target_amplicon.name]=alignment_matrix[used_rows]

class ClassificationResult:
    """Class to store classfication predictions and
    their asssesed quality

    :param name: Name for the model, amplicon name is the most useful
    :type name: str
    """
    def __init__(self, amplicon: Amplicon) -> None:
        self._amplicon: Amplicon=amplicon
        self._consensus=""
        self._mismatches={}
        self._predicted_classes=None
        self._mismatches=[]
        self._genotype_snps=[]

    def get_consensus(self, data: np.array) -> str:
        most_common=np.apply_along_axis( lambda x: Counter(x).most_common(1)[0][0], axis=0,  arr=data  )
        self.consensus="".join([number_dic[f] for f in most_common  ])
        return self.consensus

    def calculate_genotype(self, genot_snps:List[GenotypeSNP]) -> str:
        result=""
        for pos, alt in self._mismatches.items():
            known_alleles=[f for f in genot_snps if f.position==pos and f.contig_id==self._amplicon.ref_contig]
            if len(known_alleles)>0 :
                if alt in known_alleles[0].genotypes:
                    result=result + known_alleles[0].genotypes[alt]
                else:
                    result=result + f'Unknown GT at known position: {pos}'
        return result
               
    #region
    @property
    def predicted_classes(self) -> np.array:
        return self._predicted_classes

    @predicted_classes.setter
    def predicted_classes(self, value: np.array):
        self._predicted_classes = value

    @property
    def positive_cases(self) -> int:
        return Counter(self._predicted_classes)[0]

    @property
    def negative_cases(self) -> int:
        return Counter(self._predicted_classes)[1]

    @property
    def amplicon(self) -> Amplicon:
        return self._amplicon

    @amplicon.setter
    def amplicon(self, value: Amplicon):
        self._amplicon = value

    @property
    def consensus(self) -> str:
        return self._consensus

    @consensus.setter
    def consensus(self, value: str):
        self._mismatches=dict([ (pos,alt) for ref, alt, pos in zip(self.amplicon.seq, value, range(0,len(value))) if ref!=alt ])
        self._consensus = value

    @property
    def has_mismatches(self) -> bool:
        if self._consensus=="":
            raise ValueError("No consensus sequence, check if mismatches were calculated")
        elif self.num_mismatches==0:
            return False
        else:
            return True
    @property
    def num_mismatches(self) -> int:
        return len(self._mismatches)

    def result_description(self, genotype_snps: List[GenotypeSNP]=None) -> str:
        prefix=f'Total {len(self.predicted_classes)} mapping reads of which '+"{:.0%}".format(Counter(self.predicted_classes)[1]/len(self.predicted_classes))+' are from target organism. \n'
        if self.num_mismatches==0:
            return f'{self.amplicon.name}: {prefix} Consensus is identical to reference'
        else:
            if genotype_snps is None:
                genotypes="No genotype SNPs provided"
            else:
                genotypes="SNPs from genotypes: "+self.calculate_genotype(genot_snps=genotype_snps)
            snps=",".join([f'{str(pos+1)}:{self.amplicon.seq[pos]}->{alt}' for pos, alt in self._mismatches.items()])
            return f'{self.amplicon.name}: {prefix} Total SNPs: {self.num_mismatches}, {genotypes}, Mismatched alleles: {snps}'
    
    #endRegion

class Classifier:
    """Class that both trains model and applies model to classify new data

    :param name: Name for the model, amplicon name is the most useful
    :type name: str
    :param nucleotide_seq: The reference nucleotide sequence to which reads are mapped
    :type nucleotide_seq: str    
    :param genotype_snps: Information on which positions contain genotype SNPs
    :type genotype_snps: List[GenotypeSNP]
    """
    def __init__(self, name: str, nucleotide_seq:str, genotype_snps: List[GenotypeSNP]) -> None:
        self.name=name
        self._sensitivity=-1.0
        self._specificity=-1.0
        self._training_samples=-1
        self._testing_samples=-1
        self._testing_samples_misclassfied=-1
        self._nucletoide_seq=nucleotide_seq
        self._genotype_snps: List[GenotypeSNP]=genotype_snps
        self._not_trained=False
        self._model=None
        self._variable_columns:List[int]=[]

    def _nucleotide_to_binary_matrix(self, matrix: np.array) -> np.array:
        """Convert a matrix where nucleotides are represented by A:1,C:2,G:3,T:4,N:5,-:0
        to a binary matrix of dummy varialbes

        :param matrix: non-binary matrix to convert
        :type matrix: np.array

        :return binary matrix of nucleotides
        :rtype: np.array
        """
        #return  matrix
        ouput=np.zeros( (matrix.shape[0], matrix.shape[1]*len(base_dic)), dtype=int )
        for i, key in enumerate(number_dic.keys()):
            ouput[:, (0+i)*matrix.shape[1]: (1+i)*matrix.shape[1] ]=matrix==key
        return ouput

    def _generate_train_test_sets(self, matrix, n_train: int, n_test: int):
        """Splits a reads matrix into train and test sets of required sizes

        :param matrix: matrix to split
        :type matrix: np.array
        :param n_train: Required number of reads in training set
        :type n_train: int
        :param n_test: Required number of reads in test set
        :type n_test: int

        :return Matrix of reads for training, matrix of reads for testing
        :rtype: (np.array, np.array)
        """        
        binarised_matrix=self._nucleotide_to_binary_matrix(matrix)

        sample_train, sample_test = train_test_split( range(0, binarised_matrix.shape[0] ),
                                                train_size=n_train, test_size=n_test    )

        reads_train = binarised_matrix[sample_train,:]
        reads_test = binarised_matrix[sample_test,:]

        # sample_train, sample_test = train_test_split( range(0, matrix.shape[0] ),
        #                                         train_size=n_train, test_size=n_test    )

        return reads_train, reads_test, sample_test

    def generate_test_train_sets(self, positive_matrix, negative_matrix, positive_n_train, negative_n_train, positive_n_test, negative_n_test, negative_bams):
        """Splits positive and negative matrices into train and test sets

        :param positive_matrix: matrix of reads from positive cases
        :type positive_matrix: np.array
        
        :param negative_matrix: matrix of reads from negative cases
        :type positive_matrix: np.array

        :param positive_n_train: number of positive reads to use for training
        :type positive_n_train: np.array

        :param negative_n_train: number of negative reads to use for training
        :type negative_n_train: np.array

        :param positive_n_test: number of positive reads to use for testing
        :type positive_n_test: np.array

        :param negative_n_test: number of negative reads to use for testing
        :type negative_n_test: np.array

        :return Four matrices: reads for training, true class for training, 
        reads for testing, true class for testing
        :rtype: (NDArray, NDArray, NDArray, NDArray)
        """

        positive_train, positive_test, positive_test_samples = self._generate_train_test_sets(positive_matrix,positive_n_train, positive_n_test)
        negative_train, negative_test, negative_test_samples = self._generate_train_test_sets(negative_matrix,negative_n_train, negative_n_test)
        #test is the quality of differentiation between reads based on the typhi specific alleles
        positive_matrix=self._nucleotide_to_binary_matrix(positive_matrix)
        negative_matrix=self._nucleotide_to_binary_matrix(negative_matrix)

        reads_train = np.vstack( ( positive_train , negative_train ) )
        reads_test = np.vstack( ( positive_test , negative_test ) )
        true_class_train=[1]*positive_train.shape[0] +[0]*negative_train.shape[0]
        true_class_test=[1]*positive_test.shape[0] +[0]*negative_test.shape[0]

        self.test_samples_ids=["Positive"] * positive_n_test
        self.test_samples_ids=self.test_samples_ids+[f[0] for f in negative_bams[negative_test_samples] ] #for debugging
       
        return reads_train, true_class_train,  reads_test, true_class_test

    def train_classifier(self, positive_data, negative_data, variable_columns, negative_data_bams) -> None:
        """Splits a reads matrix into train and test sets of required sizes

        :param positive_data: matrix of reads from target organism/genotype
        :type positive_data: np.array
        :param negative_data: matrix of reads from non-target organism/genotype
        :type negative_data: np.array
        :param variable_columns: list of columns which will be excluded from model training.
            Captures variable positions such as AMR mutations
        :type variable_columns: List[int]
        :param negative_data_bams: names of BAM files, these will help identify cause of misclassifications
        :type negative_data_bams: np.array
        """
        self._variable_columns=variable_columns
        used_columns=np.setdiff1d(range(0, positive_data.shape[1]), self._variable_columns)

        if positive_data.shape[0]<100 or negative_data.shape[0]<100: #the train test set is too small. Model won't be trained and will rely on sequence similarity instead
            self._train_without_negative(positive_data)
            self.not_trained=True
            self.sensitivity=1
            self.specificity=1
            self.testing_samples_misclassfied=0
            self.testing_samples=0
            return
        else:
            self.not_trained=False
            train_data, train_true,  test_data, test_true = self.generate_test_train_sets(positive_data[:,used_columns], negative_data[:,used_columns],
                                                                                                            min(int(positive_data.shape[0]*0.8), 5000),
                                                                                                            min(int(negative_data.shape[0]*0.8), 5000),
                                                                                                            min(int(positive_data.shape[0]*0.2), 5000),
                                                                                                            min(int(negative_data.shape[0]*0.2), 5000), negative_data_bams )
            self._train_with_negative(train_data, train_true)

            
            failed_predictions=[k!=z for k, z in zip(test_true, self._classify(test_data)) ]

            prediction_matrix=confusion_matrix(test_true, self._classify(test_data) )
            self.sensitivity=prediction_matrix[1,1]/ np.count_nonzero(test_true)
            self.specificity=prediction_matrix[0,0]/ (len(test_true)- np.count_nonzero(test_true))
            print(prediction_matrix)
            print(f'sens: {str(self.sensitivity)}, spec:{str(self.specificity)}')
            self.misclassified_inputs=Counter([f.replace(".bam","") for f in np.asarray(self.test_samples_ids)[failed_predictions] if f!="Positive"])
            print(self.misclassified_inputs)
            self.testing_samples_misclassfied=len(Counter([f.replace(".bam","") for f in np.asarray(self.test_samples_ids)[failed_predictions] if f!="Positive"]))
            self.testing_samplestotal_negative_samples=len(Counter(self.test_samples_ids))
            del train_data, train_true,  test_data, test_true


    def _train_without_negative(self, data) -> None:
        '''In some cases, there are no whole lenght reads from
        non-target organism that map to amplicon. As a result, classifier
        model can't be trained. This alternative approach 
        Calculates the distance of each read from reference
        and determines the distribution of distances. This forms
        basis of calling reads source from unknown data.
        
        :param data: numpy array of reads from target organism
        :type data: NDArray
        '''
        distances=self._get_distances(data).reshape(-1,1)
        self._model = GaussianMixture(n_components=1, random_state=0).fit(distances)

    def _train_with_negative(self, train_data, train_classes) -> None:
        '''Trains model to classify data using an ensemble model

        :param train_data: Matrix of reads (positive and negative)
        :type train_data: NDArray

        :param train_classes: Vector of true classes
        :type train_classes: NDArray
        '''
        ### Classification
        #clf2 = DecisionTreeClassifier(criterion="entropy")
        #clf4 = RandomForestClassifier(criterion="log_loss", n_estimators=10, random_state=1)
        #clf2 = SVC(kernel="linear", probability=True, random_state=0)

        #clf2 = GaussianNB()
        #clf2=KNeighborsClassifier(n_neighbors=5)
        #clf2=MLPClassifier(learning_rate="adaptive")
        #clf3=CategoricalNB(min_categories=2)
        clf2=SGDClassifier(max_iter=1000, tol=1e-5, loss="log_loss")

        eclf = VotingClassifier(
                    estimators=[('rc', clf2)],
                    voting='hard')
        
        # eclf = VotingClassifier(
        #             estimators=[('rc', clf2),  ('gnb', clf3)],
        #             voting='soft')

        #params = {'lr__C': [1.0, 100.0], 'rf__n_estimators': [20, 200]}
        #grid = GridSearchCV(estimator=eclf, param_grid=params, cv=5)

        self._model = eclf.fit(train_data, train_classes)

    def _get_distances(self, reads: np.array) -> np.array:
        ref_vector=[base_dic[f] for f in self._nucletoide_seq[1:]]
        distances=np.apply_along_axis(lambda x: len(ref_vector)-sum(x==ref_vector), axis=1, arr=reads[:,:-1])
        return distances

    def _classify(self, data: np.array):
        if data.shape[0]==0:
            return np.asarray([])

        if self.not_trained:
            #prediction based on reads being within 50% CI of
            #positive cases
            distances=self._get_distances(data)
            upper_limit=self._model.means_[0][0]+0.68*self._model.covariances_[0][0][0] / self._model.n_features_in_ **(1/2)
            lower_limit=self._model.means_[0][0]-0.68*self._model.covariances_[0][0][0] / self._model.n_features_in_ **(1/2)
            predictions=np.asarray([1 if f<=upper_limit and f>=lower_limit else 0 for f in distances ])
        else:
            predictions=self._model.predict(data)
        return predictions

    def classify_new(self, data: np.array) -> np.array:
        """Apply model's classifier to categories data

        :param data: Matrix of reads
        :type data: NDArray
        """
        if not self.not_trained:
            used_columns=np.setdiff1d(range(0, data.shape[1]), self._variable_columns)
            data=self._nucleotide_to_binary_matrix(data[:,used_columns])
        predictions=self._classify(data)
        return predictions

    #region Properties
    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, value):
        self._model = value

    @property
    def genotype_snps(self) -> List[GenotypeSNP]:
        return self._genotype_snps

    @genotype_snps.setter
    def genotype_snps(self, value: List[GenotypeSNP]):
        self._genotype_snps = value

    @property
    def sensitivity(self) -> float:
        return self._sensitivity

    @sensitivity.setter
    def sensitivity(self, value: float):
        self._sensitivity = value

    @property
    def specificity(self) -> float:
        return self._specificity

    @specificity.setter
    def specificity(self, value: float):
        self._specificity = value

    @property
    def not_trained(self) -> bool:
        return self._not_trained

    @not_trained.setter
    def not_trained(self, value: bool):
        self._not_trained = value

    @property
    def training_samples(self) -> int:
        return self._training_samples

    @training_samples.setter
    def training_samples(self, value: int):
        self._training_samples = value

    @property
    def testing_samples(self) -> int:
        return self._testing_samples

    @testing_samples.setter
    def testing_samples(self, value: int):
        self._testing_samples = value

    @property
    def testing_samples_misclassfied(self) -> int:
        return self._testing_samples_misclassfied

    @testing_samples_misclassfied.setter
    def testing_samples_misclassfied(self, value: int):
        self._testing_samples_misclassfied = value

    @property
    def nucletoide_seq(self) -> float:
        return self._nucletoide_seq

    @nucletoide_seq.setter
    def nucletoide_seq(self, value: float):
        self._nucletoide_seq = value
    #endregion
