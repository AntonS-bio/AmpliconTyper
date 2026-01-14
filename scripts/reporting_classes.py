from pathlib import Path
from os.path import basename
from datetime import datetime
from typing import List, Dict, Tuple
from collections import Counter
from read_classifier import ClassificationResult, Classifier, GenotypeSNP, number_dic, ModelsData
import pandas as pd

OPENER='<tr><td>'
OPENER_BOLD='<tr class=table_header><td>'
OPENER_ALT='<tr class=alt_row><td>'
CLOSER='</td></tr>'
MIDDLE='</td><td>'
MIDDLE_WIDE='</td class=td_alt_wide><td>'
GT_REPORTING_MIN_READ_PERCENT=0.05
POSITIVE_CASES_THRESHOLD=50
CORRECT_TO_WRONG_LEN_RATIO=3.0
POSITIVE_CASES_MIN_AMPLICONS=2
MAX_ALLOWED_SNPS=1
EXPECTED_ERROR_RATE=0.03

class ReportingUtilities:
    def __init__(self):
        pass

    def get_styles(self):
        styles=""
        with open( Path(__file__).with_name('html_head.txt') ,"r" ) as styles_data:
            for line in styles_data:
                styles+=line
            styles+="\n"
        return styles
    
    def clean_sample(self,sample_name: str) -> str:
        return Path(sample_name).stem

    def insert_paragraph(self, number: int):
        text=""
        for i in range(0,number):
            text+="<br/>"
        text+="\n"
        return text

    def get_most_recent_model_element(self, results: List[ClassificationResult]):
        models_list=sorted(set([(f.model_timestamp) for f in results]))
        return models_list[-1].strftime("%d/%b/%y")

    def generate_header_row(self, header_values: List[str], tooltip_values: List[str]="") -> str:
        if len(tooltip_values)==0:
            table_header=OPENER_BOLD+MIDDLE.join(header_values)+CLOSER+"\n"
        else:
            #tooltops are present, so add them to header cells
            values=[value+"  "+f'<img class="help-icon" title="{hint}">' for value, hint in zip(header_values, tooltip_values)]
            table_header=OPENER_BOLD+MIDDLE.join(values)+CLOSER+"\n"
        return table_header

class AlleleInfo:
    def __init__(self, ref: str, nucl: str, position: int, depth:int, frequency: float, amplicon_name = "", implication=""):
        self._nucl=nucl
        self._ref=ref
        self._position=position
        self._frequency=frequency
        self._implication=implication
        self._amplicon_name=amplicon_name
        self._depth=depth
    
    @property
    def implication(self) -> str:
        return self._implication
    @implication.setter
    def implication(self, value: str):
        self._implication=value

    @property
    def depth(self) -> int:
        return self._depth
    @property
    def is_ref(self) -> bool:
        return self._nucl==self._ref
    @property
    def amplicon_name(self) -> str:
        return self._amplicon_name
    @property
    def nucl(self) -> str:
        return self._nucl
    @property
    def ref(self) -> str:
        return self._ref
    @property
    def pos(self) -> int:
        return self._position
    @property
    def freq(self) -> float:
        return self._frequency

class AmpliconReportValues:
    def __init__(self, result: ClassificationResult, model: Classifier):
        self._name=result.amplicon.name
        self._result: ClassificationResult=result
        self._is_transient=False
        self._all_gts: Tuple[str, bool]=("", False)
        self._num_unknown_snps=-1
        self._model: Classifier = model
        self._dominant_snps: Dict[int, str] = {}
        self._alleles: List[AlleleInfo] = self.calculate_all_alleles()

    @property
    def is_transient(self) -> bool:
        return self._model.amplicon_transient
    
    @property
    def amplicon_gts_str_format(self) -> str:
        result=""
        if self.is_amr:
            return result
        else:
            for snp in self.all_known_snps:
                result+=snp.implication+" ("+'{:.0%}'.format(snp.freq)+")"
                if int(self.target_org_reads*snp.freq)<POSITIVE_CASES_THRESHOLD:
                    result+="(?)"
                result+=", "
            result=result[0:-2]
        return result

    @property
    def amplicon_gts_list(self) -> List[AlleleInfo]:
        result=[]
        if self.is_amr:
            return result
        else:
            for snp in self.all_known_snps:
                result.append(snp)
        return result

    @property
    def amplicon_amrs(self) -> str:
        result=""
        if self.is_transient and self.target_org_reads>=POSITIVE_CASES_THRESHOLD:
            #capture detection of transient AMR amplicons
            return self.name

        if not self.is_amr:
            return result
        elif self.target_org_reads>=POSITIVE_CASES_THRESHOLD:
            for snp in self.all_known_snps:
                result+=snp.implication+" ("+'{:.0%}'.format(snp.freq)+")"
                # if self.target_org_reads<POSITIVE_CASES_THRESHOLD:
                #     result+="(?)"
        return result

    @property
    def model(self) -> Classifier:
        return self._model
    
    @model.setter
    def model(self, value: Classifier):
        self._model=value

    @property
    def not_trained(self) -> bool:
        return self.model.not_trained
    
    @property
    def consensus(self) -> str:
        return self._result.consensus
    
    @property
    def is_amr(self) -> bool:
        return True in [f.is_amr for f in self.all_defined_snps if f.contig_id==self.name]

    def calculate_all_alleles(self) -> List[AlleleInfo]:
        """Returns AlleleInfo for each position of the amplicon sequence
          Only removes the position where read depth is below POSITIVE_CASES_THRESHOLD
          or where allele % of total reads is below LOW_READ_WARNING_PERCENT
        """
        result=[]
        for pos, position_frequencies in enumerate(self._result._allele_frequencies):
            #if sum([f[1] for f in position_frequencies])> POSITIVE_CASES_THRESHOLD:
            for nt_code, read_count in position_frequencies:
                nt_freq=read_count/self.target_org_reads
                matching_known_snps=[f for f in self.all_defined_snps if f.position==pos and number_dic[nt_code] in f.genotypes]

                #check for presence of attribute is required for compatability with older models.
                if len(matching_known_snps)==1 and \
                        hasattr(matching_known_snps[0], '_reporting_threshold') and \
                        matching_known_snps[0].reporting_threshold:
                    freq_reporting_threshold=matching_known_snps[0].reporting_threshold
                else:
                    freq_reporting_threshold=GT_REPORTING_MIN_READ_PERCENT
                if nt_freq>freq_reporting_threshold:
                    result.append( AlleleInfo(self._result.amplicon.ref_seq.sequence[pos], number_dic[nt_code],
                                                pos, read_count, nt_freq, self.name ) )
        return result

    @property
    def alleles(self) -> List[AlleleInfo]:
        """Returns all sites (AlleleInfo class) 
        with variable nucleotides
        """
        return self._alleles
    
    @property
    def dominant_alleles(self) -> List[AlleleInfo]:
        result: Dict[int, AlleleInfo]={}
        for snp in self.alleles: #collect most frequent alleles across variable sites
            if snp.pos not in result:
                result[snp.pos]=snp
            elif result[snp.pos].freq<snp.freq:
                result[snp.pos]=snp
        return [f for pos, f in result.items() if not f.is_ref]

    @property
    def dominant_known_snps(self) -> List[AlleleInfo]:
        result: List[AlleleInfo] = []
        for allele in self.dominant_alleles:
            known_alleles=[f for f in self.all_defined_snps if f.position==allele.pos and allele.nucl in f.genotypes]
            if len(known_alleles)==1:
                allele.implication=known_alleles[0].genotypes[allele.nucl]
                result.append( allele )
            elif len(known_alleles)>1:
                raise ValueError(f'Allele {allele.nucl} at {self._result.amplicon.ref_contig} pos {str(allele.pos)} implies more than one genotype, check VCF supplied at training!')
        return result

    @property
    def all_known_snps(self) -> List[ AlleleInfo]:
        """Returns all SNPs detected in the data that match
          SNP defined in the amplicon classifier model.
        """
        if self.name=="1_3_4.3.1":
            a=1
        result: List[GenotypeSNP] = []
        for allele in self.alleles:
            known_alleles=[f for f in self.all_defined_snps if f.position==allele.pos and allele.nucl in f.genotypes]
            if len(known_alleles)==1:
                allele.implication=known_alleles[0].genotypes[allele.nucl]
                result.append( allele )
            elif len(known_alleles)>1:
                raise ValueError(f'Allele {allele.nucl} at {self._result.amplicon.ref_contig} pos {str(allele.pos)} implies more than one genotype, check VCF supplied at training!')
        return result

    @property
    def dominant_unknown_snps(self) -> List[AlleleInfo]:
        result: List[AlleleInfo] = []
        for allele in self.dominant_alleles:
            if not allele.is_ref:
                known_alleles=[f for f in self.all_defined_snps if f.position==allele.pos and allele.nucl in f.genotypes]
                if len(known_alleles)==0:
                    allele.implication="Unknown"
                    result.append( allele )
        return result

    @property
    def all_unknown_snps(self):
        """Returns all SNPs detected in the data that DO NOT match
          SNP defined in the amplicon classifier model.
        """
        result: List[AlleleInfo] = []
        for allele in self.alleles:
            if not allele.is_ref:
                known_alleles=[f for f in self.all_defined_snps if f.position==allele.pos and allele.nucl in f.genotypes]
                if len(known_alleles)==0:
                    allele.implication="Unknown"
                    result.append( allele )
        return result

    @property
    def all_defined_snps(self) -> List[GenotypeSNP]:
        """Returns all SNPs that are listed in the amplicon's classifier model,
        not just those found in the sample data
        """
        return self.model.genotype_snps

    @property
    def consensus_snp_implications(self) -> str:
        '''Summarises each SNP into its position, 
        allele and implication of the change based on the 
        SNPs defined in amplicon classifier model
        '''
        snp_values: List[str] = []
        consensus_snps=sorted(self.dominant_known_snps+self.dominant_unknown_snps, key=lambda x: x.pos)
        for snp in consensus_snps:
            allele_frequency=" ("+'{:.0%}'.format(self._result.consensus_frequency[snp.pos])+") "
            snp_value="Pos:"+str(snp.pos) + allele_frequency
            snp_value+=", Change:"+self._result.amplicon.ref_seq.sequence[snp.pos] + ">"+snp.nucl
            snp_value+=", Implication: "+snp.implication
            snp_values.append(snp_value)
        return '<br/>'.join(snp_values)

    @property
    def mapped_count(self) -> int:
        return self.wrong_length_count+self.right_length_count

    @property
    def target_org_share(self) -> int:
        if self.right_length_count>0:
            return '{:.2%}'.format(self.target_org_reads/self.right_length_count)
        else:
            return ""

    @property
    def target_org_reads(self) -> int:
        return self._result.positive_cases

    @property
    def wrong_length_count(self) -> int:
        return self._result.wrong_len_reads[self._name]
    
    @property
    def right_length_count(self) -> int:
        return self._result.positive_cases+self._result.negative_cases

    @property
    def name(self) -> str:
        return self._name

class ReportingGenotype:
    MISSING_VALUE="MISSING_VALUE"
    DEFAULT_HIGH_LEVEL_GT="Any"
    def __init__(self, hierarchy: Dict[str, Dict]=None):
        self._reporting_gt=""
        self._hierarchy=hierarchy
        self._levels: Dict[str, int]={} #key = genotype, value = depth of nesting
        if not hierarchy is None:
            self.generate_levels()
        
    def generate_levels(self):
        #the hierarchy data is nested dictionary {"1": {"type": "genotype", children: { "2": "type": "genotype", "children"     }    } }
        self._levels.clear()
        for key, values in self._hierarchy.items():
            self._levels[key]=1
            if "children" in values:
                self.generate_levels_helper(values["children"], 2)

    def generate_levels_helper(self, data: Dict[str, Dict], current_level) -> int:
        for key, values in data.items():
            self._levels[key]=current_level
            if "children" in values:
                self.generate_levels_helper(values["children"], current_level+1)

    def get_lowest_level_gt(self, alleles: List[AlleleInfo]) -> str:
        if len(alleles)==0:
            return ""
        else:
            lowest_gt=None
            for allele in alleles:
                if allele.implication in self._levels:
                    if lowest_gt is None or self._levels[allele.implication]>self._levels[lowest_gt.implication]:
                        lowest_gt=allele
            if lowest_gt is None: #the identified allele(s) are not genotype specific
                return ""
            else:
                return lowest_gt.implication+" ("+'{:.0%}'.format(lowest_gt.freq)+")"

    def possible_high_level_gts(self, alleles: List[AlleleInfo], genotypes_matrix: pd.DataFrame) -> List[str]:
        '''Returns possible values for high level genotypes
        Only makes sense for the organisms like S. Typhi where
        some genotypes (0, 1, 2, 3, 4 in Typhy) cannot be defined by single SNP
        Each row of genotypes_matrix corresponds to single high_level genotype allele
        index is Tuple of amplicon name and position. 
        ['1_3_4.3.1',381, 'C', 'C', 'C', 'T','T']
        columns names are ["Amplicon", "Pos","0","1","2","3","4"] where 0 - 4 are genotypes
        '''
        present_gts=pd.DataFrame(1, index=genotypes_matrix.index, columns=genotypes_matrix.columns, dtype=int)
        for allele in alleles:
            if allele.depth>POSITIVE_CASES_THRESHOLD:
                if (allele.amplicon_name, allele.pos) in genotypes_matrix.index:
                    possible_gts = genotypes_matrix.loc[ (allele.amplicon_name, allele.pos)] != allele.nucl
                    present_gts.loc[ (allele.amplicon_name, allele.pos), possible_gts] = 0
        possible_gts=list(present_gts.columns[present_gts.product()==1])
        if len(possible_gts)==len(present_gts.columns):
            return ReportingGenotype.DEFAULT_HIGH_LEVEL_GT
        else:
            return ", ".join(possible_gts)

class SampleReportValues:
    
    def __init__(self, sample_name:str, amplicon_results: List[ClassificationResult], model_data: ModelsData):
        self._amplicon_results: Dict[str, AmpliconReportValues] = {}
        self._model_data=model_data
        self._reporting_genotype=ReportingGenotype( model_data.metadata.get("hierarchy",{}) )
        self._utilities=ReportingUtilities()
        for result in amplicon_results:
            new_report_amplicon=AmpliconReportValues(result, model=model_data.classifiers[result.amplicon.name])
            self._amplicon_results[result.amplicon.name]=new_report_amplicon

        all_alleles=[allele for amplicon, result in self._amplicon_results.items() for allele in result.alleles]
        if "high_level_genotypes" in model_data.metadata:
            self._possible_highlevel_gts=self._reporting_genotype.possible_high_level_gts( alleles=all_alleles, genotypes_matrix=model_data.metadata["high_level_genotypes"] )
        else:
            self._possible_highlevel_gts=None

        self._name=sample_name
        self._sample_descrition="No description"
        self._sequencing_read_count = ""

    def snp_frequency_diagram(self):
        HTML_SPACE = "&nbsp;"
        unknown_snps=Counter( [len(f.dominant_unknown_snps) for k, f in self._amplicon_results.items()
                               if not f.is_transient and f.target_org_reads>=POSITIVE_CASES_THRESHOLD] )
        diagram=""
        for x in range(0, 5):
            if x in unknown_snps:
                diagram=diagram+str(unknown_snps[x])+","+HTML_SPACE
            else:
                diagram=diagram+"0,"+HTML_SPACE
        over_five=sum([value for key, value in unknown_snps.items() if key>=5])
        diagram=diagram+str(over_five)
        return diagram

    @property
    def model_data(self) -> ModelsData:
        '''Returns and object containing Classifier models and metadata for the models'''
        return self._model_data

    @property
    def target_org_support(self) -> str:
        '''Returns and object containing Classifier models and metadata for the models'''
        valid_amplicons=sum([1 if not f.is_amr and f.target_org_reads>50 and len(f.all_unknown_snps)<=2 else 0 for f in  self._amplicon_results.values()])
        if valid_amplicons>=4:
            return "Strong"
        elif valid_amplicons>=2:
            return "Weak"
        else:
            return "None"


    @property
    def display_order(self):
        return sorted(self.amplicon_results.items(),
                        key=lambda x: (x[1].is_transient, x[1].is_amr, x[0]))
    @property
    def amplicon_results(self) -> Dict[str, AmpliconReportValues]:
        return self._amplicon_results

    @property
    def target_positive(self) -> bool:
        '''Indicates if sample has target organism'''
        amplicon_with_sufficient_reads=0
        for name, result in  self.amplicon_results.items():
            if not result.is_amr and not result.is_transient:
                if result.target_org_reads>POSITIVE_CASES_THRESHOLD and len(result.dominant_unknown_snps)<=MAX_ALLOWED_SNPS:
                    amplicon_with_sufficient_reads+=1
        return amplicon_with_sufficient_reads>=POSITIVE_CASES_MIN_AMPLICONS

    @property
    def possible_highlevel_gts(self) -> str:
        '''Returns the string showing possible high level genotypes
        for the bacteria where some genotypes are defined by multiple SNP
        '''
        if not self.target_positive:
            return "-"
        return self._possible_highlevel_gts

    @property
    def agg_gts(self) -> str:
        if not self.target_positive:
            return "-"
        gts: List[AlleleInfo]=[]
        for amplicon_name, value in self._amplicon_results.items():
            gts+=value.amplicon_gts_list
        return self._reporting_genotype.get_lowest_level_gt(gts)

    @property
    def agg_amrs(self) -> str:
        if not self.target_positive:
            return "-"
        amrs=""
        for amplicon_name, value in self._amplicon_results.items():
            if amrs=="":
                amrs=value.amplicon_amrs
            elif value.amplicon_amrs!="":
                amrs+=", "+value.amplicon_amrs
        return amrs

    @property
    def display_details(self) -> bool:
        #Meant to accomodate hiding of the results if very little data is present
        return True

    @property
    def sample_description(self) -> str:
        return self._sample_descrition

    @sample_description.setter
    def sample_description(self, value: str):
        self._sample_descrition=value

    @property
    def name(self) -> str:
        sample_cell_value=self._utilities.clean_sample(self._name)
        return sample_cell_value
    
    @property
    def name_with_description(self) -> str:
        sample_cell_value=self._utilities.clean_sample(self._name)
        return sample_cell_value+"<br/>"+self.sample_description

    @property
    def sequencing_read_count(self) -> str:
        return self._sequencing_read_count
    @sequencing_read_count.setter
    def sequencing_read_count(self, value: int):
        self._sequencing_read_count=str(value)

    @property
    def full_name(self) -> str:
        return self._name
    @property
    def mapped_count(self) -> int:
        return sum([values.mapped_count for name, values in self._amplicon_results.items()])
    @property
    def wrong_length_count(self) -> int:
        return sum([values.wrong_length_count for name, values in self._amplicon_results.items()])
    @property
    def right_length_count(self) -> int:
        return sum([values.right_length_count for name, values in self._amplicon_results.items()])
    @property
    def amplicon_snp_profile(self) -> str:
        return ""

class SummaryTable:
    utilities=ReportingUtilities()
    def __init__(self, model_date: str, model_file: str, include_optional_columns="False" ):
        self._model_file=  basename(model_file).replace(".pkl","")
        self._model_date=model_date
        self._include_optional_columns = include_optional_columns
        self._show_len_ratio_warning=False
        self._sequencing_error_rate=-1
        
    OPTIONAL_COLUMN=["Possible High Level Genotypes"]
    FIRST_LINE='<div class=header_line><a name="Summary">Summary of results</div>\n'
    HEADER_VALUES=["Sample", "Support for Target Organism",  "Possible High Level Genotypes", "Implied Genotypes", "Implied AMR",
                   "Sequencing Reads","Mapped, correct length", "Mapped, but too short/long", "Amplicons by # of dominant unknown SNPs"]#," Highest read count","Lowest read count"]
    HEADER_TOOLTIP_TEXT=["Sample name, you can make it more informative by adding sample description using options -c and -d.", #Sample
                         "Based on number of valid amplicons: Strong is 4+, Weak is 2-3 and None 0-1. Valid amplicon has to be&#10\
                         a) not AMR related&#10b)have at least 50 target reads&#10c)have <=2 uknown SNPs",
                         "Relevant for some organisms (e.g. S. Typhi) where some genotypes cannot be determined with single SNP.&#10\
Shows high level genotypes (in S. Typhi these are 0, 1, 2, 3, 4) that are possible given detected amplicons.&#10\
Uses all alleles that are supported by >"+str(POSITIVE_CASES_THRESHOLD)+"reads.&#10'Any' means detected amplicons provide no information on the high level genotypes.", #HL genotypes
                         "Shows only the most specific genotype based on detected amplicons. Will not show multiple genotypes even if sample consists of more than one.", #Genotype
                         "Blank field means AMR amplicons not detected,&#10Amplicon name means gene detected, but doesn't have AMR mutations,&#10\
Mutation description means this mutation was detected in respective amplicon", #Implied AMR
                          "Total reads in FASTQ files. If classification is done from BAM files this field will be blank.", #Total reads
                          "Reads which mapped to some amplicon and have length that meets -l and -s parameters.", #Corrent lengths
                          "Reads which mapped to some amplicons, but have length that does not meet -l and -s parameters.&#10See values used just below 'Summary of results' table",
                            "Numbers means number of amplicons, number's position-1 indicates how many unknown SNPs these amplicons have.&#10\
                            E.g. 2,0,4,0,0 means there are 2 amplicons without unknown SNPs (position is 1, so 1-1=0), and \
                            4 amplicons with 2 unknown SNPs each (position is 3-1=2).&#10\
                            Only examines the amplicons that should be always present i.e. excludes plasmids, recombinations, etc."] #Amplicon/SNP diagram                          ] #Too long

    def calculate_error_rate(self, results: List[ClassificationResult])-> float:
        '''Calculates the expected error rate based on PHRED/Q scores and
        weighted by read counts'''
        self._sequencing_error_rate=sum( [(f.negative_cases+f.positive_cases)*f.read_error_rate for f in results] )/sum( [(f.negative_cases+f.positive_cases) for f in results] )
        return self._sequencing_error_rate
    
    def get_error_rate_blurb(self) -> str:
        if self._sequencing_error_rate==-1:
            return "<div>Something went wrong - sequencing error rate not calculated.</div>"
        else:
            if self._sequencing_error_rate<EXPECTED_ERROR_RATE*1.5:
                return f'<div>Expected error rate based on Q-scores: {self._sequencing_error_rate:.2%}. This is within the expected range.</div>'
            else:
                result='<div style="color:red">'
                result+=f'Expected error rate based on Q-scores:{self._sequencing_error_rate:.2%}. This is much higher than expected range (~{EXPECTED_ERROR_RATE:.2%})'
                result+=" and will make genotype and AMR calls very unreliable. This is a sequencing, not analysis problem."
                result+=" Please check MinKNOW to identify possible causes.</div>"
                return  result
            
    def get_main_header(self) -> str:
        main_header=SummaryTable.FIRST_LINE
        if not self._include_optional_columns:
            for optional_value in self.OPTIONAL_COLUMN:
                self.HEADER_VALUES.remove(optional_value)
        main_header+='<div>Model name: '+ self._model_file + \
                        ', model date: '+self._model_date+'</div>\n'
        main_header+=self.get_error_rate_blurb()
        main_header+=self.utilities.insert_paragraph(1)
        main_header+="<table>\n"
        main_header+="\t<tbody>\n"
        main_header+=self.utilities.generate_header_row(SummaryTable.HEADER_VALUES, SummaryTable.HEADER_TOOLTIP_TEXT)
        return main_header
    
    def get_main_tail(self, length_parameters: Dict[str, int]) -> str:
        main_tail="\t</tbody>\n"
        main_tail+="</table>\n"
        #self.output_file.write(f'<p>* Defined as at least {self._positive_amplicons_threshold*100}% of amplicons that <a href="#ModelSignatures">are part of core genome</a> having at least {POSITIVE_CASES_THRESHOLD} target organism reads</p>')
        main_tail+=f'<p>"-" fewer than {POSITIVE_CASES_MIN_AMPLICONS} amplicons with {POSITIVE_CASES_THRESHOLD} target organism reads and fewer than {MAX_ALLOWED_SNPS+1} SNPs</p>'
        if self._show_len_ratio_warning:
            main_tail+=f'<p>^ The number of reads either longer or shorter than amplicon is >{str(CORRECT_TO_WRONG_LEN_RATIO)}x higher than reads that match amplicon length. If this is unexpected, consider using options "-s" or "-l". </p>'
        main_tail+=f'<p>Classifier allowed soft clip (-s) of {length_parameters["-s"]}nt and 10nt length difference (-l) between reference and mapped read.'
        main_tail+='<p>Report generated on '+datetime.now().strftime("%d/%b/%y at %H:%M")+'. <a href="#ModelSignatures">See models list</a></p>'
        return main_tail

    def assemble_row_values(self, sample: SampleReportValues, wrong_len_suffix: str=""):
        if self._include_optional_columns:
            values=[sample.target_org_support, sample.possible_highlevel_gts, sample.agg_gts, sample.agg_amrs, sample.sequencing_read_count,
                sample.mapped_count-sample.wrong_length_count, str(sample.wrong_length_count)+wrong_len_suffix,sample.snp_frequency_diagram() ]
           
        else:
            values=[sample.target_org_support, sample.agg_gts, sample.agg_amrs, sample.sequencing_read_count,
                sample.mapped_count-sample.wrong_length_count, str(sample.wrong_length_count)+wrong_len_suffix,sample.snp_frequency_diagram() ]
        return values

    def get_table_row(self, sample: SampleReportValues, alt_style=False) -> str:
        sample_name_href='<a href="#'+sample.name+'">'+sample.name_with_description+"</a>"
        if sample.mapped_count==sample.wrong_length_count:
            lengths_ratio = CORRECT_TO_WRONG_LEN_RATIO+1 if sample.mapped_count>100 else 0
        else:
            lengths_ratio=sample.wrong_length_count/(sample.mapped_count-sample.wrong_length_count)
        if not self._show_len_ratio_warning and lengths_ratio>CORRECT_TO_WRONG_LEN_RATIO:
            self._show_len_ratio_warning=True #will display information that 3.0x read are rejected
            #and suggest the adjustment to -s and -l parameters to fix this.
        wrong_len_suffix="^" if lengths_ratio>CORRECT_TO_WRONG_LEN_RATIO else ""
        values=[sample_name_href]+self.assemble_row_values(sample, wrong_len_suffix)#+[wrong_len_suffix]
        # if self._include_optional_columns:
        #     values=[sample_name_href, sample.snp_frequency_diagram(), sample.possible_highlevel_gts, sample.agg_gts, sample.agg_amrs, sample.sequencing_read_count,
        #         sample.mapped_count-sample.wrong_length_count, str(sample.wrong_length_count)+wrong_len_suffix ]
           
        # else:
        #     values=[sample_name_href, sample.snp_frequency_diagram(), sample.agg_gts, sample.agg_amrs, sample.sequencing_read_count,
        #         sample.mapped_count-sample.wrong_length_count, str(sample.wrong_length_count)+wrong_len_suffix ]

        if alt_style: #Alternate row colours
            table_line=OPENER_ALT+MIDDLE.join([str(f) for f in values])+CLOSER
        else:
            table_line=OPENER+MIDDLE.join([str(f) for f in values])+CLOSER

        return table_line+"\n"
    



class SampleTable:

    HEADER_VALUES=["Amplicon","Tot. Reads", "Rejected on Length", "% Target Org.",
                        "# Target Org.","All known SNPs", "Dominant known SNPs",
                        "Dominant unknown SNPs", "Implied Genotypes", "Implied AMR"]
    HEADER_TOOLTIP_TEXT=["Sample name, you can make it more informative by adding sample description using options -c and -d.",
                         "Total number of reads that mapped to the amplicon regadless of length (i.e. -s and -l parameters)",
                         "Number of mapped reads that were not assesed by classifier due to not meeting -s and -l parameters.",
                         "Percentage of mapped reads that were classified as target organism and that satisfied -s and -l length parameteres",
                         "Number of reads that were classified as target organism and that satisfied -s and -l length parameteres",
                         f'Number of SNPs detected that match known SNPs defined in classifier model file you used.&#10\
Includes all alleles different from references that occur in >{POSITIVE_CASES_THRESHOLD} reads and no less than {GT_REPORTING_MIN_READ_PERCENT}% of all reads',
                         f'Shows those detected SNPs that match known SNPs in classifier model file you used and also are the most common nucleotide that position.&#10\
E.g. if at position 254 the reference is C, but in reads the most common nucleotide is T this would count in this column.&#10\
If T occcured, but the most common was C, this would not count in this column. ',
                          f'Shows those detected SNPs that DO NOT match known SNPs in classifier model file you used and also are the most common nucleotide that position.&#10\
E.g. if at position 315 the reference is A, but in reads the most common nucleotide is G this would count in this column.&#10\
If G occcured, but the most common was A, this would not count in this column. If position has three occuring nucleotides, only the most common is examined.',
                         "Genotypes implied by detected SNPs (but not reference alleles). Number in brackets indicates what % of target organism reads supports this genotype.&#10\
For high level genotypes, a range of possible gnotypes (e.g. 0-2 or 3-4) is reported.",
                         "AMR genotype implied by detected SNPs (but not reference alleles). Number in brackets indicates what % of target organism reads supports this genotype."]
                         
    def __init__(self) -> None:
        self._utilities=ReportingUtilities()
        self._reporting_genotype=ReportingGenotype()

    def get_table_header(self, sample: SampleReportValues) -> str:
        text=self._utilities.insert_paragraph(2)
        text+='<div class=header_line><a name="'+sample.name+'">'+sample.name_with_description+'</a><div><a href="#Summary">Back to Summary</a></div></div>\n'
        text+=self._utilities.insert_paragraph(1)
        text+="<table>\n"
        text+="\t<tbody>\n"
        text+=self._utilities.generate_header_row(SampleTable.HEADER_VALUES, SampleTable.HEADER_TOOLTIP_TEXT)
        text+="\n"

        return text

    def get_table_row(self, amplicon:AmpliconReportValues, sample: SampleReportValues, alt_style=False) -> str:
        first_cell='<a href="#'+sample.name+'_'+amplicon.name+'_consensus"/>'+amplicon.name
        first_cell+="^" if amplicon.model.not_trained else ""
        if amplicon.target_org_reads<POSITIVE_CASES_THRESHOLD:
            values=[first_cell]+["*"]+["-"]*8
        else:
            # high_level_gts=self._reporting_genotype.possible_high_level_gts( amplicon.alleles, sample.model_data.metadata["high_level_genotypes"] )
            # if high_level_gts==ReportingGenotype.DEFAULT_HIGH_LEVEL_GT:
            #     high_level_gts==""
            values=[first_cell, amplicon.mapped_count, amplicon.wrong_length_count, amplicon.target_org_share,
                        amplicon.target_org_reads, len(amplicon.all_known_snps), len(amplicon.dominant_known_snps),
                        len(amplicon.dominant_unknown_snps), amplicon.amplicon_gts_str_format, amplicon.amplicon_amrs]
            
        if alt_style: #Alternate row colours
            table_line=OPENER_ALT+MIDDLE.join([str(f) for f in values])+CLOSER
        else:
            table_line=OPENER+MIDDLE.join([str(f) for f in values])+CLOSER

        return table_line+"\n"


    def get_table_tail(self) -> str:
        text="\t<tbody>\n"
        text+="<table>\n"
        text+=f'* Amplicon has fewer than {POSITIVE_CASES_THRESHOLD+1} target organism reads<br>'
        text+=f'<p>(?) Prediction is based on amplicon with fewer than {POSITIVE_CASES_THRESHOLD} reads</p>'
        #text+="<p><del>AMR Amplicon</del> means no reads were found for this AMR amplicon</p>"
        text+='^ Using model that was trained without negative cases - results are likely to be less reliable'
        text+=self._utilities.insert_paragraph(5)

        text+="</div>\n"
        return text

class ConsensusSequencesTable:
    def __init__(self) -> None:
        self.utilities=ReportingUtilities()


    FIRST_LINE='<div class=header_line>Consensus sequences</div>\n'
    HEADER_VALUES=["Sample", "Amplicon", "# of target reads", "Consensus <br/> Mismatches <br/> vs Reference", "Consensus sequence"]
    HEADER_TOOLTIP_TEXT=["Sample name, you can make it more informative by adding sample description using options -c and -d.",
                         "Amplicon name",
                         "Number of reads classified as target organism",
                         "List of position where most common nucleotide does not match the reference sequence. If allele has known consequence this is also listed.&#10\
This list usually doesn't match the one in Sample's table because this column only reports cases where most common nucleotide doesn't match reference.",
                         "Amplicon consensus sequence. If the most common nucleotide is deletion this is reported as '-'. Insertions are not currently supported."]


    def get_main_header(self) -> str:
        main_header=ConsensusSequencesTable.FIRST_LINE
        main_header+=self.utilities.insert_paragraph(1)
        main_header+="<table>\n"
        main_header+="\t<tbody>\n"
        main_header+=self.utilities.generate_header_row(ConsensusSequencesTable.HEADER_VALUES, ConsensusSequencesTable.HEADER_TOOLTIP_TEXT)
        return main_header
    
    def get_table_row(self, amplicon: AmpliconReportValues, sample: SampleReportValues, alt_style=False) -> str:
        hyperlink_name=sample.name+"_"+amplicon.name+"_consensus"
        first_cell='<a name="'+hyperlink_name+'" href="#'+sample.name+'">'+sample.name_with_description+'</a>'
        values=[first_cell, amplicon.name, amplicon.target_org_reads, amplicon.consensus_snp_implications, amplicon.consensus]
        if alt_style: #Alternate row colours
            table_line=OPENER_ALT+MIDDLE.join([str(f) for f in values])+CLOSER
        else:
            table_line=OPENER+MIDDLE.join([str(f) for f in values])+CLOSER

        return table_line+"\n"
    
    def get_main_tail(self) -> str:
        main_tail="\t</tbody>\n"
        main_tail+="</table>\n"
        main_tail+="</div>\n"
        return main_tail

class DelimitedTable:
    def __init__(self):
        pass

    def get_main_header(self):
        values=["Barcode", "Amplicon", "Amplicon Description", "Tot. Reads", "Rejected on Length", "% Target Org.",
                 "# Target Org.", "Possible High Level Genotypes", "Implied Genotypes",
                   "Implied AMR", "Dominant known SNPs", "Dominant unknown SNPs",
                     "Consensus Mismatches vs Reference", "Consensus"]
        return "\t".join(values)+"\n"
    
    def get_table_row(self, amplicon:AmpliconReportValues, sample: SampleReportValues):
        values=[sample.name, amplicon.name, sample.sample_description, amplicon.mapped_count,
                amplicon.wrong_length_count, amplicon.target_org_share, amplicon.target_org_reads,
                sample.possible_highlevel_gts, amplicon.amplicon_gts_str_format, amplicon.amplicon_amrs,
                len(amplicon.dominant_known_snps), len(amplicon.dominant_unknown_snps),
                amplicon.consensus_snp_implications,  amplicon.consensus]
        
        table_line = "\t".join([str(f) for f in values])
            
        return table_line+"\n"