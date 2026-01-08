import pickle
from typing import List, Dict
from read_classifier import ClassificationResult, Classifier, ModelsData
from os.path import expanduser
from map import MappingResult
from tqdm import tqdm
from reporting_classes import SummaryTable, ReportingUtilities, SampleReportValues, SampleTable, ConsensusSequencesTable, DelimitedTable, ReportingGenotype

OPENER='<tr><td>'
OPENER_BOLD='<tr class=table_header><td>'
OPENER_ALT='<tr class=alt_row><td>'
CLOSER='</td></tr>'
MIDDLE='</td><td>'
MIDDLE_WIDE='</td class=td_alt_wide><td>'

class ClasssifierReport:

    def __init__(self, output_file: str, models_file: str, positive_amplicons_threshold: float, length_parameters: Dict[str, int],
                 sample_labels:Dict[str, str]={}, mapping_results: List[MappingResult]=[]) -> None:
        self._positive_amplicons_threshold = positive_amplicons_threshold/100
        self._report_utilities=ReportingUtilities()
        self._model_file=models_file
        self._output_file_name=output_file
        self.output_file=open( expanduser(output_file), "w" )
        with open(expanduser(models_file), "rb") as model_file:
            self.models_data: ModelsData = pickle.load(model_file)
        self.sample_labels=sample_labels
        self.sample_hyperlinks={}
        self.mapping_results={}
        self._length_parameters=length_parameters #the inputs values for -s and -l.
        for item in mapping_results:
            self.mapping_results[item.sample_name] = item.total_reads


    ####START Write models' summary
    def _write_model_summary(self, results) -> str:
        model_summary='<div class=header_line><a name="ModelSigna ures">Model Signatures</a><div><a href="#Summary">Back to Summary</a></div></div>\n'
        model_summary+=self._report_utilities.insert_paragraph(1)
        model_summary+="<table>\n"
        model_summary+="\t<tbody>\n"
        model_summary+=self._report_utilities.generate_header_row(["Amplicon", "Is Core","Signature", "Date Trained"])
        models_list=sorted(set([(f.amplicon.name, f.model_fingerprint,f.model_timestamp) for f in results]))
        for i, result in enumerate(models_list):
            is_core=not self.models_data.classifiers[result[0]].amplicon_transient
            row_values=[result[0], is_core, result[1], result[2].strftime("%d/%b/%y at %H:%M")]
            if i % 2 == 0: #Alternate row colours
                table_row=OPENER_ALT+MIDDLE.join([str(f) for f in row_values])+CLOSER
            else:
                table_row=OPENER+MIDDLE.join([str(f) for f in row_values])+CLOSER
            model_summary+="\t\t"+table_row+"\n"
        model_summary+="<table>\n"
        model_summary+="\t<tbody>\n"
        model_summary+="</div>\n"
        return model_summary
    ####END Write models' summary

    def create_report(self, results: List[ClassificationResult]):
        
        samples=set([f.sample for f in results])
        sample_report_values: Dict[str, SampleReportValues ]={}
        for sample in tqdm(samples):
            sample_report_values[sample]=SampleReportValues(sample, [ f for f in results if f.sample==sample ], self.models_data )
            if sample_report_values[sample].name in self.sample_labels:
                sample_report_values[sample].sample_description=self.sample_labels[sample_report_values[sample].name]
            if sample_report_values[sample].name in self.mapping_results:
                sample_report_values[sample].sequencing_read_count=self.mapping_results[sample_report_values[sample].name]

        self.output_file.write("<html>\n")
        self.output_file.write("<body>\n")

        summary_table=SummaryTable( self._report_utilities.get_most_recent_model_element(results),
                                    self._model_file, "high_level_genotypes" in self.models_data.metadata)
        self.output_file.write( self._report_utilities.get_styles() )
        summary_table.calculate_error_rate(results)
        self.output_file.write( summary_table.get_main_header() )
        #sort the list of samples and amplicons to ensure consistent order of results
        sample_report_values=sorted(sample_report_values.items(), key= lambda x: x[0])

        i=0
        for sample_name, sample in  sample_report_values:
            self.output_file.write(  summary_table.get_table_row(sample, i % 2 ==0) )
            i+=1
            
        self.output_file.write( summary_table.get_main_tail(self._length_parameters) )

        # self._write_gt_suppport()

        sample_table=SampleTable()
        for sample_name, sample in sample_report_values:
            self.output_file.write(sample_table.get_table_header(sample))
            i=0
            for amplicon_name, amplicon in sample.display_order:
                self.output_file.write(  sample_table.get_table_row(amplicon, sample, i % 2 ==0) )
                i+=1
            self.output_file.write(sample_table.get_table_tail())

        consensus_table=ConsensusSequencesTable()
        self.output_file.write(consensus_table.get_main_header())
        i=0
        for sample_name, sample in sample_report_values:
            for amplicon_name, amplicon in sample.display_order:
                self.output_file.write(  consensus_table.get_table_row(amplicon, sample, i % 2 ==0) )
                i+=1
        self.output_file.write(consensus_table.get_main_tail())

        self.output_file.write(self._write_model_summary(results))

        self.output_file.write("</body>\n")
        self.output_file.write("</html>\n")
        self.output_file.close()

        with open(expanduser(self._output_file_name.replace(".html", ".tsv")), "w") as delim_data:
            delimited_table=DelimitedTable()
            delim_data.write(delimited_table.get_main_header())
            for sample_name, sample in sample_report_values:
                for amplicon_name, amplicon in sample.display_order:
                    delim_data.write(  delimited_table.get_table_row(amplicon, sample) )

        with open(expanduser(self._output_file_name.replace(".html", "_summary.tsv")), "w") as summary_tsv:
            summary_tsv.write("\t".join(summary_table.HEADER_VALUES)+"\n")
            for sample_name, sample in  sample_report_values:
                row_value=sample.name+"\t"+ "\t".join([str(f) for f in summary_table.assemble_row_values(sample)] )+"\n"
                summary_tsv.write( row_value.replace("&nbsp;","")  )
