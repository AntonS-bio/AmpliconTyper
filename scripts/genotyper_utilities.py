from read_classifier import Classifier, ModelsData, GenotypeSNP
from pickle import load, dump
from typing import Dict
from os.path import expanduser, exists
import argparse
import warnings
import pandas as pd
warnings.filterwarnings("ignore")

def get_reference_fasta(model_file: str) -> None:
    with open(model_file, "rb") as input_model:
        model_manager: Dict[str, Classifier] =load(input_model)
        for name, model in model_manager.classifiers.items():
            print(">"+str(model.name))
            print(str(model.nucletoide_seq))

def get_snps(model_file: str) -> None:
    with open(model_file, "rb") as input_model:
        model_manager: Dict[str, Classifier] =load(input_model)
        for name, model in model_manager.classifiers.items():
            for genotype_snp in model.genotype_snps:
                print(genotype_snp.to_string)

def load_hierarchy(model_file: str, heirarchy_file:str) -> None:
    with open(model_file, "rb") as input_model:
        model_manager: Dict[str, Classifier] =load(input_model)



def test(model_file: str) -> None:
    pass
    # with open("/home/lshas17/AmpliconTyper/models/paratyphi_A_v3.pkl", "rb") as input_model:
    #     model_manager_para: Dict[str, Classifier] =load(input_model)

    # with open("/home/lshas17/paratyphi/model/gt_snps.tsv") as snps:
    #     snps.readline()
    #     for line in snps:
    #         contig_id, reference_nucl, Alt_nucl, GT, position, is_amr=line.strip().split("\t")
    #         if GT=="2.3":
    #             snp=GenotypeSNP(contig_id, position,Alt_nucl,False)
    #             snp.genotypes[reference_nucl]=GT
    #         else:
    #             snp=GenotypeSNP(contig_id, position,reference_nucl,False)
    #             snp.genotypes[Alt_nucl]=GT
    #         model_manager_para.classifiers[contig_id].genotype_snps.append(snp)

    # new_model_file=expanduser("/home/lshas17/AmpliconTyper/models/paratyphi_A_v4.pkl")
    # with open(new_model_file, "wb") as new_model_file_output:
    #     dump(model_manager_para, new_model_file_output)

def rename_mdr_loci(model_manager: ModelsData) -> None:
    new_names={"chr_4.3.1.1_none_LT904852.1":"chromosomal_MDR_yidA",
                    "chr_4.3.1.1_none_LT904894.1":"chromosomal_MDR_cyaA",
                    "plasmid_2.2_none_LT904892.1":"plasmid_MDR_LT904892.1",
                    "plasmid_3.2.1_non-PST6_AL513383.1":"plasmid_MDR_AL513383.1",
                    "plasmid_4.3.1.1_PST6_CP029645.1":"plasmid_MDR_CP029645.1",
                    "plasmid_4.3.1.1_PST6_LT904879.1":"plasmid_MDR_LT904879.1",
                    "plasmid_4.3.1.3_PST6_CP029879.1":"plasmid_MDR_CP029879.1",
                    "plasmid_4.3.1.3_PST6_CP029924.1":"plasmid_MDR_CP029924.1",
                    "plasmid_4.3.1.3_PST6_CP029957.1":"plasmid_MDR_CP029957.1",
                    }
    # for key, value in new_names.items():
    #     model_manager.classifiers[value]=model_manager.classifiers[key]
    #     model_manager.classifiers[value].name=value
    #     for snp in model_manager.classifiers[value].genotype_snps:
    #         snp.contig_id=value
    #     model_manager.classifiers.pop(key)

    # new_model_file=expanduser("~/HandyReadGenotyper/models/typhi_v7.pkl")
    # with open(new_model_file, "wb") as new_model_file_output:
    #     dump(model_manager, new_model_file_output)


def main():

    parser = argparse.ArgumentParser(description='Various utilities for manipulating classification model')
    parser.add_argument('--reference', action="store_true",
                        help='Use to get reference sequences used in classification model. Must specify model (-m)', required=False)
    parser.add_argument('--snps', action="store_true",
                        help='Use to get SNPs in reference sequences used in genotyping. Must specify model (-m)', required=False)
    parser.add_argument('--test', action="store_true",
                        help='Used for testing only', required=False)
    parser.add_argument('-m','--model', metavar='', type=str,
                        help='Pickle (.pkl) file containing pretrained model.', required=False)

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        exit(0)

    if args.model is None:
        print("To get reference sequence a model file must be specified!")
        exit(0)
    model_file=expanduser(args.model)
    if not exists(model_file):
        print(f'Model file {model_file} does not exist!')
        exit(0)

    if args.test:
        test(model_file)
    if args.reference or args.snps:
        if args.reference:
            get_reference_fasta(model_file)
        elif args.snps:
            get_snps(model_file)

if __name__=="__main__":
    main()