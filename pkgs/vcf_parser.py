import pandas as pd
import numpy as np
import os
from pathlib import Path

from .vcf_parser_utils import convert_variants_to_df_by_mb, convert_sample_variants_to_df_by_mb

SUBSTITUTION_TYPES_MAPPING = {
        "C>A":"C:G > A:T", "C>G":"C:G > G:C", "C>T":"C:G > T:A", "T>A":"T:A > A:T", "T>C":"T:A > C:G", "T>G":"T:A > G:C",
        "G>T":"C:G > A:T", "G>C":"C:G > G:C", "G>A":"C:G > T:A", "A>T":"T:A > A:T", "A>G":"T:A > C:G", "A>C":"T:A > G:C"
    }


class VCFParser:
    ''' Basic VCF file parser class '''
    def __init__(self, vcf_path=None):
        if not os.path.exists(vcf_path):
            raise OSError(f'VCF file not found: {vcf_path}')
        else:
            self.vcf_path = vcf_path
        self.metadata = None
        self.variants = None
        self.samples = None
    
    def parse_metadata(self):
        ''' 
        Method parses metadata lines of the input VCF file line by line and updates metadata attribute.
        metadata is a dict
        '''
        metadata = {"fileformat":None, "reference":"unknown", "metadata_size":None, "header_line":None, "samples":None}
        with open(self.vcf_path, 'r') as vcfin:
            for line_idx, line in enumerate(vcfin):
                # parse metadata lines
                if line[:2] == "##":
                    metadata_key, metadata_value = line.strip()[2:].split("=", 1)
                    if metadata_key in metadata:
                        metadata[metadata_key] = metadata_value
                # parse header line
                elif line[0] == "#":
                    metadata["header_line"] = line.strip().split('\t')
                    metadata["samples"] = metadata["header_line"][9:]
                    metadata["metadata_size"] = line_idx
                else:
                    if not metadata["fileformat"]:
                        raise UnknownVCFFormatError("Failed to detect VCF format in the metadata: fileformat is required in VCF")
                    else:
                        self.metadata = metadata


    def parse_vcf(self):
        ''' 
        Method parses the input VCF file line by line and updates variants and samples attributes.
        variants is a nested dict with the following structure: {contigs: {positions: [ref, alt, vt]}}
        samples is a nested dict with the following structure: {samples: {contigs: {positions: gt}}}
        '''
        if not self.metadata:
            self.parse_metadata()

        variants = dict()
        samples = dict((sample, {}) for sample in self.metadata["samples"])

        with open(self.vcf_path, 'r') as vcfin:
            for line_idx, line in enumerate(vcfin):
                if line_idx > self.metadata["metadata_size"]:
                    # parse vcf record
                    record_parser = VCFRecordParser(line, line_idx,
                                                  self.metadata["header_line"], 
                                                  self.metadata["samples"])
                    # get variant features
                    contig = record_parser.get_contig()
                    pos = record_parser.get_pos()
                    ref = record_parser.get_ref_allele()
                    alt = record_parser.get_alt_allele()
                    variant_type = record_parser.get_variant_type()

                    # add contig if it is new
                    if contig not in variants:
                        variants[contig] = {}
                        for sample in self.metadata["samples"]:
                            samples[sample][contig] = {}

                    # raise error if variant already exists in variants dict (and thus in samples dict) 
                    if pos in variants[contig]:
                        raise VariantDuplicationError(f'At least two variants with POS {pos} on chromosome {contig} found.')
                    # update variants and samples dicts if the variant is new
                    else:
                        variants[contig][pos] = [ref, alt, variant_type]
                        for sample in self.metadata["samples"]:
                            samples[sample][contig][pos] = record_parser.get_genotype(sample)

        self.variants = variants
        self.samples = samples

    def get_samples_number(self):
        return len(self.metadata["samples"])

    def get_contigs_number(self):
        if self.variants:
            return len(self.variants)
        else:
            raise Exception("Can't get contigs number: VCF is not parsed yet.")

    def get_variants_number(self):
        if self.variants:
            return sum([len(self.variants[contig]) for contig in self.variants])
        else:
            raise Exception("Can't get variants number: VCF is not parsed yet.")

    def get_variant_type_counts_per_mb(self, contig):
        ''' 
        Method is used for processing of parsed VCF 
        for further creation of plots in the Variant summary section
        Returns DataFrame with the following columns:
        window:float, SNP:float, DEL:float, INS:float
        '''
        # convert contig variants dict to dataframe
        variants_contig_df = convert_variants_to_df_by_mb(self.variants[contig])

        # specify variant type: separate INDELs into INSs and DELs
        variants_contig_df["VT"] = np.where(variants_contig_df.VT != "INDEL", 
                                            variants_contig_df.VT, 
                                            np.where(variants_contig_df.REF.map(len) > variants_contig_df.ALT.map(len), 
                                                "DEL", 
                                                "INS"
                                                )
                                            )
        # group by window and variant type and pivot dataframe so that variant types are in columns
        variants_contig_df = pd.pivot_table(
                                    variants_contig_df.groupby(["window", "VT"], as_index=False)["POS"].count(), 
                                    values='POS', 
                                    index=['window'],
                                    columns=['VT'], 
                                    fill_value=0
                                        ).reset_index()

        return variants_contig_df

    def get_substitution_type_prc_per_mb(self, contig):
        ''' 
        Method is used for processing of parsed VCF 
        for further creation of plots in the Variant summary section.
        Returns DataFrame with the following columns:
        window:float, C:G > A:T:str... and other 5 substitution types
        '''
        # convert contig variants dict to dataframe
        substitutions_contig_df = convert_variants_to_df_by_mb(self.variants[contig])

        # keep only SNPs
        substitutions_contig_df = substitutions_contig_df[substitutions_contig_df.VT == "SNP"]

        # map substitution type
        substitutions_contig_df["substitution_type"] = substitutions_contig_df.apply(lambda x: SUBSTITUTION_TYPES_MAPPING[f"{x.REF}>{x.ALT}"], axis=1)
        
        # group by window and variant type
        substitutions_contig_df = substitutions_contig_df.groupby(["window", "substitution_type"], as_index=False)['POS'].count()

        # convert variants counts to percentages
        substitutions_contig_df["POS"] = substitutions_contig_df["POS"] / substitutions_contig_df.groupby("window")["POS"].transform('sum') * 100

        # pivot dataframe so that substitution types are in columns
        substitutions_contig_df = pd.pivot_table(
                                            substitutions_contig_df, 
                                            values='POS', 
                                            index=['window'],
                                            columns=['substitution_type'], 
                                            fill_value=0
                                                ).reset_index()

        return substitutions_contig_df

    def get_sample_variant_type_counts_per_mb(self, sample, contig):
        ''' 
        Method is used for processing of parsed VCF 
        for further creation of plots in the Samples variant profiles summary section.
        Returns DataFrame with the following columns:
        window:float, VT:str, HET:float, HOM:float
        '''
        # convert contig variants dict to dataframe
        sample_variants_contig_df = convert_sample_variants_to_df_by_mb(self.samples[sample][contig], self.variants[contig])

        # specify genotype: keep only HOM = HOMA (homozygous by alternative allele) and HET (heterozygous)
        sample_variants_contig_df = sample_variants_contig_df[sample_variants_contig_df.GT.isin(["HOMA", "HET"])]
        sample_variants_contig_df["GT"] = sample_variants_contig_df["GT"].apply(lambda x: x[:3])

        # specify variant type: separate INDELs into INSs and DELs
        sample_variants_contig_df["VT"] = np.where(sample_variants_contig_df.VT != "INDEL", 
                                            sample_variants_contig_df.VT, 
                                            np.where(sample_variants_contig_df.REF.map(len) > sample_variants_contig_df.ALT.map(len), 
                                                "DEL", 
                                                "INS"
                                                )
                                            )

        # group by window, variant type and genotype and pivot dataframe so that genotypes (HOM / HET) are in columns
        sample_variants_contig_df = pd.pivot_table(
                                    sample_variants_contig_df.groupby(["window", "VT", "GT"], as_index=False)["POS"].count(), 
                                    values='POS', 
                                    index=['window', 'VT'],
                                    columns=['GT'], 
                                    fill_value=0
                                        ).reset_index()

        return sample_variants_contig_df

    def get_sample_substitution_type_prc_per_mb(self, sample, contig):
        ''' 
        Method is used for processing of parsed VCF 
        for further creation of plots in the Samples variant profiles summary section.
        Returns DataFrame with the following columns:
        window:float, C:G > A:T:str... and other 5 substitution types
        '''
        # convert contig variants dict to dataframe
        sample_substitutions_contig_df = convert_sample_variants_to_df_by_mb(self.samples[sample][contig], self.variants[contig])

        # keep only SNPs
        sample_substitutions_contig_df = sample_substitutions_contig_df[sample_substitutions_contig_df.VT == "SNP"]

        # map substitution type
        sample_substitutions_contig_df["substitution_type"] = sample_substitutions_contig_df.apply(lambda x: SUBSTITUTION_TYPES_MAPPING[f"{x.REF}>{x.ALT}"], axis=1)
        
        # group by window and variant type
        sample_substitutions_contig_df = sample_substitutions_contig_df.groupby(["window", "substitution_type"], as_index=False)['POS'].count()

        # convert variants counts to percentages
        sample_substitutions_contig_df["POS"] = sample_substitutions_contig_df["POS"] / sample_substitutions_contig_df.groupby("window")["POS"].transform('sum') * 100

        # pivot dataframe so that substitution types are in columns
        sample_substitutions_contig_df = pd.pivot_table(
                                            sample_substitutions_contig_df, 
                                            values='POS', 
                                            index=['window'],
                                            columns=['substitution_type'], 
                                            fill_value=0
                                                ).reset_index()

        return sample_substitutions_contig_df


class VCFRecordParser:
    ''' Type for parsing variant record from a VCF file.'''
    def __init__(self, record=None, record_idx=None, colnames=None, samples=None):
        self.record = record
        self.record_idx = record_idx
        self.colnames = colnames
        self.samples = samples
        self.parsed_record = self.parse_record()

    def parse_record(self):
        # check that the number of items in the record equals the number of VCF colums 
        # and split record and convert it to dict with column names as keys
        record_items = self.record.strip().split("\t")
        if len(self.colnames) != len(record_items):
            raise VCFRecordTokenizingError(f'Number of items in the VCF record does not equal to the number of items in the header line (see line {self.record_idx + 1})')
        else:
            parsed_record = dict(zip(self.colnames, record_items))

        return parsed_record

    def get_contig(self):
        contig = self.parsed_record["#CHROM"]
        if contig[:3] != "chr":
            return f'chr{contig}'
        else:
            return str(contig)

    def get_pos(self):
        pos = self.parsed_record["POS"]
        if not pos.isdigit():
            raise TypeError(f'Incorrect POS type: POS value must be numeric, but "{pos}" found (see line {self.record_idx + 1})')
        else:
            return int(pos)

    def get_ref_allele(self):
        ref = self.parsed_record["REF"]
        if "," in ref:
            raise ValueError(f'Variant is not biallelic: multiple REF alleles found (see line {self.record_idx + 1})')
        else:
            return ref

    def get_alt_allele(self):
        alt = self.parsed_record["ALT"]
        if "," in alt:
            raise ValueError(f'Variant is not biallelic: multiple ALT alleles found (see line {self.record_idx + 1})')
        else:
            return alt
    
    def get_variant_type(self):
        ''' Returns value in VT key in the INFO field or UNKNOWN if VT not found '''
        info = self.parsed_record["INFO"]
        info_dict = dict([info_metric.split("=") for info_metric in info.split(";")])
        return info_dict.setdefault("VT", "UNKNOWN")

    def get_genotype(self, sample):
        ''' 
        Returns HOMR for homozygous variant with reference alleles,
        HOMA for homozygous variant with alternative alleles,
        HET for heterozygous variant.    
        '''
        gt = self.parsed_record[sample]

        ref_allele_count = gt.count("0")
        alt_allele_count = gt.count("1")
        if ref_allele_count + alt_allele_count != 2:
            raise ValueError(f'Incorrect genotype value: must be a combination of 0 and 1, but "{gt}" found (see line {self.record_idx + 1})')
        else:
            if ref_allele_count == 2:
                return "HOMR"
            elif alt_allele_count == 2:
                return "HOMA"
            else:
                return "HET"


class UnknownVCFFormatError(Exception):
    ''' Exception class for unknown VCF file format error '''
    pass

class VCFRecordTokenizingError(Exception):
    ''' Exception class for VCF variant record tokenization error '''
    pass

class VariantDuplicationError(Exception):
    ''' Exception class for VCF variant duplication error '''
    pass