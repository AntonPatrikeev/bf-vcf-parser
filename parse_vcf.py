'''

'''

import argparse
from argparse import RawTextHelpFormatter
from pkgs.vcf_report_generator import VCFReportGenerator
from pkgs.vcf_parser import VCFParser



def main(vcf_path, outdir):
    # parse VCF
    vcf_parser = VCFParser(vcf_path)
    vcf_parser.parse_vcf()

    # generate report
    vcf_report_generator = VCFReportGenerator(vcf_parser, outdir)
    vcf_report_generator.create_docx_report()
    vcf_report_generator.output_metadata_metrics()
    vcf_report_generator.output_variants_metrics()
    vcf_report_generator.output_samples_metrics()
    vcf_report_generator.save_report()




if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description=\
        '''
        This script parses the input VCF file and outputs .docx report with a summary on VCF contents.

        Requirements to the input VCF file:
        1/ must have fileformat field in the metadata (otherwise is considered a non-VCF file)
        2/ must contain only unique biallelic variants based on hg38 reference sequence
        3/ all genotypes must be called: missing genotypes (./.) are not allowed; FORMAT consists of GT only
        4/ it is assumed that VT INFO field contain variant types (SNP / INDEL)

        Output report file is a .docx document with the same base name as the input VCF and .report.docx extension.
        ''')

    parser.add_argument("--vcf", type=str, help="path to the input VCF", required=True)
    parser.add_argument("--out", type=str, help="path to output directory", required=True)

    args = parser.parse_args()

    vcf_path = args.vcf
    outdir = args.out

    main(vcf_path, outdir)
