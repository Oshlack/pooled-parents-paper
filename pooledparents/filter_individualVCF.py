import argparse
import sys
import vcf
from vcf.parser import _Filter
import os
from collections import Counter
# shared functions
from filtervcf import *

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description=('Filter proband vcfs based on'
        ' a pooled parent vcf and report recall per variant'))
    parser.add_argument(
        '--individual_vcfs', type=str, required=True, nargs='+',
        help='VCFs containing variant calls from all probands (one for each individual)')
    parser.add_argument(
        '--pool_vcf', type=str, required=True,
        help='A single VCF containing variant calls for the pool')
    parser.add_argument(
        '--out_csv', type=str, required=False, default='pools_probands.compare.csv',
        help='Output filename for csv (default: %(default)s)')
    parser.add_argument(
        '--suffix', type=str, required=False, default='.pool_filter.vcf',
        help=('Output suffix for filtered proband VCFs. Filenames will be in '
            'the form inferred_sample_ID + suffix (default: %(default)s)'))
    parser.add_argument(
        '--falsepos', action='store_true',
        help=('Report false positives as additional lines in the output csv '
            '(only relevant for comparing pooled vs. individual sequencing of '
            'the same individuals).'))
    parser.add_argument(
        '--exclude_filtered', action='store_true',
        help=('Do not output filtered variants to the individual VCFs. By '
            'default filtered variants are reported with the value "inPool" in '
            'the FILTER field.'))
    return parser.parse_args()

def parse_pool_vcf(pool_vcf_file):
    """Parse pool VCF and save all variants found in them
    Args:
        pool_vcf_file (str): path to pool VCF file
    Returns:
        set: variant ids for all variants in the pool VCF
    """
    pool_vars = set()

    pool = sample_id_from_fname(pool_vcf_file)

    with open(pool_vcf_file, 'r') as this_vcf:
        for record in vcf.Reader(this_vcf):
            variants = variant_id_split(record)
            for variant in variants: # usually one, but could be multiple
                pool_vars.add(variant)
    return pool_vars

def R_bool(py_bool):
    """Convert python bool (True/False) to a string that can be read as a bool
    by R (TRUE/FALSE/NA)"""
    if py_bool == None:
        return 'NA'
    if py_bool:
        return 'TRUE'
    elif not py_bool:
        return 'FALSE'
    else:
        return 'NA'

def main():
    # Parse command line arguments
    args = parse_args()
    individual_vcf_files = args.individual_vcfs
    pool_vcf_file = args.pool_vcf
    outfile = args.out_csv
    out_vcf_suffix = args.suffix
    report_falsepos = args.falsepos
    output_filtered = not args.exclude_filtered

    outstream = open(outfile, 'w')

    # Write header
    outstream.write(('proband,variant,recovered_proband,falsepos,'
                    'QD,AF_EXOMESgnomad'
                    '\n'))

    # Parse vcfs for pools
    # Simply record which variants were found in the pool
    pool_vars = parse_pool_vcf(pool_vcf_file)

    # Parse vcfs of individuals
    individual_vars = set()
    for vcf_file in individual_vcf_files:
        proband = sample_id_from_fname(vcf_file)
        with open(vcf_file, 'r') as this_vcf:
            vcf_reader = vcf.Reader(this_vcf)
            # Add an aditional filter that will be inherited by the vcf writer
            vcf_reader.filters['InPool'] = _Filter('InPool',
                'All alleles found in the probands are also found in the pool.')
            # Create vcf writer based on the header from the input vcf
            vcf_writer = vcf.Writer(open(proband + out_vcf_suffix, 'w'), vcf_reader)

            for record in vcf_reader:
                falsepos = 'FALSE'
                qual = record.QUAL
                QD = qual/record.INFO['DP']

                variants = variant_id_split(record)
                AF_EXOMESgnomad_all = extract_record_info_multi(record, 'AF_EXOMESgnomad')
                if len(AF_EXOMESgnomad_all) != len(variants):
                    if not AF_EXOMESgnomad_all == ['NA']:
                        sys.stderr.write(('WARNING: Number of variant allelese and gnomAD '
                            'records do not match. Writing AF_EXOMESgnomad = NA '
                            'for all. Variants: {} AF_EXOMESgnomad {} \n'
                            ).format(variants, AF_EXOMESgnomad_all))
                    AF_EXOMESgnomad_all = ['NA'] * len(variants)

                all_variants_in_pool = True
                for variant, AF_EXOMESgnomad in zip(variants, AF_EXOMESgnomad_all):
                    individual_vars.add(variant)
                    variant_in_pool = variant in pool_vars
                    # If any variant is not in the pool, then set to false
                    if not variant_in_pool:
                        all_variants_in_pool = False

                    variant_in_pool = R_bool(variant_in_pool)
                    outstream.write(','.join([str(x) for x in [
                        proband, variant, variant_in_pool, falsepos,
                        QD, AF_EXOMESgnomad
                        ]]) + '\n')

                # Either report variants as filtered in the VCF or skip them completely
                if all_variants_in_pool:
                    if output_filtered:
                        record.FILTER = 'InPool' # Set in_pool vcf filter
                        vcf_writer.write_record(record)
                else:
                    vcf_writer.write_record(record)

    if report_falsepos:
        proband = 'NA'
        variant_in_pool = 'TRUE'
        falsepos = 'TRUE'
        QD = 'NA'
        AF_EXOMESgnomad = 'NA'
        with open(pool_vcf_file, 'r') as this_vcf:
            for record in vcf.Reader(this_vcf):
                variants = variant_id_split(record)
                for variant in variants: # usually one, but could be multiple
                    if not variant in individual_vars:
                        outstream.write(','.join([str(x) for x in [
                            proband, variant, variant_in_pool, falsepos,
                            QD, AF_EXOMESgnomad
                            ]]) + '\n')


    outstream.close

if __name__ == '__main__':
    main()
