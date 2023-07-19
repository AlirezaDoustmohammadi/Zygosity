import argparse
import os
import pandas as pd
import filteration


def calc_zygosity(df,  vcf_file_path, include_info, output_format, output_file, sample_name=None, verbosity=False):
    """
    Calculates the zygosity for the variants of a specific sample in a VCF file.

    @param df Pandas DataFrame representing a VCF file.
    @param vcf_file_path The path of the original VCF file.
    @param include_info A comma-separated string of additional INFO fields to include in the output.
    @param output_format The format of the output file ('csv' or 'json').
    @param output_file Boolean indicating whether to output a file.
    @param sample_name The name of the sample for which zygosity will be calculated.
    @param verbosity A boolean flag that when True, prints out diagnostic messages to console.

    """

    variants = {}

    if sample_name:
        # Get column names that start with 'SAMPLE' but not equal to sample_name
        cols_to_drop = df.columns[df.columns.str.startswith('SAMPLE') & (df.columns != sample_name)]
        # Drop these columns
        df = df.drop(columns=cols_to_drop)

        # Apply the lambda function to count '1/1', '1/0', '0/1', '0/0' and sum the counts
        variants['g 1/1'] = df[[sample_name]].applymap(lambda x: '1/1' in x.split(':')[0]).sum().sum()
        variants['g 1/0'] = df[[sample_name]].applymap(lambda x: '1/0' in x.split(':')[0]).sum().sum()
        variants['g 0/1'] = df[[sample_name]].applymap(lambda x: '0/1' in x.split(':')[0]).sum().sum()
        variants['g 0/0'] = df[[sample_name]].applymap(lambda x: '0/0' in x.split(':')[0]).sum().sum()
    else:
        # Get all columns starting with 'SAMPLE'
        sample_columns = [col for col in df.columns if col.startswith('SAMPLE')]

        # Apply the lambda function to count '1/1', '1/0', '0/1', '0/0' and sum the counts
        variants['g 1/1'] = df[sample_columns].applymap(lambda x: '1/1' in x.split(':')[0]).sum().sum()
        variants['g 1/0'] = df[sample_columns].applymap(lambda x: '1/0' in x.split(':')[0]).sum().sum()
        variants['g 0/1'] = df[sample_columns].applymap(lambda x: '0/1' in x.split(':')[0]).sum().sum()
        variants['g 0/0'] = df[sample_columns].applymap(lambda x: '0/0' in x.split(':')[0]).sum().sum()

    if include_info:
        for new_col in include_info.split(','):
            df[new_col.strip()] = df.INFO.str.extract(r'{}=([\d.]+)'.format(new_col.strip()), expand=False)

    if verbosity:
        print('Final filtered VCF file:')
        print(df)
        print('\n')
        print('Zygosity breakdown: ')
        for item in variants.items():
            print(str(item[0]) + ': ' + str(item[1]))

    if output_file:
        writing_path = os.path.dirname(vcf_file_path)
        file_name = os.path.basename(vcf_file_path).split('.vcf')[0]

        # convert to dataframe
        zygosity_df = pd.DataFrame(variants, index=[0]).T
        # Reset the index
        zygosity_df.reset_index(inplace=True)
        # Rename the columns
        zygosity_df.columns = ['Genotype', 'Count']

        if output_format == 'csv':
            df.to_csv(writing_path + '/' + file_name + '.filtered.csv', index=False)
            zygosity_df.to_csv(writing_path + '/' + file_name + '.zygosity.csv', index=False)
        elif output_format == 'json':
            df.to_json(writing_path + '/' + file_name + '.filtered.json', orient='records')
            zygosity_df.to_json(writing_path + '/' + file_name + '.zygosity.json', orient='records')
        else:
            print("WARNING: The output format is not specified!")

    print('Process has been successfully done!')


def process_vcf_file(vcf_file_path, sample, depth, af, chromosome, position_range, include_filter, exclude_filter,
                     include_info, output_format, output_file, verbosity):

    """
    Processes a VCF file, filters it based on various parameters, calculates zygosity, and outputs the results.

    @param vcf_file_path The path to the VCF file to process.
    @param sample The name of the sample for which zygosity will be calculated.
    @param depth The depth (DP) threshold for filtering.
    @param af The allele frequency (AF) threshold for filtering.
    @param chromosome The selected chromosome for filtering.
    @param position_range The selected position range for filtering.
    @param include_filter The selected value to filter the FILTER column by, for inclusion.
    @param exclude_filter The selected value to filter the FILTER column by, for exclusion.
    @param include_info A comma-separated string of additional INFO fields to include in the output.
    @param output_format The format of the output file ('csv' or 'json').
    @param output_file Boolean indicating whether to output a file.
    @param verbosity A boolean flag that when True, prints out diagnostic messages to console.
    """

    # read vcf file
    df_vcf = pd.read_csv(vcf_file_path, sep='\t')
    # remove comments (start with ##)
    df_vcf = df_vcf[~df_vcf.iloc[:, 0].str.contains('##')]
    df_vcf.columns = df_vcf.iloc[0]
    df_vcf = df_vcf.iloc[1:].reset_index(drop=True)

    df_vcf = filteration.filter_dp(df_vcf, depth, verbosity)

    if af:
        df_vcf = filteration.filter_af(df_vcf, af, verbosity)
    if chromosome:
        df_vcf = filteration.filter_chr(df_vcf, chromosome, verbosity)

    if position_range:
        df_vcf = filteration.filter_position(df_vcf, position_range, verbosity)

    if include_filter:
        df_vcf = filteration.include_filter(df_vcf, include_filter, verbosity)

    if exclude_filter:
        df_vcf = filteration.exclude_filter(df_vcf, exclude_filter, verbosity)

    calc_zygosity(df_vcf, vcf_file_path, include_info, output_format, output_file, sample_name=sample,
                  verbosity=verbosity)


if __name__ == '__main__':
    # Create the parser
    parser = argparse.ArgumentParser(description='Process some genomic data.')

    # Add the arguments
    parser.add_argument('vcf_file', help='Path to the VCF file')

    parser.add_argument('--sample', type=str, help='Specify the samples of interest.'
                                                   'For example, `--sample SAMPLE2` would only count variants in '
                                                   'SAMPLE2')

    parser.add_argument('--depth', type=float, help='Minimum depth threshold (default: 10)')

    parser.add_argument('--allele-frequency', type=float,
                        help='Specify an allele frequency threshold for considering variants. '
                             'For instance, `--allele-frequency 0.05` would only count variants with an allele '
                             'frequency greater than 0.05.')

    parser.add_argument('--chromosome', type=str, help='Filter variants based on a specific chromosome. '
                                                       'For example, `--chromosome 1` would only process '
                                                       'variants on chromosome 1.')

    parser.add_argument('--position-range', type=str, help='Specify a range of positions to focus on. '
                                                           'For instance, `--position-range 1000-5000` would '
                                                           'only consider variants with positions between 1000'
                                                           ' and 5000.')

    parser.add_argument('--include-filter', type=str, help='Include variants with specific filter flags. '
                                                           'For example, `--include-filter PASS` would only count '
                                                           'variants marked as "PASS" in the FILTER column')

    parser.add_argument('--exclude-filter', type=str, help='Exclude variants with specific filter flags. '
                                                           'For instance, `--exclude-filter LowQual` would exclude '
                                                           'variants marked as "LowQual" in the FILTER column.')

    parser.add_argument('--include-info', type=str, help='Include additional INFO fields in the output. '
                                                         'For example, `--include-info DP` would display the "DP" '
                                                         'INFO field along with other variant information. '
                                                         'In addition, Multiple fields can be specified by placing '
                                                         '"," between them, such as DP,AF')

    parser.add_argument('--output-format', type=str, choices=['csv', 'json'], default='csv',
                        help='Specify the output format for the parsed data.')

    parser.add_argument('--output-file', action='store_true', help='Save the parsed data to a specified output file.')

    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose mode to display detailed progress and additional information.')

    # Parse the arguments
    args = parser.parse_args()

    # run processing
    process_vcf_file(vcf_file_path=args.vcf_file, sample=args.sample, depth=args.depth, af=args.allele_frequency,
                     chromosome=args.chromosome, position_range=args.position_range, include_filter=args.include_filter,
                     exclude_filter=args.exclude_filter, include_info=args.include_info,
                     output_format=args.output_format, output_file=args.output_file, verbosity=args.verbose)
