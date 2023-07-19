

def filter_af(df, af_threshold, verbosity):
    """
    Filters a dataframe (representing a VCF file) based on a specific allele frequency threshold.

    @param df Pandas DataFrame representing a VCF file.
    @param af_threshold The allele frequency (AF) threshold to filter by. This should be a float between 0 and 1.
    @param verbosity A boolean flag that when True, prints out diagnostic messages to console.

    return: A dataframe filtered by the given allele frequency threshold.
    """

    if 0 <= af_threshold <= 1:
        # Filter rows where the "DP" column is greater than 10
        df_filtered_by_af = df[df.INFO.str.extract(r'AF=([\d.]+)', expand=False).astype(float) > af_threshold]

        if df_filtered_by_af.shape[0] == 0:
            print('WARNING: The VCF file does not contain any rows that accept allele frequency threshold')
        elif verbosity:
            print('allele frequency (AF) = ' + str(af_threshold))
            print('The output of filtration according to the allele frequency (AF) is shown here: ')
            print(df_filtered_by_af)
            print('\n\n')
    else:
        print("WARNING: The allele frequency (AF) is not between 0 and 1. Therefore, this filtration does not apply!")
        df_filtered_by_af = df

    return df_filtered_by_af


def filter_dp(df, dp_threshold, verbosity):
    """
    Filters a dataframe (representing a VCF file) based on a specific depth (DP) threshold.

    @param df Pandas DataFrame representing a VCF file.
    @param dp_threshold The depth (DP) threshold to filter by (default: 10).
    @param verbosity A boolean flag that when True, prints out diagnostic messages to console.

    return: A dataframe filtered by the given depth threshold.
    """

    if not dp_threshold:
        print('WARNING: The depth threshold was not specified. '
              'The filtering process will proceed with a default depth value of 10.')
        dp_threshold = 10
    elif dp_threshold < 0:
        print("WARNING: The depth threshold is negative. "
              "Therefore, The filtering process will proceed with a default depth value of 10.")
        dp_threshold = 10
    # Filter rows where the "DP" column is greater than threshold
    df_filtered_by_dp = df[df.INFO.str.extract(r'DP=(\d+)', expand=False).astype(float) > dp_threshold]

    if df_filtered_by_dp.shape[0] == 0:
        print('WARNING: The VCF file does not contain any rows that accept depth threshold')
    elif verbosity:
        print('depth threshold = ' + str(dp_threshold))
        print('The output of filtration according to the depth threshold is shown here: ')
        print(df_filtered_by_dp)
        print('\n\n')

    return df_filtered_by_dp


def filter_chr(df, selected_chr, verbosity):
    """
    Filters a dataframe (representing a VCF file) based on a selected chromosome.

    @param df  Pandas DataFrame representing a VCF file.
    @param selected_chr The selected chromosome to filter by.
    @param verbosity A boolean flag that when True, prints out diagnostic messages to console.

    return: A dataframe filtered by the selected chromosome.
    """

    df_filtered_by_chr = df[df['#CHROM'].str.extract(r'chr(\d+)', expand=False).astype(int) == selected_chr]

    if df_filtered_by_chr.shape[0] == 0:
        print('WARNING: No rows exist in the VCF file for selected chromosome')
    elif verbosity:
        if verbosity:
            print('chromosome = ' + str(selected_chr))
            print('The output of filtration according to the chromosome is shown here: ')
            print(df_filtered_by_chr)
            print('\n\n')

    return df_filtered_by_chr


def filter_position(df, selected_pos_range, verbosity):
    """
    Filters a dataframe (representing a VCF file) based on a selected position range.

    @param df Pandas DataFrame representing a VCF file.
    @param selected_pos_range The selected position range to filter by. This should be a string in the format "start-end".
    @param verbosity A boolean flag that when True, prints out diagnostic messages to console.

    return: A dataframe filtered by the selected position range.
    """

    try:
        selected_pos_range = [int(val.strip()) for val in selected_pos_range.split('-')]

        if len(selected_pos_range) == 2:
            df_filtered_by_pos = df[
                (df.POS.astype(int) > selected_pos_range[0]) & (df.POS.astype(int) < selected_pos_range[1])]

            if df_filtered_by_pos.shape[0] == 0:
                print('WARNING: No rows exist in the VCF file for selected position range')
            elif verbosity:
                print('position range = ' + '-'.join([str(v) for v in selected_pos_range]))
                print('The output of filtration according to the position range is shown here: ')
                print(df_filtered_by_pos)
                print('\n\n')
        else:
            print('WARNING: the position range is wrong! Therefore, this filtration does not apply!')
            df_filtered_by_pos = df
    except ValueError:
        print("Error: The position range should consist of integers only. Therefore, this filtration does not apply!")
        df_filtered_by_pos = df

    return df_filtered_by_pos


def include_filter(df, selected_val, verbosity):
    """
    Filters a dataframe (representing a VCF file) to include rows where the FILTER column matches a selected value.

    @param df Pandas DataFrame representing a VCF file.
    @param selected_val The selected value to filter the FILTER column by.
    @param verbosity A boolean flag that when True, prints out diagnostic messages to console.

    return: A dataframe filtered to include rows where the FILTER column matches the selected value.
    """

    if selected_val not in df.FILTER.values:
        print('WARNING: This value does not exist in the FILTER columns for including it!')

    df_filtered_by_filter_col = df[df.FILTER == selected_val]

    if verbosity:
        print('included filter = ' + selected_val)
        print('The output of filtration according to the included filter is shown here: ')
        print(df_filtered_by_filter_col)
        print('\n\n')

    return df_filtered_by_filter_col


def exclude_filter(df, selected_val, verbosity):
    """
    Filters a dataframe (representing a VCF file) to exclude rows where the FILTER column matches a selected value.

    @param df Pandas DataFrame representing a VCF file.
    @param selected_val The selected value to filter the FILTER column by.
    @param verbosity A boolean flag that when True, prints out diagnostic messages to console.

    return: A dataframe filtered to exclude rows where the FILTER column matches the selected value.
    """

    if selected_val not in df.FILTER.values:
        print('WARNING: This value does not exist in the FILTER columns for excluding it!')
    df_filtered_by_filter_col = df[df.FILTER != selected_val]

    if verbosity:
        print('Excluded filter = ' + selected_val)
        print('The output of filtration according to the excluded filter is shown here: ')
        print(df_filtered_by_filter_col)
        print('\n\n')

    return df_filtered_by_filter_col
