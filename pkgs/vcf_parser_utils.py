import pandas as pd




def convert_variants_to_df_by_mb(variants_contig_dict):
    ''' 
    Function is used for processing of variants from parsed VCF 
    for further creation of plots in the Variant summary section
    Returns DataFrame with the following columns:
    POS:float, REF:str, ALT:str, VT:str, window:float
    '''
    variants_contig_df = pd.DataFrame.from_dict(variants_contig_dict, 
        orient='index', columns=["REF", "ALT", "VT"])
    variants_contig_df["POS"] = variants_contig_df.index
    variants_contig_df.reset_index(drop=True, inplace=True)

    # assign Mb window to each variant
    variants_contig_df["window"] = variants_contig_df["POS"] // 1000000 + 0.5
    
    return variants_contig_df


def convert_sample_variants_to_df_by_mb(samples_contig_dict, variants_contig_dict):
    ''' 
    Function is used for processing of variants from parsed VCF 
    for further creation of plots in the Samples variant profiles summary section
    Returns DataFrame with the following columns:
    GT:str, POS:float, REF:str, ALT:str, VT:str, window:float
    '''
    samples_contig_df = pd.DataFrame.from_dict(samples_contig_dict, 
        orient='index', columns=["GT"])
    samples_contig_df["POS"] = samples_contig_df.index
    samples_contig_df.reset_index(drop=True, inplace=True)

    # map variants attributes by POS
    for idx, var_attr in enumerate(["REF", "ALT", "VT"]):
        samples_contig_df[var_attr] = samples_contig_df.POS.apply(lambda x: variants_contig_dict[x][idx])

    # assign Mb window to each variant
    samples_contig_df["window"] = samples_contig_df["POS"] // 1000000 + 0.5
    
    return samples_contig_df