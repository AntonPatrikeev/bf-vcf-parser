import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pkg_resources




def load_cytobands():
    '''
    Process cytoBandIdeo UCSC annotation file with chromosomes Giemsa staining regions
    downloaded from: http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBandIdeo.txt.gz 
    the file is required for getting chromosomes lengths
    and for plotting chromosomes in a nice way

    Reurns DataFrame with columns: 
    chrom:float, start:float, end:float, name:str, gieStain:str, band_len:float, band_color:tuple
    '''
    
    # Giemsa staining regions to color mapping    
    band_palette = {
        'gneg': (1., 1., 1.),
        'gpos25': (.6, .6, .6),
        'gpos50': (.4, .4, .4),
        'gpos75': (.2, .2, .2),
        'gpos100': (0., 0., 0.),
        'acen': (.8, .4, .4),
        'gvar': (.8, .8, .8),
        'stalk': (.9, .9, .9),
    }

    # load the original cytoBandIdeo file
    stream = pkg_resources.resource_stream(__name__, 'data/cytoBandIdeo.txt.gz')
    cytobands = pd.read_csv(stream, 
                            sep='\t', compression="gzip", header=None, 
                            names=['chrom', 'start', 'end', 'name', 'gieStain'])
    
    cytobands['end'] = cytobands['end'] / 1e6
    cytobands['start'] = cytobands['start'] / 1e6
    cytobands["band_len"] = cytobands['end'] - cytobands['start']
    cytobands["band_color"] = cytobands['gieStain'].map(band_palette)
    
    return cytobands


def plot_variants_distribution(variant_counts_per_mb, substitution_type_per_mb, contig, out_path):
    '''
    Creates variant distribution plot for Variants summary section of the report.
    Generates .png figure in the temporary directory
    '''
    # load cytobands
    cytobands = load_cytobands()
    cytobands = cytobands[cytobands.chrom == contig]

    # get contig length
    contig_len = cytobands['end'].max()

    # generate variants distribution plot 
    fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True, 
                              figsize=(15, 6), gridspec_kw={'height_ratios': [6, 4, 0.5]})

    plot_variants_per_mb(ax[0], variant_counts_per_mb, contig_len)
    plot_substitution_type_per_mb(ax[1], substitution_type_per_mb, contig_len)
    plot_cytobanded_chromosome(ax[2], cytobands, contig)

    fig.tight_layout()
    fig.savefig(out_path, dpi=600)


def plot_variants_per_mb(ax, variant_counts, contig_len):
    '''
    Creates a stacked bar plot with variant types counts per Mb as a part of variant distribution plot 
    for Variants summary section of the report.
    Return axes object.
    '''
    variant_types = ["SNP", "DEL", "INS"]
    
    for idx, variant_type in enumerate(variant_types):
        if idx == 0:
            ax.bar(x = variant_counts["window"], label=variant_type, 
                   height = variant_counts[variant_type], 
                   width = 0.8, 
                   bottom = None)
            cumulative_variant_type_values = np.array(variant_counts[variant_type].tolist())
        else:
            ax.bar(x = variant_counts["window"], label=variant_type,
                   height = variant_counts[variant_type],
                   width = 0.8, 
                   bottom = cumulative_variant_type_values)
            cumulative_variant_type_values += np.array(variant_counts[variant_type].tolist())

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    left, right = 0, contig_len
    ax.set_xlim(left=left, right=right + 0.1)
    stepsize = 5
    ax.xaxis.set_ticks(np.arange(left, right, stepsize))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel("# variants per 1 Mb")
    ax.set_xticks([])
    
    ax.legend(loc='upper center', ncol=3, frameon=False)

    return ax


def plot_substitution_type_per_mb(ax, snp_type_per_mb, contig_len):
    '''
    Creates a stacked bar plot with substitution types % per Mb as a part of variant distribution plot 
    for Variants summary and Samples variant profiles summary sections of the report.
    Return axes object.
    '''
    substitution_types = snp_type_per_mb.columns[1:]
    
    for idx, substitution_type in enumerate(substitution_types):
        if idx == 0:
            ax.bar(x = snp_type_per_mb["window"], label=substitution_type, 
                   height = snp_type_per_mb[substitution_type], 
                   width = 0.8, 
                   bottom = None)
            cumulative_snp_type_per_mb_values = np.array(snp_type_per_mb[substitution_type].tolist())
        else:
            ax.bar(x = snp_type_per_mb["window"], label=substitution_type,
                   height = snp_type_per_mb[substitution_type],
                   width = 0.8, 
                   bottom = cumulative_snp_type_per_mb_values)
            cumulative_snp_type_per_mb_values += np.array(snp_type_per_mb[substitution_type].tolist())

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    ax.set_ylim(bottom=0, top=124)
    ax.yaxis.set_ticks(np.arange(0, 124, 25))
    
    left, right = 0, contig_len
    ax.set_xlim(left=left, right=right + 0.1)
    stepsize = 5
    ax.xaxis.set_ticks(np.arange(left, right, stepsize))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel("% substitution type\nper 1 Mb")
    ax.set_xticks([])
    
    ax.legend(loc='upper center', ncol=6, frameon=False)

    return ax


def plot_cytobanded_chromosome(ax, cytobands, contig):
    '''
    Creates a chromosome visualization based on cytoBandIdeo UCSC annotation as a part of variant distribution plot 
    for Variants summary and Samples variant profiles summary sections of the report.
    Return axes object
    ''' 
    ax.barh(y=[contig], 
            width=cytobands.band_len, 
            left=cytobands.start, 
            color=cytobands.band_color,
            edgecolor="black")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    left, right = 0, cytobands['end'].max()
    ax.set_xlim(left=left, right=right + 0.1)
    stepsize = 5
    ax.xaxis.set_ticks(np.arange(left, right, stepsize))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel(f"position on {contig}, Mb")
    ax.set_yticks([])
    
    return ax


def plot_sample_variants_profiles(sample_variant_type_counts_per_mb, sample_substitution_type_prc_per_mb, contig, out_path):
    '''
    Creates sample variant profiles plot for Samples variant profiles summary section of the report.
    Generates .png figure in the temporary directory
    '''
    # load cytobands
    cytobands = load_cytobands()
    cytobands = cytobands[cytobands.chrom == contig]

    # get contig length
    contig_len = cytobands['end'].max()

    # generate sample genetic variants profiles 
    fig, ax = plt.subplots(nrows=5, ncols=1, sharex=True, 
                              figsize=(15, 8), gridspec_kw={'height_ratios': [3, 3, 3, 3, 0.5]})

    variant_types = ["SNP", "DEL", "INS"]
    for idx, variant_type in enumerate(variant_types):
        # show legend only for the first plot
        show_legend = idx == 0
        # plot HET / HOM stacked bar plot for each Mb
        plot_genotypes_per_mb(ax[idx], 
            sample_variant_type_counts_per_mb[sample_variant_type_counts_per_mb.VT == variant_type], 
            variant_type, contig_len, legend=show_legend)

    plot_substitution_type_per_mb(ax[3], sample_substitution_type_prc_per_mb, contig_len)
    plot_cytobanded_chromosome(ax[4], cytobands, contig)

    fig.tight_layout()
    fig.savefig(out_path, dpi=600)


def plot_genotypes_per_mb(ax, sample_variant_type_counts, variant_type, contig_len, legend):
    '''
    Creates a stacked bar plot with genotype (HOM = homozygous, HET = heterozygous) counts per Mb 
    as a part of variant distribution plot for Samples variant profiles summary section of the report.
    Return axes object.
    '''
    genotypes = ["HOM", "HET"]
    for idx, genotype in enumerate(genotypes):
        if idx == 0:
            ax.bar(x = sample_variant_type_counts["window"], label=genotype, 
                   height = sample_variant_type_counts[genotype], 
                   width = 0.8, 
                   bottom = None)
            cumulative_sample_variant_type_values = np.array(sample_variant_type_counts[genotype].tolist())
        else:
            ax.bar(x = sample_variant_type_counts["window"], label=genotype,
                   height = sample_variant_type_counts[genotype],
                   width = 0.8, 
                   bottom = cumulative_sample_variant_type_values)
            cumulative_sample_variant_type_values += np.array(sample_variant_type_counts[genotype].tolist())

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    left, right = 0, contig_len
    ax.set_xlim(left=left, right=right + 0.1)
    stepsize = 5
    ax.xaxis.set_ticks(np.arange(left, right, stepsize))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel(f'# {variant_type} per 1 Mb')
    ax.set_xticks([])
    
    if legend:
        ax.legend(loc='upper center', ncol=2, frameon=False)

    return ax
