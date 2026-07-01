import argparse

import os
import string
import sys
import math
import random

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
import matplotlib.patches as patches

import matplotlib.colors as mcolors
from pathlib import Path
from matplotlib import colormaps
from matplotlib import ticker

import numpy as np
from sklearn.linear_model import LinearRegression

from scipy.stats.distributions import chi2
from scipy.optimize import curve_fit

from modules.HweResultFile import HweResultFile
from modules.GraphParameters import GraphParameters


def func(x, a, b):
    y = 0
    
    denom = (2*math.pi*x*b)**0.5 * math.exp(b*x/2)
    if np.any(denom <= 0):
        y = 0
    else:
        y = a / denom

    return  y
  
def ChiPdf(x):
    y = 0

    if x < 100:
        denom = (2*math.pi*x)**0.5 * math.exp(x/2)
    
        if (denom > 0):
            y = 1 / denom
        
    return  y

def ChiD2Pdf(x):
    denom = 2 * math.exp(x/2)

    y = 1 / denom if (denom > 0) else 0
        
    return  y

def PlotHistogram(params, in_file, out_file):
    """
    Plot histogram for the selected HWE variable
    """
    fig, ax = plt.subplots(figsize=(params.width, params.height), layout='constrained')

    pf = None
    try:
        pf = pd.read_csv(in_file, sep="\t")
    except Exception as error:
        print(f'\nERROR: Cannot open {in_file}: {error}\n')
        exit()
        
    num_rows = len(pf[params.xvar])
    num_snps = num_rows - 1
    print(f'Reading {num_rows} rows from {in_file}\n')

    xvals = []
    for i in range(num_snps):
        x = pf[params.xvar][i+1]
        xvals.append(x)

    print(f'Total {num_snps} SNPs for plotting.\n')

    min_chi = params.xmin
    max_chi = params.xmax

    histo_vals = []        
    num_bins = round((max_chi - min_chi) / params.bin_size) + 1
    
    bin_starts = []
    bin_cnts = []

    for bin_no in range(num_bins):
        bin_st = min_chi + bin_no * params.bin_size
        bin_starts.append(bin_st)
        bin_cnts.append(0)

    bin_starts.append(max_chi + params.bin_size)

    for i in range(num_snps):
        x = xvals[i]
        bin_id = 0
            
        if x > min_chi and x < max_chi:
            bin_id = int((x * 1.0 - min_chi) / params.bin_size)
            if bin_id < num_bins:
                bin_cnts[bin_id] += 1

    max_bin_id = 0
    max_bin_cnt = 0
    tot_cnts = 0
    for bin_no in range(num_bins):
        cnt = bin_cnts[bin_no]
        tot_cnts += cnt
        if cnt > max_bin_cnt:
            max_bin_id = bin_no
            max_bin_cnt = cnt
#        print(f'BinNo. {bin_no} cnt {cnt}')

    if params.show_bin_cnts:            
        print(f'No. {max_bin_id} cnt {max_bin_cnt}')
                
    ax.set_xlim(params.xmin, params.xmax)
    ax.set_ylim(params.ymin, params.ymax)

    bins = np.array(bin_starts) 
    cnts = np.array(bin_cnts) 
    
    ax.hist(bins[:-1], bins, weights=cnts, rwidth=1.0)
    if params.chi_dist:
        x_pdf = []
        y_pdf = []

        x = 0
        y = 0
            
        for bin_no in range(num_bins):
            bin_st = min_chi + bin_no * params.bin_size
            cnt = bin_cnts[bin_no]
            tot_cnts += cnt

            if bin_st > params.xmin and bin_st < params.xmax:
                x = bin_st
                y = ChiPdf(x) * num_snps * params.bin_size
    
                x_pdf.append(x-params.bin_size*0) # Why shifting to right by two bins
                y_pdf.append(y)

        print(f'Total {tot_cnts} SNPs in {num_bins} bins')
        ax.plot(x_pdf, y_pdf, linewidth=params.line_width, color="red")

    ax.set_xlabel(params.disp_xvar, fontsize=params.axis_lbl_font)
    ax.set_ylabel("Number of SNPs", fontsize=params.axis_lbl_font)

    ax.tick_params(axis='x', labelsize=params.tick_lbl_size)

    if params.title != "":
        ax.set_title(params.title, fontsize=params.title_font)

    try:
        plt.savefig(out_file)
    except Exception as error:
        print(f'\nERROR: Cannot save results to {out_file}: {error}\n')
    else:
        print(f'\nSaved results to {out_file}\n')


def PlotScatterplot(params, in_file, out_file):
    """
    Plot scatterplot showing relationship of two selected HWE variables
    """

    fig, ax = plt.subplots(figsize=(params.width, params.height), layout='constrained')

    pf = None
    try:
        pf = pd.read_csv(in_file, sep="\t")
    except Exception as error:
        print(f'\nERROR: Cannot open {in_file}: {error}\n')
        exit()
        
    num_rows = len(pf[params.xvar])
    num_snps = num_rows - 1
    print(f'Reading {num_rows} rows from {in_file}\n')

    # All variants or all except excluded variants when specified
    xvals = []
    yvals = []
    colors = []

    # All heterozygous excess variants when specified
    xhvals = []
    yhvals = []
    hcolors = []

    max_xval, max_yval = 0, 0

    for i in range(num_snps):

        x = pf[params.xvar][i+1]
        y = pf[params.yvar][i+1]

        if x > max_xval:
            max_xval = x
        if y > max_yval:
            max_yval = y

        if params.show_hetexc:
            is_het_exc = pf[params.het_column][i+1]
            if is_het_exc:
                xhvals.append(x)
                yhvals.append(y)
                hcolors.append(params.het_color)
            else:
                xvals.append(x)
                yvals.append(y)
                colors.append(params.dot_color)
        else:
            xvals.append(x)
            yvals.append(y)
            colors.append(params.dot_color)

    print(f'Total {num_snps} SNPs for plotting.')        
    print(f'x: 0 - {max_xval} y: 0 - {max_yval}')        
        
    ax.scatter(xvals, yvals, s=params.dot_size, c=colors, data=pf)
    ax.scatter(xhvals, yhvals, s=params.dot_size, c=hcolors, data=pf)
    ax.plot([params.xmin, params.xmax],[params.ymin, params.ymax], linewidth=params.line_width, color=params.line_color)

    ax.set_xlabel(params.disp_xvar, fontsize=params.axis_lbl_font)
    ax.set_ylabel(params.disp_yvar, fontsize=params.axis_lbl_font)
    ax.set_xlim(params.xmin, params.xmax)
    ax.set_ylim(params.ymin, params.ymax)

    if params.title != "":
        ax.set_title(params.title, fontsize=params.title_font)

    try:
        plt.savefig(out_file)
    except Exception as error:
        print(f'\nERROR: Cannot save results to {out_file}: {error}\n')
    else:
        print(f'\nSaved results to {out_file}\n')


def PlotHweResults(params, in_file, out_file):
    """
    Plot HWE results either in a histogram or a scatterplot, specified by the user
    """

    if params.plot_type == 'histo':
        PlotHistogram(params, in_file, out_file)
    else:
        PlotScatterplot(params, in_file, out_file)


def main():
    """
    Main function to parse arguments and generate the plot.
    """
    parser = argparse.ArgumentParser(description="This script plots HWE test results in either histogram (default) \
        or scatter plot (with --scatter option). PopHWE generates 3 HWE scores: Chi2, EChi2, and PChi2. When plotting\
        a histogram, you can select a score to plot with --xcol, and specify the bin size with --bin_size. \
        When plotting a scatter plot, you can specify which two scores to plot with --xcol and --ycol options. ")
    parser.add_argument("input_file", help="HWE sim result")
    parser.add_argument("output_file", help="The .png file to save the scatterplot")
    parser.add_argument("--gw", help="Graph width in inches", metavar="", default=10)
    parser.add_argument("--gh", help="Graph height in inches", metavar="", default=10)
    parser.add_argument("--xmin", help="Min x value to plot", metavar="")
    parser.add_argument("--xmax", help="Max x value to plot", metavar="")
    parser.add_argument("--ymin", help="Min y value to plot", metavar="")
    parser.add_argument("--ymax", help="Max y value to plot", metavar="")
    parser.add_argument("--xvar", help="Variable to be plotted on x-axis", metavar="", default="PChi2")
    parser.add_argument("--yvar", help="Variable to be plotted on y-axis", metavar="", default="EChi2")
    parser.add_argument("--bin_size", help="Bin size for chi", metavar="")
    parser.add_argument("--chi_dist", help="Show chi2 distribution with df=1", action="store_true")
    parser.add_argument("--dot_size", help="Dot size (0.5, 1, 3, etc) of the scatter plot", metavar="", default=1)
    parser.add_argument("--title", help="The title of the scatterplot", metavar="", default="")
    parser.add_argument("--scatter", help="Show scatterplot", action="store_true")
    parser.add_argument("--show_bin_cnts", help="Show bins and counts", action="store_true")
    parser.add_argument("--show_hetexc", help="Highlight HetExc variants", action="store_true")
    parser.add_argument("--no_outlier", help="Do not show outliers", action="store_true")
    parser.add_argument("--tight_plot", help="Generate a tight plot without margins", action="store_true")
    parser.add_argument("--nolimits", help="Let python select axis limits", action="store_true")

    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        print(f"\nERROR: didn't find file {args.input_file}\n")
        exit()

    hwe_result_file = HweResultFile(args.input_file)
    if not hwe_result_file.is_valid:
        print(f"\nERROR: The input file {args.input_file} is not a valid HWE result file.\n")
        exit() 

    params = GraphParameters(args, hwe_result_file)
    params.ShowSummary()

    print(f'\nPlotting HWE-test results')        
    PlotHweResults(params, args.input_file, args.output_file)

if __name__ == "__main__":
    main()
    