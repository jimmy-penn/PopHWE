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

    y = 0
    if (denom > 0):
        y = 1 / denom
        
    return  y

def PlotChi2Hist(args):
    print(f'in {args.input_file} out {args.output_file} width {args.gwidth} height {args.gheight}')
        
    fig, ax = plt.subplots(figsize=(args.gwidth, args.gheight), layout='constrained')
            
    pf = None
    try:
        pf = pd.read_csv(args.input_file, sep="\t")
    except Exception as error:
        print(f'\nERROR: Cannot open {args.input_file}: {error}\n')
        exit()
        
    num_rows = len(pf[args.xcol])
    print(f'Reading {num_rows} rows from {args.input_file}\n')

    exit()

        
def main():
    """
    """
    parser = argparse.ArgumentParser(description="This script plots HWE test results.")
    parser.add_argument("input_file", help="HWE sim result")
    parser.add_argument("output_file", help="The .png file to save the scatterplot")
    parser.add_argument("--gwidth", help="Graph width in inches", metavar="", default=8)
    parser.add_argument("--gheight", help="Graph height in inches", metavar="", default=8)
    parser.add_argument("--xmin", help="Min x value to plot", metavar="")
    parser.add_argument("--xmax", help="Max x value to plot", metavar="")
    parser.add_argument("--ymin", help="Min y value to plot", metavar="")
    parser.add_argument("--ymax", help="Max y value to plot", metavar="")
    parser.add_argument("--xcol", help="X column to plot", metavar="", default="VariantID")
    parser.add_argument("--ycol", help="Y column to plot", metavar="", default="PChi2")

    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        print(f"\nERROR: didn't find file {args.input_file}\n")
        exit()

    hwe_result_file = HweResultFile(args.input_file)
    if not hwe_result_file.is_valid:
        print(f"\nERROR: The input file {args.input_file} is not a valid HWE result file.\n")
        exit() 

    if args.scatter:
        print(f'Plotting scatterplot')
#        PlotChi2Hist(args)
    else:
        print(f'Plotting histogram with bin size {args.bin_size} w {args.gwidth} h {args.gheight} xcol {args.xcol} ycol {args.ycol}')
        PlotChi2Hist(args)
    
if __name__ == "__main__":
    main()
    