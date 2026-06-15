import sys
import argparse
import time

from modules.GenotypeCounts import GenotypeCounts

def main():
    parser = argparse.ArgumentParser(description="This script computes HWE-test results from genotype counts.")
    parser.add_argument("geno_cnt_file", help="Name (with full path) of the genotype count file")
    parser.add_argument("output_file", help="Name of the output file")

    args = parser.parse_args()
    st_time = time.time()
    geno_cnts = GenotypeCounts(args.geno_cnt_file)
 
    geno_cnts.ShowSummary()


if __name__ == "__main__":
    main()
    
