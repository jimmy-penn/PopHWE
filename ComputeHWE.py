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

    num_pops = geno_cnts.num_anc_groups

    geno_cnts.ShowSummary()

    with open(args.geno_cnt_file, 'r') as file:
        line_no = 0
        for line in file:
            if line_no == 0:
                line_no += 1
                continue
            
            line = line.strip()

            vals = line.split(geno_cnts.delimiter)
            num_vals = len(vals)

            if num_vals != num_pops + 1:
                print(f"ERROR: Line {line_no} has {num_vals} values, expected {num_pops + 1}.")
                exit()

            var_id = vals[0].strip()
            for col in range(num_pops):
                val = vals[col+1].strip()
                cnt_str = str(val)
                aa, ab, bb = val.split('/')
                print(f"SNP: {var_id}, Population: {geno_cnts.anc_groups[col]}, Genotype Counts: {aa}+{ab}+{bb}")
                
            line_no += 1

if __name__ == "__main__":
    main()
    
