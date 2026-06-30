import argparse

import os
import string
import sys
import math
import random
import time

import pandas as pd
import gzip
import modules.MyUtil as MyUtil

from modules.VcfSamples import VcfSamples
from modules.SampleAncestryGroups import SampleAncestryGroups

def ComputePopGenoCntsFromVcf(smp_anc_grps, args):
    """
    Read genotypes from a VCF file and calculate genotype counts for all ancestry groups.
    """

    vcf_file = smp_anc_grps.vcf_smp_obj.vcf_file

    print(f"Adm groups {args.num_mix_groups}") 
    if args.use_sub_anc:
        print(f"Use sub") 

    exit()

    ofh = open(args.output_file, "w")
    if not ofh:
        print(f"ERROR: cannot create output file {args.output_file} for writing!")
        exit()  

    out_header = "VariantID"
    for anc_id in smp_anc_grps.all_anc_groups:
        out_header += f"\tAncID_{anc_id}"

    print(f"head {out_header}")

    ofh.write(out_header + "\n")

    # VCF file might be very large. Genotype counts are calculated line by line and saved to output file 
    # in blocks of lines. This is to save memory and also to save some results to output file while 
    # calculating genotype counts from the input file. The block size can be adjusted as needed.
    block_size = 10000         
    num_lines_in_block = 0  
    save_lines = [] 
    num_save_lines = 0
    num_samples = 0
    mtime0 = time.time()
    mtime1 = time.time()
    mtime2 = time.time()

    with gzip.open(vcf_file, 'rt') as vfh:
        head_samples = []
        line_count = 0
        num_snps = 0
        
        for line in vfh:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('#'):
                # Parse column headers to find sample names
                if line.startswith('#CHROM'):
                    columns = line.split('\t')
                    samples = columns[9:]  # Samples start at column 10
                    num_samples = len(samples)

                    print(f"Total {num_samples} samples in vcf file\n")
                    if num_samples != smp_anc_grps.num_vcf_samples:
                        print(f"\nERROR: Wrong number of samples found in VCF. Expected {smp_anc_grps.num_vcf_samples}\n")
                        exit()
                continue
            
            # Parse variant lines
            values = line.split('\t')
            chrom, pos, snp_id, ref, alt = values[0], values[1], values[2], values[3], values[4]

            use_var_id = f"{chrom}_{pos}_{ref}_{alt}"
            if snp_id and snp_id.lower().startswith("rs"):
                use_var_id = snp_id
                
            # Save ancestry group genotype counts in 3 arrays                
            g00s = []
            g01s = []
            g11s = []
            for anc_no in range(smp_anc_grps.num_anc_groups):
                g00s.append(0)
                g01s.append(0)
                g11s.append(0)
            
            # Read genotypes for each variant and count for each ancestry group
            num_genos = len(values) - 9

            for smp_no in range(num_samples):
                anc_grp_no = smp_anc_grps.smp_anc_group_pos[smp_no]

                gt = values[9+smp_no]
                a1, a2 = gt[0], gt[2]
 
                if   a1 == "0" and a2 == "0":
                    g00s[anc_grp_no] += 1
                elif a1 == "0" and a2 == "1":
                    g01s[anc_grp_no] += 1
                elif a1 == "1" and a2 == "0":
                    g01s[anc_grp_no] += 1
                elif a1 == "1" and a2 == "1":
                    g11s[anc_grp_no] += 1
                                            
            save_line = f"{use_var_id}"
            for anc_no in range(smp_anc_grps.num_anc_groups):
                save_line += f"\t{g00s[anc_no]}/{g01s[anc_no]}/{g11s[anc_no]}"

            save_lines.append(save_line)                
            num_lines_in_block += 1

            if num_lines_in_block >= block_size:
                ofh.write("\n".join(save_lines) + "\n")
                save_lines = []
                num_lines_in_block = 0
                num_save_lines += block_size

                mtime2 = time.time()
                tot_time = MyUtil.GetTimeDiffString(mtime0, mtime2)
                diff_time = MyUtil.GetTimeDiffString(mtime1, mtime2)
                mtime1 = mtime2
                print(f"\tSaved {num_save_lines} lines to output file {args.output_file}. Time used: {diff_time}. Total {tot_time}")

            num_snps += 1
#            if num_snps > 40000:
#                break

        vfh.close()

        # Write any remaining lines in the block to the output file
        if save_lines:
            ofh.write("\n".join(save_lines) + "\n")
            num_save_lines += len(save_lines)

            mtime2 = time.time()
            tot_time = MyUtil.GetTimeDiffString(mtime0, mtime2)
            diff_time = MyUtil.GetTimeDiffString(mtime1, mtime2)
            mtime1 = mtime2
            print(f"Done. Saved total {num_save_lines} lines to output file {args.output_file}. Time used: {diff_time}. Total {tot_time}")

    ofh.close()


def main():
    """
    Read GrafAnc assigned anc IDs. Check subject counts in each group.
    For admixed populations LT1 and AAM, split them into more groups using RawPf values
    For admixed populations LT2 and NAM, split them into more groups using RawPa values
    """
    parser = argparse.ArgumentParser(description="This script calculates variant allele frequencies by ancestry groups.")
    parser.add_argument("vcf_file", help="Vcf file with genotypes")
    parser.add_argument("grafanc_file", help="File generated by GrafAnc showing sample ancestry groups")
    parser.add_argument("output_file", help="The .png file to save the scatterplot")
    parser.add_argument("--use_sub_anc", help="Use sub-continental ancestry groups", action="store_true")
    parser.add_argument("--num_mix_groups", help="Further divide admixed groups into xx subgroups", metavar="", default=1)

    args = parser.parse_args()

    if not os.path.isfile(args.vcf_file):
        print(f"\nERROR: didn't find file {args.vcf_file}\n")
        exit()

    st_time = time.time()

    vcf_smps = VcfSamples(args.vcf_file)
    if vcf_smps.file_is_valid:
        print(f"\nFound {vcf_smps.num_samples} samples in vcf file")
    else:
        print(f"Found errors in vcf file {args.vcf_file}\n")
        exit()

    smp_anc_grps = SampleAncestryGroups(vcf_smps, args)
    smp_anc_grps.ShowSummary()
    
    ComputePopGenoCntsFromVcf(smp_anc_grps, args)

    ed_time = time.time()
    diff_time = MyUtil.GetTimeDiffString(st_time, ed_time)
    print(f'\nDone. Time used: {diff_time}\n')
    
if __name__ == "__main__":
    main()
    
