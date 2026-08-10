import sys
import argparse
import time
import math

from modules.GenotypeCounts import GenotypeCounts

def CalcHwe(aa, ab, bb):
    """
    This function computes the chi-square statistic for the HWE test, 
    given the counts of the three genotypes (aa, ab, bb).
    """

    chi = 0
    (eaa, eab, ebb) = (0, 0, 0)    
    tot = aa + ab + bb
    
    if tot > 0:
        p = (2.0 * aa + ab) / (2 * tot)
        q = 1 - p
        eaa = tot * p**2
        eab = tot * p * q * 2
        ebb = tot * q**2
        
        if (eaa * eab * ebb > 0):
            chi = (aa - eaa)**2/eaa + (ab - eab)**2/eab + (bb - ebb)**2/ebb;
    
    return (chi, eaa, eab, ebb)


def main():
    parser = argparse.ArgumentParser(description="This script computes HWE-test results from genotype counts.")
    parser.add_argument("geno_cnt_file", help="Name (with full path) of the genotype count file")
    parser.add_argument("output_file", help="Name of the output file")
    parser.add_argument("--min_pop_size", help="Minimum population size for chi-square test", metavar="", type=int, default=10)
    parser.add_argument("--min_MAF", help="Minimum minor allele frequency", metavar="", type=float, default=0.01)

    args = parser.parse_args()
    st_time = time.time()

    min_pop_size = args.min_pop_size
    min_maf = args.min_MAF

    geno_cnt_file = args.geno_cnt_file
    output_file = args.output_file

    geno_cnts_obj = GenotypeCounts(geno_cnt_file)
    num_pops = geno_cnts_obj.num_anc_groups
    
    geno_cnts_obj.ShowSummary()

    out_file = open(output_file, "w")
    if not out_file:
        print(f"ERROR: cannot create output file {output_file} for writing!")
        exit()  

    out_header = "VariantID\tNaa\tNab\tNbb\tSaa\tSab\tSbb\tHetExc\tChi2\tEChi2\tPChi2\tz-value\tDChi2"
    out_file.write(out_header + "\n")

    # Input file might be very large. HWE results are calculated line by line and saved to output file 
    # in blocks of lines. This is to save memory and also to save some results to output file while 
    # calculating HWE from the input file. The block size can be adjusted as needed.
    block_size = 100000         
    num_lines_in_block = 0  
    save_lines = [] 
    num_save_lines = 0

    # Open the genotype count file and process each line to compute HWE statistics
    # The expected format of each line is: var_id, geno_cnt_pop1, geno_cnt_pop2, ..., geno_cnt_popN
    # where geno_cnt_popX is in the format "aa/ab/bb" representing the counts of each genotype in that population.
    # Cells with invalid format will be treated as genotype counts = 0 (aa + ab + bb = 0) for that population.
    # The script uses similar variable names as in the paper, such as Naa, Saa, Eaa, etc. 
    with open(geno_cnt_file, 'r') as cnt_file:
        line_no = 0
        for line in cnt_file:
            # Skip the header line
            if line_no == 0:
                line_no += 1
                continue
            
            line = line.strip()

            vals = line.split(geno_cnts_obj.delimiter)
            num_vals = len(vals)
            if num_vals != num_pops + 1:
                print(f"ERROR: Line {line_no} has {num_vals} values, expected {num_pops + 1}.")
                exit()

            var_id = vals[0].strip()

            Naa, Nab, Nbb = 0, 0, 0  # genotype counts of all populations for aa, ab, bb respectively
            Saa, Sab, Sbb = 0, 0, 0  # sum of expected genotype counts for aa, ab, bb respectively
            for col in range(num_pops):
                val = vals[col+1].strip()
                cnt_str = str(val)
                geno_cnts = val.split('/')
                num_genos = len(geno_cnts)

                Niaa, Niab, Nibb = 0, 0, 0  # genotype counts of population i for aa, ab, bb respectively
                Ni = 0                      # total genotype count of population i       
                p = 1
                q = 0

                if num_genos == 3 and all(cnt.isdigit() for cnt in geno_cnts):
                    Niaa, Niab, Nibb = map(int, geno_cnts)
                    Ni = sum([Niaa, Niab, Nibb])

                    Naa += Niaa
                    Nab += Niab
                    Nbb += Nibb

                    if Ni > 0:
                        p = (Niaa + 0.5 * Niab) / Ni
                        q = 1 - p

                        Eiaa = Ni * p**2
                        Eiab = Ni * 2 * p * q
                        Eibb = Ni * q**2

                        Saa += Eiaa
                        Sab += Eiab
                        Sbb += Eibb       

            num_tot_smps = Naa + Nab + Nbb
            tot_p = (2.0 * Naa + Nab) / (2 * num_tot_smps) if num_tot_smps > 0 else 1
            maf = min(tot_p, 1 - tot_p)

            if num_tot_smps >= min_pop_size and maf > min_maf:
                chi2, Enaa, Enab, Enbb = CalcHwe(Naa, Nab, Nbb)
                echi2, Esaa, Esab, Esbb = CalcHwe(Saa, Sab, Sbb)
                pchi2 = 0
                z_val = 0
                if Saa * Sab * Sbb > 0:
                    pchi2 = (Naa-Saa)**2/Saa + (Nab-Sab)**2/Sab + (Nbb-Sbb)**2/Sbb
                    zval = math.sqrt(pchi2)

                is_het_exc = True 
                if Nab < Enab:
                    is_het_exc = False
                    zval *= -1

                dchi2 = abs(chi2 - echi2)

                save_line = f"{var_id}\t{Naa}\t{Nab}\t{Nbb}\t{Saa:.2f}\t{Sab:.2f}\t{Sbb:.2f}\t{is_het_exc}\t{chi2:.4f}\t{echi2:.4f}\t{pchi2:.4f}\t{zval:.4f}\t{dchi2:.4f}"
                save_lines.append(save_line)
                num_lines_in_block += 1

                if num_lines_in_block >= block_size:
                    out_file.write("\n".join(save_lines) + "\n")
                    save_lines = []
                    num_lines_in_block = 0
                    num_save_lines += block_size
                    print(f"\tSaved {num_save_lines} lines to output file {output_file}")

            line_no += 1
            
        cnt_file.close()

        # Write any remaining lines in the block to the output file
        if save_lines:
            out_file.write("\n".join(save_lines) + "\n")
            num_save_lines += len(save_lines)
            print(f"Done. Saved total {num_save_lines} lines to output file {output_file}")

    out_file.close()


if __name__ == "__main__":
    main()
    
