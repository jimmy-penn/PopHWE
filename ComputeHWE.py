import sys
import argparse
import time

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

    args = parser.parse_args()
    st_time = time.time()

    min_pop_size = args.min_pop_size

    geno_cnts_obj = GenotypeCounts(args.geno_cnt_file)
    num_pops = geno_cnts_obj.num_anc_groups
    num_pops = 5

    geno_cnts_obj.ShowSummary()

    # Open the genotype count file and process each line to compute HWE statistics
    # The expected format of each line is: var_id, geno_cnt_pop1, geno_cnt_pop2, ..., geno_cnt_popN
    # where geno_cnt_popX is in the format "aa/ab/bb" representing the counts of each genotype in that population.
    # Cells with invalid format will be treated as genotype counts = 0 (aa + ab + bb = 0) for that population.
    # The script uses similar variable names as in the paper, such as Naa, Saa, Eaa, etc. 
    with open(args.geno_cnt_file, 'r') as file:
        line_no = 0
        for line in file:
            # Skip the header line
            if line_no == 0:
                line_no += 1
                continue
            
            line = line.strip()

            line = "a0	0/4/9	4/0/6	0/0/3	5/0/0	0/5/0"

            vals = line.split(geno_cnts_obj.delimiter)
            num_vals = len(vals)
            print(f"vals: {vals}")
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

                if num_genos == 3 and all(cnt.isdigit() for cnt in geno_cnts):
                    Niaa, Niab, Nibb = map(int, geno_cnts)
                    Ni = sum([Niaa, Niab, Nibb])

                    Naa += Niaa
                    Nab += Niab
                    Nbb += Nibb

                    p = (Niaa + 0.5 * Niab) / Ni
                    q = 1 - p

                    Eiaa = Ni * p**2
                    Eiab = Ni * 2 * p * q
                    Eibb = Ni * q**2

                    Saa += Eiaa
                    Sab += Eiab
                    Sbb += Eibb       

            if Naa + Nab + Nbb >= min_pop_size:
                chi2, Eaa, Eab, Ebb = CalcHwe(Naa, Nab, Nbb)
                eChi, Esaa, Esab, Esbb = CalcHwe(Saa, Sab, Sbb)
                pChi = (Naa-Saa)**2/Saa + (Nab-Sab)**2/Sab + (Nbb-Sbb)**2/Sbb

                print(f"{var_id}:") 
                print(f"\tN: {Naa}\t{Nab}\t{Nbb}") 
                print(f"\tS: {Saa}\t{Sab}\t{Sbb}") 
                print(f"\tE: {Eaa}\t{Eab}\t{Ebb}") 
                print(f"\tChi2: {chi2}\teChi: {eChi}\tpChi: {pChi}")

            line_no += 1
            break

if __name__ == "__main__":
    main()
    
