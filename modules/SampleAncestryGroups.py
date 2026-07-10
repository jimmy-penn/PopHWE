import os
import argparse

from modules.VcfSamples import VcfSamples

class SampleAncestryGroups:
    """
    Assign VCF samples to ancestry groups derived from a GrafAnc output file

    2026-07-09: Revised to also take a population file with pre-assigned population groups for samples. 
    The population file can be any tab-delimited text file with the following columns: Sample, Population
    """

    def __init__(self, vcf_smps, args):
        self.smp_pop_file = args.grafanc_file  # Name of the Sample Population file
        self.is_grafanc = True                 # Flag indicating if the sample pop file is from GrafAnc
        self.use_sub_anc = args.use_sub_anc    # Whether to assign subjects to subcontiental groups or not
        self.use_groups = args.use_groups      # Whether the current ancestry grouping is acceptable or not
        self.merge_file = args.merge_list      # File to contain ancestry group merge list
        
        self.num_mix_groups = 1        # Number of subgroups to divide admixed groups into
        self.pct_anc_bin = 100         # Bin to sub groups based on pe, pf, or pa values

        self.vcf_smp_obj = vcf_smps    # VcfSamples object
        self.anc_groups = []           # Sorted ancestry groups derived from GrafAnc file
        self.vcf_samples = []          # List of samples in VcfSamples object
        self.grafanc_samples = []      # List of samples in the GrafAnc result file
        self.smp_anc_ids = []          # List of ancestry groups (indices) assigned to every sample
        self.num_vcf_samples = 0       # Number of samples in the VcfSamples object
        self.num_ga_samples = 0        # Number of samples in the GrafAnc output file
        self.num_rm_samples = 0        # Number of samples excluded for HWE test

        self.no_anc_samples = []       # Vcf samples not found in GrafAnc output file
        self.no_vcf_samples = []       # Samples in GrafAnc file but not in the VCF
        self.num_no_anc_samples = 0    # Number of Vcf samples not found in GrafAnc output file
        self.num_no_vcf_samples = 0    # Number of samples in GrafAnc file but not in the VCF
        self.num_overlap_samples = 0   # Number of samples in both GrafAnc file and the VCF

        self.ga_score_names = []       # Names of all GrafAnc scores 
        self.smp_ga_scores = {}        # GrafAnc scores of all samples in the GrafAnc file 

        self.all_anc_groups = []       # All ancestry groups 
        self.anc_group_sizes = {}      # Number of samples assigned to each ancestry group 
        self.num_anc_groups = 0        # Number of ancestry groups 

        self.anc_id_col = -1           # Index of ancestry group ID column in GrafAnc output file
        self.raw_pe_col = -1           # Index of RawPe column in GrafAnc output file
        self.raw_pf_col = -1           # Index of RawPf column in GrafAnc output file
        self.raw_pa_col = -1           # Index of RawPa column in GrafAnc output file

        # Need to find ancestry groups quickly for samples
        # A dict from sample ID to ancestry group ID would be convenient
        # An array can be used to find ancestry group IDs more quickly
        # The following array has same size as vcf_samples, each pos coresponding to a sample
        # The value of each element is the pos in all_anc_groups, showing the anc group of the sample
        self.smp_anc_group_pos = []    # Array showing ancestry group index of each sample

        self.vcf_samples = self.vcf_smp_obj.samples
        self.num_vcf_samples = self.vcf_smp_obj.num_samples

        num_mix_groups = int(args.num_mix_groups)
        if num_mix_groups > 1 and num_mix_groups <= 100:
            self.num_mix_groups = num_mix_groups
            self.pct_anc_bin = int(100.0/self.num_mix_groups)

        self.CheckSamplePopFile()

        self.ReadSampleGrafAncScores()
        self.GetSampleAncestryGroups()
        self.GetAllAncestryGroups()

        self.ShowSummary()

        if self.merge_file:
            self.MergeAncestryGroups()

            print("\nAfter merging ancestry groups, the new ancestry groups are:")
            self.ShowSummary()


    def GetIndexGivenGaColumnName(self, chk_col):
        """
        Given a GrafAnc column name, eg., EU1, return the column index
        """

        chk_col_no = -1
        col_no = 0
        for col_name in self.ga_score_names:
            if col_name == chk_col:
                chk_col_no = col_no

            col_no += 1

        return chk_col_no


    def GetSampleAncestryGroups(self):
        """
        Read or derive sample ancestry groups from GrafAnc output file

        If it is not a GrafAnc output file, it assigns ancestry groups read from the population column.  
        """

        self.smp_anc_ids = []

        if self.is_grafanc:
            self.anc_id_col = self.GetIndexGivenGaColumnName("AncGroupID")
            self.raw_pe_col = self.GetIndexGivenGaColumnName("RawPe")
            self.raw_pf_col = self.GetIndexGivenGaColumnName("RawPf")
            self.raw_pa_col = self.GetIndexGivenGaColumnName("RawPa")

            if self.anc_id_col < 0:
                print(f"\nERROR: didn't find AncGroupID column in {self.smp_pop_file}\n")
                exit()

        ga_anc_id = 0
        for smp_no in range(self.num_vcf_samples):
            smp_id = self.vcf_samples[smp_no]

            if smp_id in self.smp_ga_scores:
                scores = self.smp_ga_scores.get(smp_id)

                if scores:
                    if self.is_grafanc:
                        ga_anc_id = self.AssignAncestryGroupsUsingScores(scores)
                    else:
                        # If it is not a GrafAnc output file, directly assign ancestry group read from the population column
                        ga_anc_id = scores[self.anc_id_col]

                    self.smp_anc_ids.append(ga_anc_id)
                else:
                    self.no_anc_samples.append(smp_id)
            else:
                self.no_anc_samples.append(smp_id)                      


    def MergeAncestryGroups(self):
        """
        Merge ancestry groups based on the user-specified merge list file
        """

        merge_anc_grps = {} # A dict to store the mapping from old ancestry group IDs to new merged ancestry group IDs

        with open(self.merge_file, 'r') as mfh:
            try:
                header = mfh.readline().strip()
            except ValueError as error:
                print(f'\nERROR: failed to read header from {self.merge_file}: {error}\n')
                exit()
            except Exception as error:
                print(f'\nError: {error}\n')
                exit()
                
            cols = header.split("\t")
            exp_columns = ['Ancestries', 'MergeTo']

            if cols != exp_columns:
                print(f"\nERROR: Invalid header in merge list file {self.merge_file}: {header}\n")
                print(f"Expected header is: {'\t'.join(exp_columns)}\n")
                exit()  

            line_no = 0
            for line in mfh:
                line = line.strip()
                if not line:
                    continue

                cols = line.split("\t")
                if len(cols) != 2:
                    print(f"\nERROR: Invalid line in merge list file {self.merge_file}: {line}\n")
                    exit()  

                if self.is_grafanc:
                    if not cols[1].strip().isdigit():
                        print(f"\nERROR: Invalid ancestry group ID in merge list file {self.merge_file}: {cols[1]}\n")
                        exit() 

                is_rm_anc = False
                main_anc_id  = cols[1].strip()
                if self.is_grafanc:
                    main_anc_id  = int(cols[1])
                    if main_anc_id == 0:
                        is_rm_anc = True
                else:
                    if main_anc_id.lower() == "none" or main_anc_id.lower() == "remove" or main_anc_id.lower() == "0":
                        is_rm_anc = True
                    
                if is_rm_anc:
                    if self.is_grafanc:
                        main_anc_id = 0
                    else:
                        main_anc_id = '0'
                else:
                    if not main_anc_id in self.all_anc_groups:
                        print(f"\nERROR: Ancestry group {main_anc_id} in merge list file {self.merge_file} is not found in the current ancestry groups\n")
                        exit()     

                mrg_ancs = []

                ancestries = cols[0].strip('').replace('\"', '').split(",")
                
                for mrg_anc_str in ancestries:
                    if mrg_anc_str.strip() == "":
                        continue

                    if self.is_grafanc:
                        if not mrg_anc_str.strip().isdigit():
                            print(f"\nERROR: Invalid ancestry group ID in merge list file {self.merge_file}: {mrg_anc_str}\n")
                            exit() 

                    mrg_anc = mrg_anc_str.strip()
                    if self.is_grafanc:
                        mrg_anc = int(mrg_anc_str.strip())
    
                    if not mrg_anc in self.all_anc_groups:
                        print(f"\nERROR: Ancestry group {mrg_anc} in merge list file {self.merge_file} is not found in the current ancestry groups\n")
                        exit()  

                    if mrg_anc in merge_anc_grps:
                        print(f"\nERROR: Ancestry group {mrg_anc} is already assigned to ancestry group {merge_anc_grps[mrg_anc]}\n")
                        exit()

                    if mrg_anc != main_anc_id:
                        merge_anc_grps[mrg_anc] = main_anc_id

                line_no += 1
        mfh.close()                               

        for smp_no in range(self.num_vcf_samples):
            anc_id = self.smp_anc_ids[smp_no]

            if anc_id in merge_anc_grps:
                new_anc_id = merge_anc_grps[anc_id]
                self.smp_anc_ids[smp_no] = new_anc_id 

        self.GetAllAncestryGroups()


    def AssignAncestryGroupsUsingScores(self, scores):
        """
        Assign a sample into an ancestry group using its GrafAnc scores and user specifications
        """

        smp_anc_id = 1

        sub_anc_id = int(scores[self.anc_id_col])
        sup_anc_id = int(sub_anc_id/100) * 100

        # For super populations, we still need to separte AMR into AFA, LAT1 and LAT2
        if not self.use_sub_anc:
            if sub_anc_id in (107, 601, 602, 603):
                sup_anc_id = sub_anc_id

        smp_anc_id = sup_anc_id
        if self.use_sub_anc:
            smp_anc_id = sub_anc_id

        # For 107 = AFA, 601 = LAT 1, 602 = LAT 2, 603 = Native Amr, split them into more groups
        if self.num_mix_groups > 1: 
            if sub_anc_id in (107, 601, 602, 603): 
                raw_pf = float(scores[self.raw_pf_col])
                raw_pa = float(scores[self.raw_pa_col])

                bin_id = 0
                if smp_anc_id == 107 or smp_anc_id == 601:
                    bin_id = int(raw_pf / self.pct_anc_bin)
                elif smp_anc_id == 603 or smp_anc_id == 603:
                    bin_id = int(raw_pa / self.pct_anc_bin)
        
                if bin_id < 0:
                    bin_id = 0
                if bin_id > 99:
                    bin_id = 99
        
                bin_id += 1

                smp_anc_id = int(smp_anc_id * 100 + bin_id)
            else:
                smp_anc_id = int(smp_anc_id * 100)

        return smp_anc_id


    def CheckSamplePopFile(self):
        """
        Check the sample population file and decide if it is a GrafAnc output file
        """

        with open(self.smp_pop_file, 'r') as mfh:
            try:
                header = mfh.readline().strip()
            except ValueError as error:
                print(f'\nERROR: failed to read header from {self.smp_pop_file}: {error}\n')
                exit()
            except Exception as error:
                print(f'\nError: {error}\n')
                exit()
                
            columns = header.split("\t")
            num_columns = len(columns)
            if not columns[0].lower().startswith("sample"):
                print(f"\nERROR: Invalid header in sample population file {self.smp_pop_file}: {header}\n")
                print(f"The first column should be a Sample column.\n")
                exit()

            if num_columns > 26 and columns[2] == "GD1" and columns[-1] == "AncGroupID":
                self.is_grafanc = True
            else:
                self.is_grafanc = False

                for col_no in range(1, num_columns):
                    if columns[col_no].lower().startswith("population"):
                        # Ignore the first column, which is Sample ID, so subtract 1 to get the index in the scores array
                        self.anc_id_col = col_no-1  
                        break   

                if self.anc_id_col < 0:
                    print(f"\nERROR: Invalid header in sample population file {self.smp_pop_file}: {header}\n")
                    print(f"Cannot find a Population column.\n")
                    exit()  

        mfh.close()                               


    def ReadSampleGrafAncScores(self):
        """
        Read sample GrafAnc scores and save them into memory.

        If the sample population file is not a GrafAnc output file, it saves some junk into memory, 
        but it will not be used later.
        """

        print(f"Reading sample ancestry groups from {self.smp_pop_file}")

        with open(self.smp_pop_file, 'rt') as afh:
            line_no = 0
            variant_count = 0

            num_columns = 0
            num_values = 0

            for line in afh:
                line = line.strip()
                if not line:
                    continue

                if line_no == 0:
                    columns = line.split('\t')
                    num_columns = len(columns)
                    self.ga_score_names = columns[1:] 
                else:
                    vals = line.split('\t')
                    num_values = len(vals)
                    if num_values != num_columns:
                        print(f"\nERROR: Row {line_no} has {num_values} values, expected {num_columns}")
                        exit()

                    smp_id = vals[0]
                    scores = vals[1:] 
                    self.grafanc_samples.append(smp_id)
                    self.smp_ga_scores[smp_id] = scores
                    self.num_ga_samples += 1

                line_no += 1

        afh.close()

        print(f"Read GrafAnc scores for {self.num_ga_samples} samples")


    def GetAllAncestryGroups(self):
        """
        Summarize ancestry groups assigned to all samples
        """

        smp_anc_grps = {}
        self.anc_group_sizes = {}
        self.all_anc_groups = []

        for i in range(self.num_vcf_samples):
            smp_id = self.smp_anc_ids[i]
            anc_id = self.smp_anc_ids[i]
            smp_anc_grps[smp_id] = anc_id

            is_rm_anc = False
            if (self.is_grafanc and anc_id == 0) or (not self.is_grafanc and anc_id == '0'):
                is_rm_anc = True

            if is_rm_anc:
                self.num_rm_samples += 1
            else:
                if anc_id in self.anc_group_sizes:
                    self.anc_group_sizes[anc_id] += 1
                else:
                    self.anc_group_sizes[anc_id] = 1

        self.all_anc_groups = sorted(self.anc_group_sizes)
        self.num_anc_groups = len(self.all_anc_groups)

        # Get the array from sample indexes to ancestry group indexes
        anc_grp_to_pos = {}
        anc_pos = 0
        for anc in self.all_anc_groups:
            anc_grp_to_pos[anc] = anc_pos
            anc_pos += 1

        self.smp_anc_group_pos = []
        for smp_no in range(self.num_vcf_samples):
            smp_id = self.vcf_samples[smp_no]
            anc_id = self.smp_anc_ids[smp_no]
            anc_pos = -1
            if anc_id in anc_grp_to_pos:
                anc_pos = anc_grp_to_pos[anc_id]
            self.smp_anc_group_pos.append(anc_pos)


    def ShowSummary(self):
        """
        Show summary of the sample ancestry groups
        """

        print(f'\nSummary of sample ancestry groups')
        print(f'\tSamples in VCF:     {self.num_vcf_samples}')
        print(f'\tSamples excluded:   {self.num_rm_samples}')
        print(f'\tFile is GrafAnc:    {self.is_grafanc}')
        print(f'\tSamples in GrafAnc: {self.num_ga_samples}')
        print(f'\tTotal {self.num_anc_groups} ancestry groups:')

        for smp_no in range(self.num_vcf_samples):
            smp_id = self.vcf_samples[smp_no]
            scores = self.smp_ga_scores[smp_id]
            anc_id = self.smp_anc_ids[smp_no]
            anc_pos = self.smp_anc_group_pos[smp_no]

            if smp_no > 100:
                break

        anc_grp_no = 1
        for anc_id in self.all_anc_groups:
            num_smps = self.anc_group_sizes[anc_id]
            print(f'\t\tNo. {anc_grp_no}\tanc_group_id {anc_id}:  {num_smps}  samples')

            anc_grp_no += 1

        print('')
