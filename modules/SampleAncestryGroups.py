import os

from modules.VcfSamples import VcfSamples

class SampleAncestryGroups:
    """
    Assign VCF samples to ancestry groups derived from a GrafAnc output file
    """

    def __init__(self, vcf_smps, anc_file):
        self.grafanc_file = anc_file   # Name of the GrafAnc output file
        self.vcf_smp_obj = vcf_smps    # VcfSamples object
        self.anc_groups = []           # Sorted ancestry groups derived from GrafAnc file
        self.vcf_samples = []          # List of samples in VcfSamples object
        self.grafanc_samples = []      # List of samples in the GrafAnc result file
        self.smp_anc_ids = []          # List of ancestry groups (indices) assigned to every sample
        self.num_vcf_samples = 0       # Number of samples in the VcfSamples object
        self.num_ga_samples = 0        # Number of samples in the GrafAnc output file

        self.no_anc_samples = []       # Vcf samples not found in GrafAnc output file
        self.no_vcf_samples = []       # Samples in GrafAnc file but not in the VCF
        self.num_no_anc_samples = 0    # Number of Vcf samples not found in GrafAnc output file
        self.num_no_vcf_samples = 0    # Number of samples in GrafAnc file but not in the VCF
        self.num_overlap_samples = 0   # Number of samples in both GrafAnc file and the VCF

        self.ga_score_names = []       # Names of all GrafAnc scores 
        self.smp_ga_scores = {}        # GrafAnc scores of all samples in the GrafAnc file 

        self.vcf_samples = self.vcf_smp_obj.samples
        self.num_vcf_samples = self.vcf_smp_obj.num_samples

        self.ReadSampleGrafAncScores()
        self.GetSampleAncestryGroups()


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
        """

        self.smp_anc_ids = []

        for smp_no in range(self.num_vcf_samples):
            anc_id_col = self.GetIndexGivenGaColumnName("AncGroupID")
            if anc_id_col < 0:
                print(f"\nERROR: didn't find AncGroupID column in {self.grafanc_file}\n")
                exit()

            smp_id = self.vcf_samples[smp_no]
            ga_anc_id = 0
            scores = self.smp_ga_scores.get(smp_id)
            if scores:
                ga_anc_id = scores[anc_id_col]
                print(f"sample {smp_id} scores {scores} col {anc_id_col} ANC ID {ga_anc_id}")
            else:
                self.no_anc_samples.append(smp_id)

            self.smp_anc_ids.append(ga_anc_id)

            if smp_no > 10:
                break

        

    def ReadSampleGrafAncScores(self):
        """
        Read sample GrafAnc scores and save them into memory
        """

        print(f"Reading sample ancestry groups from {self.grafanc_file}")

        with open(self.grafanc_file, 'rt') as afh:
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


    def SetSampleAncestryGroups(self):
        """
        Set ancestry groups based on GrafAnc-assigned groups and GrafAnc scores
        """

        print(f"Setting sample ancestry groups")


    def ShowSummary(self):
        """
        Show summary of the sample ancestry groups
        """

        print(f'\nSummary')
        print(f'\tSamples in VCF:     {self.num_vcf_samples}')
        print(f'\tSamples in GrafAnc: {self.num_ga_samples}')
        print(f'\tGA scores: {', '.join(self.ga_score_names)}.')
        print(f'\tAncestry groups:')

        for smp_no in range(self.num_vcf_samples):
            smp_id = self.vcf_samples[smp_no]
            scores = self.smp_ga_scores[smp_id]
            anc_id = self.smp_anc_ids[smp_no]

            print(f'\t\tsmp {smp_id} anc_id {anc_id}')

            if smp_no > 10:
                break

        print('')
