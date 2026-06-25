import os
import gzip

class VcfSamples:
    """
    Read samples from a vcf.gz file
    """

    def __init__(self, filename):
        self.vcf_file = filename     # Name of the vcf file
        self.file_is_valid = True    # If the vcf file is valid
        self.samples = []            # All samples in the vcf file
        self.num_samples = 0         # Number of samples
        self.min_snps = 10           # Min required variables for a vcf file
        self.min_samples = 10        # Min required samples for a vcf file

        self.ReadSamplesFromVcf()


    def ReadSamplesFromVcf(self):
        """
        Read samples from the vcf file.
        """
        if not os.path.isfile(self.vcf_file):
            print(f"\nERROR: didn't find VCF file {self.vcf_file}\n")
            self.file_is_valid = False

            return

        with gzip.open(self.vcf_file, 'rt') as f:
            # Skip header lines
            num_variants = 0
            num_columns = 0

            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('#'):
                    # Parse column headers to find sample names
                    if line.startswith('#CHROM'):
                        columns = line.split('\t')
                        self.samples = columns[9:] 
                        self.num_samples = len(self.samples)

                        num_columns = len(columns)

                    continue
                
                values = line.split('\t')
                num_vals = len(values)

                num_variants += 1

                if num_vals != num_columns:
                    print(f"\nERROR: header has {num_columns} columns, but line {num_variants} has {num_vals} values\n")
                    self.file_is_valid = False
                    break
                                    
                if num_variants >= self.min_snps:
                    break

        if num_variants < self.min_snps:
            if num_variants == 0:
                print(f"\nERROR: VCF file {self.vcf_file} is empty\n")
            else:
                print(f"\nERROR: only {num_variants} variants are found in VCF file {self.vcf_file}\n")

            self.file_is_valid = False
    
        if self.num_samples < self.min_samples:
            print(f"\nERROR: only {self.num_samples} samples found in VCF file, not enough for HWE test.\n")

            self.file_is_valid = False

