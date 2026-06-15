import os
import pandas as pd

import re
import math

class GenotypeCounts:
    """
    Read population genotype counts and save them into memory

    Requirements: 
    1. Genotype counts are stored in a .txt file containing a table with column header
    2. The first column has header "Variant ID", containing non-redundant variant IDs 
    3. The headers of other columns are names of ancestry groups, each containing 3 integers
       separted by '/', e.g., 12/231/827, representing counts of the 3 genotypes AA, AB, BB
    4. Ancestry groups should be unique
    """

    def __init__(self, filename):
        ### Member variables
        self.filename = filename     # Name of the SSM file
        self.delimiter = '\t'        # Delimiter of columns, tab or comma

        self.num_groups = 0          # Number of ancestry groups
        self.num_variants = 0        # Number of variants
        self.anc_groups = []         # All ancestry groups
        self.num_anc_groups = 0      # Number of ancestry groups
        self.dup_groups = []         # All duplicate ancestry group names
        self.dup_var_ids = []        # All duplicate variant IDs
        self.num_dup_groups = 0      # Number of duplicate ancestry group names
        self.num_dup_var_ids = 0     # Number of duplicate variant IDs

        self.CheckFileFormat()

    def CheckFileFormat(self):
        """
        Check if genotype count file has a correct format
        """
        try:
            file = open(self.filename, "r")
        except Exception as error:
            print(f'\nERROR: Cannot open {self.filename}: {error}\n')
            exit()

        header = file.readline()

        vars = header.split(self.delimiter)
        num_vars = len(vars)

        id_var = ''
        for col in range(num_vars):
            var = str(vars[col])
            var = var.strip()
            if col == 0:
                id_var = var
            else:
                self.anc_groups.append(var)

        file.close()

      
    def ShowSummary(self):
        """
        Show summary of the genotype count file
        """
        print(f'Genotype count file {self.filename}')
        print(f'\tNumber of ancestry groups: {self.num_anc_groups}')
        for anc_grp in self.anc_groups:
            print(f'\t\t{anc_grp}')

        print('')
