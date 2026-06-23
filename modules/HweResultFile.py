import os
import pandas as pd

class HweResultFile:
    """
    Read HWE results from the input file and save them into memory

    Input file format is like the following, with a header line and tab-delimited columns. 
 
    VariantID       Naa     Nab     Nbb     Saa     Sab     Sbb     HetExc   Chi2    EChi2   PChi2
    rs0             1       53      1949    1.21    52.58   1949.21 1        1.0549  1.8864  3.97e-02
    rs1             306     738     959     333.29  683.41  986.29  0        61.673  111.99  7.3512
    ...
    """

    def __init__(self, filename):
        ### Member variables
        self.filename = filename     # Name of the input file
        self.delimiter = '\t'        # Delimiter of columns, tab or comma

        self.id_column = "VariantID"
        self.chi2_column = "Chi2"
        self.echi2_column = "EChi2"
        self.pchi2_column = "PChi2"
        self.het_column = "HetExc"

        self.is_valid = True
        self.err_msg = ""
        self.CheckFileFormat()

        self.chi_columns = [self.chi2_column, self.echi2_column, self.pchi2_column]

    def CheckFileFormat(self):
        """
        Check if the input file is valid
        """
        try:
            file = open(self.filename, "r")
        except Exception as error:
            print(f'\nERROR: Cannot open {self.filename}: {error}\n')
            exit()

        header = file.readline()

        vars = header.split(self.delimiter)
        num_vars = len(vars)

        exp_col_num = 11
        if num_vars != exp_col_num \
            or vars[0].strip().upper() != self.id_column.upper() \
            or vars[-3].strip().upper() != self.chi2_column.upper() \
            or vars[-2].strip().upper() != self.echi2_column.upper() \
            or vars[-1].strip().upper() != self.pchi2_column.upper() :

            self.is_valid = False
            self.err_msg += f"Error in #columns or 1st or last 3 column headers\n"

        if not self.is_valid:
            err_msg = f'\nERROR: The input file {self.filename} is not a valid HWE result file.\n' + \
                f'\nExpected a tab-delimited file with the following columns:\n' + \
                f'  {self.id_column}\n' + \
                f'  Naa\n' + \
                f'  Nab\n' + \
                f'  Nbb\n' + \
                f'  Saa\n' + \
                f'  Sab\n' + \
                f'  Sbb\n' + \
                f'  HetExc\n' + \
                f'  {self.chi2_column}\n' + \
                f'  {self.echi2_column}\n' + \
                f'  {self.pchi2_column}\n' + \
                f'\nError message: {self.err_msg}\n'

            print(err_msg)
            exit()

        self.id_column    = vars[0].strip()
        self.chi2_column  = vars[-3].strip()
        self.echi2_column = vars[-2].strip()
        self.pchi2_column = vars[-1].strip()

        file.close()
