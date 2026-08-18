from modules.HweResultFile import HweResultFile

def is_number(str):
    """
    Check if a tring is a number
    """
    if str is None:
        return False

    try:
        float(str)
        return True
    except ValueError:
        return False


class GraphParameters:
    """
    Set parameters for the graphs to be plotted    
    """

    def __init__(self, args, hwe_file):
        self.plot_type = "histo"  # type of graph to be plotted, i.e., histogram or scatterplot.

        self.xvar = args.xvar     # score to be plotted on x-axis
        self.yvar = args.yvar     # score to be plotted on y-axis

        self.width = 14           # width of graph in inches
        self.height = 14          # height of graph in inches
        
        if args.scatter:
            self.plot_type = "scatter"

        self.CheckAxisVariables(hwe_file)

        self.disp_xvar = self.GetAxisDisplayValues(self.xvar, hwe_file)
        self.disp_yvar = self.GetAxisDisplayValues(self.yvar, hwe_file)
        
        min_graph_size, max_graph_size, def_graph_size = 3, 100, 10
        self.width  = self.SetValueLimits(args.gw, min_graph_size, max_graph_size, def_graph_size)
        self.height = self.SetValueLimits(args.gh, min_graph_size, max_graph_size, def_graph_size)
        
        chi2_min_val, chi2_max_val = -1000, 100000
        default_chi2_min, default_chi2_max   = -50, 1000
        default_pchi2_min, default_pchi2_max = -20, 200
        default_zval_min, default_zval_max   = -10, 10
        default_zval_bin = 0.1

        self.xmin  = self.SetValueLimits(args.xmin, chi2_min_val, chi2_max_val, default_chi2_min)
        self.xmax  = self.SetValueLimits(args.xmax, chi2_min_val, chi2_max_val, default_chi2_max)
        self.ymin  = self.SetValueLimits(args.ymin, chi2_min_val, chi2_max_val, default_chi2_min)
        self.ymax  = self.SetValueLimits(args.ymax, chi2_min_val, chi2_max_val, default_chi2_max)
        
        # When y min and max are not specified, set them to the same values as x min and max. This is to make sure that the scatter plot is square.
        if args.ymin == None or args.ymax == None:
            self.ymin = self.xmin
            self.ymax = self.xmax   

        # For a histogram, x-axis is chi2, y-axis is number of variants
        if self.plot_type == 'histo':
            min_cnts, max_cnts, default_cnts = 100, 10000000, 10000
            min_bin_size, max_bin_size, default_chi2_bin_size, default_pchi2_bin_size = 0.01, 1000, 5, 1

            default_bin_size = default_chi2_bin_size 
            default_min, default_max = default_chi2_min, default_chi2_max

            self.xmin, self.ymin = 0, 0     # Chi2 values and counts are all positive 

            if self.xvar == hwe_file.zval_column:
                # z-values can be negative
                self.xmin  = self.SetValueLimits(args.xmin, chi2_min_val, chi2_max_val, default_min)
                default_min = default_zval_min
                default_max = default_zval_max
                default_bin_size = default_zval_bin 
            elif self.xvar in [hwe_file.pchi2_column, hwe_file.dchi2_column]:
                default_min = default_pchi2_min
                default_max = default_pchi2_max
            
            self.bin_size = self.SetValueLimits(args.bin_size, min_bin_size, max_bin_size, default_bin_size)
            self.xmax  = self.SetValueLimits(args.xmax, chi2_min_val, chi2_max_val, default_max)
            self.ymax  = self.SetValueLimits(args.ymax, min_cnts, max_cnts, default_cnts)

#            print(f"bin {self.bin_size}")
#            exit()
        else:
            self.dot_size = self.SetValueLimits(args.dot_size, 1, 100, 5)

        if self.xmin >= self.xmax:
            print(f"\nERROR: wrong x-value range: {self.xmin} - {self.xmax}")
            exit()

        if self.ymin >= self.ymax:
            print(f"\nERROR: wrong y-value range: {self.ymin} - {self.ymax}")
            exit()

        self.plot_title = ""          # graph title
        if args.title:
            self.plot_title = args.title

        self.show_bin_cnts = args.show_bin_cnts
        self.chi_dist = args.chi_dist
        self.norm_dist = args.norm_dist
        self.het_column = hwe_file.het_column
        self.hetexc_only = args.hetexc_only
        self.hetdef_only = args.hetdef_only
        self.hetexc_color = args.hetexc_color
        self.hetdef_color = args.hetdef_color   
        self.title = args.title
        self.tight_plot = args.tight_plot
        self.nolimits = args.nolimits
        self.title_size = float(args.title_size)
        self.axis_lbl_size = float(args.axis_lbl_size)

        self.chi_adj = 0
        if args.chi_adj:
            try:
                self.chi_adj = float(args.chi_adj)
            except Exception as error:
                self.chi_adj = 0

        if (self.chi_adj < 0):
            self.chi_adj = 0

        # Currently not user-speficied
        self.dot_color = 'black'
        self.line_color = 'gray'
        self.line_width = 0.8
        self.tick_lbl_size = 12

        if self.xvar in (hwe_file.chi2_column, hwe_file.echi2_column, hwe_file.pchi2_column, hwe_file.dchi2_column):
            self.norm_dist = False

        if self.xvar in (hwe_file.zval_column):
            self.chi_dist = False

    def GetAxisDisplayValues(self, var, hwe_file):
        """
        Get display values on an axis
        """

        disp_var = var

        if var == hwe_file.chi2_column:
            disp_var = '$\\chi^2$'
        elif var == hwe_file.pchi2_column:
            disp_var = 'p$\\chi^2$'
        elif var == hwe_file.dchi2_column:
            disp_var = 'd$\\chi^2$'
        elif var == hwe_file.echi2_column:
            disp_var = 'Expected $\\chi^2$'

        return disp_var


    def CheckAxisVariables(self, hwe_file):
        """
        Make sure the user-selected axis variables are correct
        """

        error = ''
        if self.xvar not in hwe_file.chi_columns:
            error += f"X-axis variable {self.xvar} is not a valid column in the input file. \n"
            error += f"Expected variables are: {', '.join(hwe_file.chi_columns)}. \n"

        if self.plot_type == 'scatter':
            if self.yvar not in hwe_file.chi_columns:
                error += f"Y-axis variable {self.yvar} is not a valid column in the input file. \n"
                error += f"Expected variables are: {', '.join(hwe_file.chi_columns)}. \n"

        if error:
            print(f'\nERROR: {error}\n')
            exit()


    def SetValueLimits(self, val, min_val, max_val, default_val):
        """
        Bound a value to a range
        """

        new_val = val

        if val is None or not is_number(val):
            new_val = default_val
        else:
            new_val = float(val)

        if new_val < min_val:
            new_val = min_val
        elif new_val > max_val:
            new_val = max_val

        return new_val
  

    def ShowSummary(self):
        """
        Show all parameters
        """

        print(f"\nParameters:")
        print(f"\tPlot type:    {self.plot_type}")
        print(f"\tPlot size:    {self.width} x {self.height}")
        if self.plot_type == 'histo':
            print(f"\tScore:        {self.xvar}")
            print(f"\tShow range:   {self.xmin} - {self.xmax}")
            print(f"\tMax Count:    {self.ymax}")
            print(f"\tBin size:     {self.bin_size}")

            if self.chi_dist:
                print()
                print(f"\tShow Chi-square distribution: yes")
            if self.norm_dist:
                print()
                print(f"\tShow normal distribution: yes")
        else:
            print(f"\tX-variable:   {self.xvar}")
            print(f"\tY-variable:   {self.yvar}")
            print(f"\tX-range:      {self.xmin} - {self.xmax}")
            print(f"\tY-range:      {self.ymin} - {self.ymax}")
            print(f"\tDot size:     {self.dot_size}")
            if self.hetexc_color:
                print(f"\tHetExc color: {self.hetexc_color}")

            if self.hetdef_color:
                print(f"\tHetDef color: {self.hetdef_color}")

            if self.hetexc_only:
                print(f"\tShow only heterozygote excess variants: yes")
 
            if self.hetdef_only:
                print(f"\tShow only heterozygote deficient variants: yes")

                
 