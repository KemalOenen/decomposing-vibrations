import argparse
from argparse import RawTextHelpFormatter


def get_args():
    """
    This function defines all the possible arguments of the nomodeco python script using the argparse package
    """
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
  #  parser.add_argument("output",
  #                      help='output file containing '
  #                           'geometry optimized coordinates and hessian')
    parser.add_argument("--log", action='store_true',
                        help='set if you want an additional .log file that contains the results for '
                             'every tested internal coordinate set')
    parser.add_argument('--matrix_opt', nargs='?', default='contr',
                        metavar='matrix', help='choose which matrix to use for optimization: i.e., '
                                               'VED matrix (keyword: ved), Diagonal elements of PED matrix '
                                               '(keyword: diag) and / or Contribution Table (keyword: contr) '
                                               '(default: %(default)s)')
    parser.add_argument('--penalty1', nargs='?', type=int,
                        default=0, metavar='INTFREQ-PENALTY', help='penalty value for asymmetric intrinsic'
                                                                   'frequency values, which can be helpful for cyclic '
                                                                   'systems'
                                                                   ' (default: %(default)s)')
    parser.add_argument('--penalty2', nargs='?', type=int,
                        default=0, metavar='INTFC-PENALTY', help='penalty value for unphysical contributions'
                                                                 ' per internal coordinate (default: %(default)s)')
    parser.add_argument('--heatmap', nargs='+', metavar='matrix',
                        help='return a heatmap for the specified matrix, i.e., VED matrix (keyword: ved), '
                             'Diagonal elements of PED matrix (keyword: diag) '
                             'and / or Contribution Table (keyword: contr)')
    parser.add_argument('--csv', nargs='+', metavar='matrix',
                        help='return a csv for the specified matrix, i.e., VED matrix (keyword: ved), '
                             'Diagonal elements of PED matrix (keyword: diag) '
                             'and / or Contribution Table (keyword: contr)')
    parser.add_argument('--latex_tab', action='store_true',
                         help='generate additional latex table file which displayes the Contribution Table')
       
    parser.add_argument('--molpro', nargs = 1, metavar=('file1'),
                         help='Molpro.out file can be specified and used for Nomodeco input')

    parser.add_argument('--gv', nargs=1, metavar=('file1'),
                        help='Gaussian.log file can be specified and used for Nomodeco input')

    parser.add_argument('--orca', nargs=1, metavar=('file1'),
                        help = 'Orca.property.txt file can be specified and used for Nomodeco input')
   
    parser.add_argument("--pymolpro", action="store_true",
                        help="Use pymolpro as an integration to automatically run molpro calculations and subsequently Nomodeco.py")

    parser.add_argument("--comb", type = int, default = 3, choices= range(1,4),
            help="""Adjust which coordinates get added for IC set generation (for hydrogen bonded clusters only):
                        1  --> all hydrogen bond coordinates  
                        2 --> all hydrogen bond coordinates + acceptor donor bonds and angles  
                        3 (default) --> all hydrogen bond coordinates + acceptor donor bonds, angles, dihedrals and oop""")

    args = parser.parse_args()
   

    return args
