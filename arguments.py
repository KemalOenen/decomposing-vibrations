import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("output",
                        help='output file containing '
                             'geometry optimized coordinates and hessian')
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
    args = parser.parse_args()
    return args
