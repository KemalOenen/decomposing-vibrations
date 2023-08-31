import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("output",
                        help='output file containing '
                             'geometry optimized coordinates and hessian')
    parser.add_argument('--penalty1',  nargs='?', type=int,
                        default=0, metavar='INTFREQ-PENALTY', help='penalty value for asymmetric intrinsic'
                                        ' frequency values, which can be helpful for cyclic systems'
                                                                   ' (default: %(default)s)')
    parser.add_argument('--penalty2',  nargs='?', type=int,
                        default=1, metavar='INTFC-PENALTY', help='penalty value for unphysical contributions'
                                        ' per internal coordinate (default: %(default)s)')
    args = parser.parse_args()
    return args
