#!/usr/bin/env python
from __future__ import print_function

import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("-m", "--matrix", dest='matrix', help="Input PWM", required=True)
parser.add_argument("-o", "--output", dest='output', help="Output directory", required=True)

args = parser.parse_args()

for i in args.matrix:
    pwm_basename = os.path.splitext(i)[0]
    pwm_outname = str(args.output) + "/" + pwm_basename + "_JASPAR.txt"
    
    df = pd.read_csv(str(i), header=0, index_col=0, sep="\t")
    
    df_trans = df.transpose()
    
    df_trans.to_csv(pwm_outname, sep="\t", index=False, header=False)