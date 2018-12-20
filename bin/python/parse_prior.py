#!/usr/bin/env python
from __future__ import print_function

import pandas as pd
import seaborn as sns
import os
import sys
import argparse
import pybedtools
import numpy as np
import snippets
parser = argparse.ArgumentParser()

parser.add_argument("-p", "-PRIOR", dest='IN_PRIOR', help="Input PRIOR file: sep = tab", required=True, default=[])

# set the arguments from the command line to variables in the args object
args = parser.parse_args()
