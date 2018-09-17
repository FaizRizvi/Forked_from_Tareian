#!/usr/bin/env python

import sys
import argparse
import logging
import os.path
import pandas as pd
import numpy as np
import snippets
from datetime import datetime

"""Define the function main-will create parser, logging, load file, process"""
def main(arguments):
    # program options
    parser = argparse.ArgumentParser(description='calculate CDR3 length statistics from part-tables')
    parser.add_argument('part_tables', metavar='part-table', nargs='+', help='the part tables to process')
    parser.add_argument('--runs', '-r', metavar='R', nargs='+', help='only process reads from the given run(s)')
    parser.add_argument('--cutoffs', '-c', metavar='M', nargs='*', type=float, default=[30, 45, 60, 75], help='calculate the fraction of clones at CDR3 length cutoffs')
    # add log level arguments
    snippets.add_log_level_arg(parser)
    args = parser.parse_args(arguments)
    snippets.set_log_level(args)
    file_counter = 0
    # columns to load and their type
    column_types = {'participant_label':      object, #needed for part info
                    'amplification_template': object, #needed for in frame exclusion
                    'specimen_tissue':        object, #needed for tissue summary
                    'v_segment':              object, #needed for v-seg gene type
                    'amplification_template': object,
                    'run_label':              object,
                    'forward_primer':         object,
                    'igh_clone_id':          np.int64, #needed to group clones
                    'v_score':                np.float64, #needed for v-score exclusion
                    'v_j_in_frame':           np.bool, #needed for in frame exclusion
                    'cdr3_seq_nt_q':          object, #filter cdr3 length
                    'v_sequence':             object, #needed for sequence analysis
                    'cdr3_seq_nt_q':          object, #needed for sequence analysis
                    'cdr3_seq_nt_v':          object, #needed for sequence analysis
                    'isosubtype':          object} #needed for sequence analysis

    # load all the part tables
    cd3_len_data = []
    for filename in args.part_tables:
        file_counter = file_counter+1
        file_counter_total=len(args.part_tables)

        logging.info(str(datetime.now()) + ': '+ 'Loading data from: ' + filename + " " + str(file_counter) + '/' + str(file_counter_total))
        data = pd.read_table(filename, header=0, true_values=['t'], false_values=['f'], usecols=column_types.keys(), dtype=column_types)

        logging.info(str(datetime.now()) + ': '+ 'Excluding reads that are out of frame ')
        data = snippets.filter_out_of_frame(data)

        logging.info(str(datetime.now()) + ': '+ 'Excluding reads based on v-score')
        data = snippets.filter_v_score(data, min_v_score=50)

        logging.info(str(datetime.now()) + ': '+ 'Excluding reads that have no CDR3')
        data = snippets.filter_no_cdr3(data)

        logging.info(str(datetime.now()) + ': '+ 'Removing allele information from the v-segment labels')
        data = snippets.remove_allele(data)

        logging.info(str(datetime.now()) + ': '+ 'Calculating CDR3 length')
        data = snippets.calculate_cdr3_length(data)

        logging.info(str(datetime.now()) + ': '+ 'Defining Isotype')
        data = snippets.calculate_type(data)

        logging.info(str(datetime.now()) + ': '+ 'Calculating CDR3 hydrophobicity')
        data = snippets.calculate_cdr3_hydrophobicity(data)

        logging.info(str(datetime.now()) + ': '+ 'Calculating CDR3 isoelectric point')
        data = snippets.calculate_cdr3_isoelectric_point(data)

        logging.info(str(datetime.now()) + ': '+ 'Calculating CDR3 n-glycosylation')
        data = snippets.calculate_any_nglycosylation(data)

        logging.info(str(datetime.now()) + ': '+ 'Calculating CDR3 instability')
        data = snippets.calculate_cdr3_instability_index(data)

        logging.info(str(datetime.now()) + ': '+ 'Calculating CDR3 aromaticity')
        data = snippets.calculate_cdr3_aromaticity(data)

        logging.info(str(datetime.now()) + ': '+ 'Calculating v-segment mutation level')
        data = snippets.calculate_v_segment_mutation_level(data)

        cd3_len_data.append(data)

    # concatenate all the data into one big dataframe
    logging.info(str(datetime.now()) + ': '+ 'Concatenating Libraries')
    cd3_len_data = pd.concat(cd3_len_data)

    logging.info(str(datetime.now()) + ': '+ 'Outputing to CSV')
    cd3_len_data.to_csv(sys.stdout, columns=['participant_label', 'igh_clone_id',
                                               'specimen_tissue', 'cdr3_length',
                                                'cdr3_hydrophobicity', 'v_segment',
                                                'cdr3_isoelectric_point', 'cdr3_any_nglycosylation',
                                                'cdr3_instability', 'cdr3_aromaticity', 'isosubtype', 'type', 'cdr3_mut_freq'])

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
