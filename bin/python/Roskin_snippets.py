#!/usr/bin/env python

import re
import sys
import itertools
from collections import OrderedDict
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import logging

"""This is a compilation of useful functions"""

igh_regions = ['cdr3']

run_of_dashes = re.compile('-+')

"""remove alleles from strings"""
def remove_allele(data):
    data['v_segment'] = data['v_segment'].str.split('*').str.get(0)

    return data

def remove_allele_igh(s):
    """Remove the allele from a V, D, J, or constant region call

    :param s: string with name of V, D, J, or constant region
    :return: V, D, J, or constant region name without allele designation
    """
    if type(s) is str or type(s) is unicode:    # if string or unicode string
        return s.split('*')[0]
    else:   # non-strings, e.g. NaNs, NULL, None are passed through
        return s

def capped_seq_mutation_percent(s, default=None):
    """Calculate the mutation level

    Calculates the fraction of mutated bases by calculating the fraction of capital letters vs. lower case. If no bases
    are found, the given default is returned instead.
    :param s: sequence of the V-segment
    :param default: values to return in case of no bases
    :return: fraction of mutated bases in s
    :type s: str
    :type default: str, NoneType
    :rtype float, NoneType
    """
    if type(s) is str or type(s) is unicode:
        diff = s.count('A') + s.count('C') + s.count('G') + s.count('T')
        same = s.count('a') + s.count('c') + s.count('g') + s.count('t')
        if same + diff == 0:
            return default
        else:
            return float(diff) / (same + diff)
    else:
        return default

def capped_seq_mutation_count(s, default=None):
    if type(s) is str or type(s) is unicode:
        diff = s.count('A') + s.count('C') + s.count('G') + s.count('T')
        return diff
    else:
        return default

def nt_sequence_hydrophobicity(s):
    if '-' in s:
        s = s.replace('-', '')

    if len(s) % 3 != 0:
        return None
    else:
        s = str(Seq(s).translate())
        return aa_sequence_hydrophobicity(s)

def aa_sequence_hydrophobicity(s):
    if '*' in s:
        s = s.replace('*', '')

    if 'X' in s:
        s = s.replace('X', '')

    if len(s) == 0:
        return None

    return ProteinAnalysis(s).gravy()

isoelectric_point_table = {
        'A': 6.00,
        'L': 5.98,
        'R': 10.76,
        'K': 9.74,
        'N': 5.41,
        'M': 5.74,
        'D': 2.77,
        'F': 5.48,
        'C': 5.05,
        'P': 6.30,
        'Q': 5.65,
        'S': 5.68,
        'E': 3.22,
        'T': 5.66,
        'G': 5.97,
        'W': 5.89,
        'H': 7.59,
        'Y': 5.66,
        'I': 6.02,
        'V': 5.96}

def nt_sequence_isoelectric_point(s):
    if '-' in s:
        s = s.replace('-', '')

    if len(s) % 3 != 0:
        return None
    else:
        s = str(Seq(s).translate())
        return aa_sequence_isoelectric_point(s)

def aa_sequence_isoelectric_point(s):
    if '*' in s:
        s = s.replace('*', '')

    if 'X' in s:
        s = s.replace('X', '')

    if len(s) == 0:
        return None

    return sum(isoelectric_point_table[c] for c in s)/len(s)

def nt_any_nglycosylation(s):
    if '-' in s:
        s = s.replace('-', '')

    if len(s) % 3 != 0:
        return None
    else:
        s = str(Seq(s).translate())
        return aa_any_nglycosylation(s)

n_linked_nglycosylation_re = re.compile('N[^P](S|T)')
def aa_any_nglycosylation(s):
    if '*' in s:
        s = s.replace('*', '')

    if 'X' in s:
        s = s.replace('X', '')

    if len(s) == 0:
        return None

    if n_linked_nglycosylation_re.search(s):
        return 1
    else:
        return 0

def nt_sequence_instability(s):
    if '-' in s:
        s = s.replace('-', '')

    if len(s) % 3 != 0:
        return None
    else:
        s = str(Seq(s).translate())
        return aa_sequence_instability(s)

def aa_sequence_instability(s):
    if '*' in s:
        s = s.replace('*', '')

    if 'X' in s:
        s = s.replace('X', '')

    if len(s) == 0:
        return None

    return ProteinAnalysis(s).instability_index()

def nt_sequence_aromaticity(s):
    if type(s) is str or type(s) is unicode:
        if '-' in s:
            s = s.replace('-', '')

        if len(s) % 3 != 0:
            return None
        else:
            s = str(Seq(s).translate())
            return aa_sequence_aromaticity(s)
    return None

def aa_sequence_aromaticity(s):
    if '*' in s:
        s = s.replace('*', '')

    if 'X' in s:
        s = s.replace('X', '')

    if len(s) == 0:
        return None

    return ProteinAnalysis(s).aromaticity()

def v_mutation_spectrum(row):
    A_to_C = 0
    A_to_G = 0
    A_to_T = 0

    C_to_A = 0
    C_to_G = 0
    C_to_T = 0

    G_to_A = 0
    G_to_C = 0
    G_to_T = 0

    T_to_A = 0
    T_to_C = 0
    T_to_G = 0

    for prefix in ['cdr1_seq_nt_', 'fr2_seq_nt_', 'cdr2_seq_nt_', 'fr3_seq_nt_']:
        q_sequence = row[prefix + 'q']
        v_sequence = row[prefix + 'v']

        if type(v_sequence) is str:
            assert type(q_sequence) is str
            for v_base, q_base in itertools.izip_longest(v_sequence, q_sequence):
                if   v_base == 'A':
                    if   q_base == 'C': A_to_C += 1
                    elif q_base == 'G': A_to_G += 1
                    elif q_base == 'T': A_to_T += 1
                elif v_base == 'C':
                    if   q_base == 'A': C_to_A += 1
                    elif q_base == 'G': C_to_G += 1
                    elif q_base == 'T': C_to_T += 1
                elif v_base == 'G':
                    if   q_base == 'A': G_to_A += 1
                    elif q_base == 'C': G_to_C += 1
                    elif q_base == 'T': G_to_T += 1
                elif v_base == 'T':
                    if   q_base == 'A': T_to_A += 1
                    elif q_base == 'C': T_to_C += 1
                    elif q_base == 'G': T_to_G += 1

    A_total =          A_to_C + A_to_G + A_to_T
    C_total = C_to_A +          C_to_G + C_to_T
    G_total = G_to_A + G_to_C +          G_to_T
    T_total = T_to_A + T_to_C + T_to_G

    if A_total == 0:
        A_to_C = float('nan') ; A_to_G = float('nan') ; A_to_T = float('nan')
    else:
        A_to_C = float(A_to_C) / A_total
        A_to_G = float(A_to_G) / A_total
        A_to_T = float(A_to_T) / A_total

    if C_total == 0:
        C_to_A = float('nan') ; C_to_G = float('nan') ; C_to_T = float('nan')
    else:
        C_to_A = float(C_to_A) / C_total
        C_to_G = float(C_to_G) / C_total
        C_to_T = float(C_to_T) / C_total

    if G_total == 0:
        G_to_A = float('nan') ; G_to_C = float('nan') ; G_to_T = float('nan')
    else:
        G_to_A = float(G_to_A) / G_total
        G_to_C = float(G_to_C) / G_total
        G_to_T = float(G_to_T) / G_total

    if T_total == 0:
        T_to_A = float('nan') ; T_to_C = float('nan') ; T_to_G = float('nan')
    else:
        T_to_A = float(T_to_A) / T_total
        T_to_C = float(T_to_C) / T_total
        T_to_G = float(T_to_G) / T_total

    result = (        A_to_C, A_to_G, A_to_T,
              C_to_A,         C_to_G, C_to_T,
              G_to_A, G_to_C,         G_to_T,
              T_to_A, T_to_C, T_to_G)
    return result


def calculate_type(data, column_name='type', add_collapsed=False, collapsed_suffix='all', add_mut_split=None, mutated_cutoff=0.01, compact=False):
    assert column_name not in data

    # add a type column that gives DNA (with frame) for gDNA and the isotype for RNA
    data[column_name] = data['isosubtype']
    # remove allele
    data[column_name] = data[column_name].map(remove_allele_igh)

    # set values for gDNA
    data.loc[(data['amplification_template'] == 'gDNA') & data['v_j_in_frame'], column_name] = 'DNA:in-frame'
    data.loc[(data['amplification_template'] == 'gDNA') & ~data['v_j_in_frame'], column_name] = 'DNA:out-of-frame'

    if add_collapsed:
        # add types fro all IgA and IgG sub-isotypes combined
        collaped_data = data[data[column_name].isin(['IGHA1', 'IGHA2',
                              'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4'])].copy()
        collaped_data[column_name] = collaped_data[column_name].map(lambda i: i[:4] + collapsed_suffix)

        data = pd.concat([data, collaped_data])

    if add_mut_split is not None:
        mut_split = data[data[column_name].isin(add_mut_split)].copy()
        mut_split[column_name] += mut_split['v_sequence'].map(capped_seq_mutation_percent).map(lambda i: 'mutated' if i > mutated_cutoff else 'unmutated')
        data = pd.concat([data, mut_split])

    if compact:
        del data['isosubtype']
        del data['amplification_template']

    return data

""" add log level flags to the given args parser """
def add_log_level_arg(parser):
    parser.add_argument('--log-level', metavar='L', default='INFO', help='logging level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
                 'debug', 'info', 'warning', 'error', 'critical'])


def calculate_cdr3_length(data, column_name='cdr3_length', compact=False):
    assert column_name not in data

    # calculate CDR3 length
    data[column_name] = data['cdr3_seq_nt_q'].apply(lambda s: np.nan if type(s) == float else len(s.replace('.', '').replace('-', '')))

    if compact:
        del data['cdr3_seq_nt_q']

    return data

def calculate_v_segment_mutation_level(data, column_suffix='_mut_freq', compact=False):
    # for each region
    for s in igh_regions:
        # if the nt sequence is there
        input_column = '%s_seq_nt_q' % s
        if input_column in data:
            output_column = s + column_suffix
            assert output_column not in data
            # calculate the mutation level for that region
            data[output_column] = data[input_column].map(capped_seq_mutation_percent)
            if compact:
                del data[input_column]

    # calculate the mutation level for the whole V segment
    output_column = 'all' + column_suffix
    assert output_column not in data
    data[output_column] = data['v_sequence'].map(capped_seq_mutation_percent)

    if compact:
        del data['v_sequence']

    return data

def calculate_inserts(data, column_name='bases_inserted', compact=False):
    assert column_name not in data

    for c in 'insertions_fr1', 'insertions_cdr1', 'insertions_fr2', 'insertions_cdr2', 'insertions_fr3':
        data[c].fillna(0.0, inplace=True)

    data[column_name] = data['insertions_fr1']  + \
                        data['insertions_cdr1'] + \
                        data['insertions_fr2']  + \
                        data['insertions_cdr2'] + \
                        data['insertions_fr3']

    if compact:
        for c in 'insertions_fr1', 'insertions_cdr1', 'insertions_fr2', 'insertions_cdr2', 'insertions_fr3':
            del data[c]

    return data

def calculate_indel_events(data, compact=False):
    # convert the nan's to empty strings
    for r in igh_regions:
        data[r + '_seq_nt_v'].fillna('', inplace=True)
        data[r + '_seq_nt_q'].fillna('', inplace=True)

    data['insert_events'] = (data['fr1_seq_nt_v'] + data['cdr1_seq_nt_v'] +
                             data['fr2_seq_nt_v'] + data['cdr2_seq_nt_v'] +
                             data['fr3_seq_nt_v']).apply(lambda s: len(re.split(run_of_dashes, s)) - 1)

    data['delete_events'] = (data['fr1_seq_nt_q'] + data['cdr1_seq_nt_q'] +
                             data['fr2_seq_nt_q'] + data['cdr2_seq_nt_q'] +
                             data['fr3_seq_nt_q']).apply(lambda s: len(re.split(run_of_dashes, s)) - 1)

    if compact:
        for r in igh_regions:
            del data[r + '_seq_nt_v']
            del data[r + '_seq_nt_q']

    return data

def calculate_cdr3_hydrophobicity(data, column_name='cdr3_hydrophobicity', compact=False):
    assert column_name not in data

    data[column_name] = data['cdr3_seq_nt_q'].map(nt_sequence_hydrophobicity)

    if compact:
        del data['cdr3_seq_nt_q']

    return data

def calculate_cdr3_isoelectric_point(data, column_name='cdr3_isoelectric_point', compact=False):
    assert column_name not in data

    data[column_name] = data['cdr3_seq_nt_q'].map(nt_sequence_isoelectric_point)

    if compact:
        del data['cdr3_seq_nt_q']

    return data

def calculate_any_nglycosylation(data, column_name='cdr3_any_nglycosylation', source='cdr3_seq_nt_q', compact=False):
    assert column_name not in data

    data[column_name] = data[source].map(nt_any_nglycosylation)

    if compact:
        del data[source]

    return data

def calculate_cdr3_aromaticity(data, column_suffix='_aromaticity', compact=False):
    for s in igh_regions:
        input_column = '%s_seq_nt_q' % s
        if input_column in data:
            output_column = s + column_suffix
            assert output_column not in data
            data[output_column] = data[input_column].map(nt_sequence_aromaticity)
            if compact:
                del data[input_column]

    return data

def calculate_cdr3_instability_index(data, column_name='cdr3_instability', compact=False):
    assert column_name not in data

    data[column_name] = data['cdr3_seq_nt_q'].map(nt_sequence_instability)

    if compact:
        del data['cdr3_seq_nt_q']

    return data

def calculate_v_mutation_spectrum(data, compact=False):
    (                        data['A_to_C_freq'], data['A_to_G_freq'], data['A_to_T_freq'], \
        data['C_to_A_freq'],                      data['C_to_G_freq'], data['C_to_T_freq'], \
        data['G_to_A_freq'], data['G_to_C_freq'],                      data['G_to_T_freq'], \
        data['T_to_A_freq'], data['T_to_C_freq'], data['T_to_G_freq']                     ) = zip(*data.apply(v_mutation_spectrum, axis=1))

    if compact:
        for prefix in ['cdr1_seq_nt_', 'fr2_seq_nt_', 'cdr2_seq_nt_', 'fr3_seq_nt_']:
            del data[prefix + 'v']
            del data[prefix + 'q']

    return data

""" filter reads that have no cdr3"""
def filter_no_cdr3(data):
    # drop reads that don't have a CDR3
    data = data[~data['cdr3_seq_nt_q'].isnull()]

    return data


""" filter reads that are out of frame"""
def filter_out_of_frame(data, keep_dna=True):
    if keep_dna:
        data = data[(data['amplification_template'] == 'gDNA') | data['v_j_in_frame']]
    else:
        data = data[data['v_j_in_frame']]

    return data


""" filter reads that are below the desired v-score """
def filter_v_score(data, min_v_score=140, compact=False):
    # only keep reads with a good enough V-score
    data = data[data['v_score'] >= min_v_score]

    if compact:
        del data['v_score']

    return data


def fraction_over(cutoff, name=None, name_prefix='frac_over_'):
    function = lambda i: (i > cutoff).astype(float).mean()

    if name is None:
        function.__name__ = name_prefix + str(cutoff)
    else:
        function.__name__ = name

    return function


def fraction_under(cutoff, name=None, name_prefix='frac_under_'):
    function = lambda i: (i <= cutoff).astype(float).mean()

    if name is None:
        function.__name__ = name_prefix + str(cutoff)
    else:
        function.__name__ = name

    return function


"""Takes a table and returns statistics on the groups:param data: dataframe with data: param subgrouping: additional subgroups to use:return: dataframe with grouping statistics"""
def groupby_summary(grouped_data, column_names, statistics=[np.mean, np.median, np.std, np.size],main_group=None, subgrouping=None, subgrouping2=None, subgrouping3=None, drop_index=True):
    summary = grouped_data.groupby(main_group + subgrouping + subgrouping2 + subgrouping3).aggregate(
        OrderedDict([(c, statistics) for c in column_names]))
    summary.columns = ['.'.join(i) for i in summary.columns]

    if drop_index:
        return summary.reset_index()

    else:
        return summary


"""Calculate percentiles"""
def percentile(p, name=None, name_prefix='percentile_'):
    function = lambda i: i.quantile(p)

    if name is None:
        function.__name__ = name_prefix + str(p)
    else:
        function.__name__ = name

    return function


""" set the log level using the args from argparse """
def set_log_level(args):
    log_level = getattr(logging, args.log_level.upper(), None)
    logging.basicConfig(level=log_level)
