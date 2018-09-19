#!/usr/bin/env python
from __future__ import print_function

import MOODS.scan
import MOODS.tools
import MOODS.parsers

import os
import sys
import argparse
from itertools import groupby, chain
import pandas as pd
import seaborn as sns
import pybedtools

##########################-----------ARGUMENT PARSER------------##############################
#Set up the argument parser with all of the options and defaults set up.
##############################################################################################

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbosity", action="count", default=0, help='verbose (-vv, -vvv for more)')
parser.add_argument("-b", "--BED",  help="Input DHS Bed file")
parser.add_argument("-z", "--MOODS", help="Run Moods", action="store")
parser.add_argument('-T','--tf', action='store', dest='tf_name_file', help='A tab-delimited file that has TF Name and Motif Name)')

#MOODS OPTIONS - input files
input_group = parser.add_argument_group("input files (at least one matrix and sequence file required)")
input_group.add_argument('-m','--matrices', metavar='M', nargs='+', action='store', dest='matrix_files', help='matrix files (count/frequency, will be converted to PWM before matching)', default = [])
input_group.add_argument('-S','--score-matrices', metavar='M', nargs='+', action='store', dest='lo_matrix_files', help='matrix files (PWM/other score matrix, will be matched directly)', default = [])
input_group.add_argument('-s','--sequences', metavar='S', nargs='+', action='store', dest='sequence_files', help='sequence files', default = [])

#MOODS OPTIONS - thresholds
th_group = parser.add_argument_group("threshold selection (exactly one required)")
th_group.add_argument('-p','--p-value', metavar='p', action='store', dest='p_val', type=float, help='compute threshold from p-value')
th_group.add_argument('-t','--threshold', metavar='T', action='store', dest='t', type=float, help='use specified absolute threshold')
th_group.add_argument('-B','--best-hits', metavar='n', action='store', dest='max_hits', type=int, help='return at least the specified amount of best matches')

#MOODS OPTIONS - output
out_group = parser.add_argument_group("output (optional)")
out_group.add_argument('-o','--outfile', metavar='outfile', action='store', dest='output_file', help='output to file instead of standard output', required=True)
out_group.add_argument('--sep', metavar='S', action='store', dest='sep', default=",", help='set field separator in out (default ",")')

#MOODS OPTIONS - behaviour
option_group = parser.add_argument_group("search and model behaviour (optional)")
option_group.add_argument('-R','--no-rc', dest='no_rc', action="store_true", help='disable matching versus reverse complement strand')
option_group.add_argument('--no-snps', dest='no_snps', action="store_true", help='ignore IUPAC symbols coding multiple nucleotides')
option_group.add_argument('--batch', dest='batch', action="store_true", help='do not recompute thresholds from p-values for each sequence separately (recommended when dealing with lots of short sequences)', default = "batch")
option_group.add_argument('--bg', metavar=('pA', 'pC', 'pG', 'pT'), nargs=4, action='store', type=float, dest='bg', default = [0.25,0.25,0.25,0.25], help='background distribution for computing thresholds from p-value with --batch (default is 0.25 for all alleles)')
option_group.add_argument('--ps', metavar='p', action='store', dest='ps', type=float, help='specify pseudocount for log-odds conversion (default = 0.1)', default = 0.01)# bg
option_group.add_argument('--log-base', metavar='x', action='store', dest='log_base', type=float, help='logarithm base for log-odds conversion (default natural logarithm)')
option_group.add_argument('--lo-bg', metavar=('pA', 'pC', 'pG', 'pT'), nargs=4, action='store', type=float, dest='lo_bg', default = [0.25,0.25,0.25,0.25], help='background distribution for log-odds conversion (default is 0.25 for all alleles)')

#CHIP OPTIONS
parser.add_argument("-c", "--CHIP", help="Folder containing ChIP BED files", action="store_true")

#set the arguments from the command line to variables in the args object
args = parser.parse_args()

##########################-----------CHECKPOINT 1------------#################################
# this If statement will determine whether MOODS needs to be run based on the input arguments.
##############################################################################################

if args.MOODS == "T":

	###########################################------MOODS-------#############################################
	###### This was lifted from the dna_MOODS.py file that came with MOODS 1.9.3. All I did was copy and paste
	###### and integrate into my pipeline.
	##########################################################################################################

	# sanity check for parameters, apply effects

	if len(args.matrix_files) == 0 and len(args.lo_matrix_files) == 0:
		print("{}: error: no matrix files given (use -h for help)".format(os.path.basename(__file__)), file=sys.stderr)
		sys.exit(1)
	if len(args.sequence_files) == 0:
		print("{}: error: no sequence files given (use -h for help)".format(os.path.basename(__file__)), file=sys.stderr)
		sys.exit(1)
	if (args.p_val is None) and (args.t is None) and (args.max_hits is None):
		print("{}: error: no threshold given (use -h for help)".format(os.path.basename(__file__)), file=sys.stderr)
		sys.exit(1)
	if (args.p_val is not None and args.t is not None) or (args.t is not None and args.max_hits is not None) or (args.p_val is not None and args.max_hits is not None):
		print("{}: error: only one threshold specification allowed (use -h for help)".format(os.path.basename(__file__)), file=sys.stderr)
		sys.exit(1)
	if args.max_hits is not None and args.batch:
		print("{}: warning: ignoring --batch when used with -B".format(os.path.basename(__file__)), file=sys.stderr)
	if args.t is not None and args.batch:
		print("{}: warning: --batch is redundant when used with -t".format(os.path.basename(__file__)), file=sys.stderr)


	if args.log_base is not None and (args.log_base <= 1):
			print("{}: error: --log-base has to be > 1".format(os.path.basename(__file__)), file=sys.stderr)
			sys.exit(1)

	# redirect output
	output_target = sys.stdout
	if args.output_file is not None:
		try:
			outfile = open(args.output_file, 'w')
		except:
			print("{}: error: could not open output file {} for writing".format(os.path.basename(__file__), args.output_file), file=sys.stderr)
			sys.exit(1)
		output_target = outfile


	# --- Helper functions for IO ---

	def pfm_to_log_odds(filename):
		if args.log_base is not None:
			mat = MOODS.parsers.pfm_to_log_odds(filename, args.lo_bg, args.ps, args.log_base)
		else:
			mat = MOODS.parsers.pfm_to_log_odds(filename, args.lo_bg, args.ps)
		if len(mat) != 4:
			return False, mat
		else:
			return True, mat

	def adm_to_log_odds(filename):
		if args.log_base is not None:
			mat = MOODS.parsers.adm_to_log_odds(filename, args.lo_bg, args.ps, 4, args.log_base)
		else:
			mat = MOODS.parsers.adm_to_log_odds(filename, args.lo_bg, args.ps, 4)
		if len(mat) != 16:
			return False, mat
		else:
			return True, mat

	def pfm(filename):
		mat = MOODS.parsers.pfm(filename)
		if len(mat) != 4:
			return False, mat
		else:
			return True, mat

	def adm(filename):
		mat = MOODS.parsers.adm_1o_terms(filename, 4)
		if len(mat) != 16:
			return False, mat
		else:
			return True, mat

	def read_text_or_fasta(filename):
		with open(filename, "r") as f:
			line = ""
			while len(line) == 0:
				line = f.next().strip()
		if line[0] == '>':
			return iter_fasta(filename)
		else:
			return iter_text(filename)

	def iter_text(filename):
		with open(filename, "r") as f:
			seq = "".join([line.strip() for line in f])
		yield os.path.basename(filename), seq


	def iter_fasta(filename):
		with open(filename, "r") as f:
			iterator = groupby(f, lambda line: line[0] == ">")
			for is_header, group in iterator:
				if is_header:
					header = group.next()[1:].strip()
				else:
					yield header, "".join(s.strip() for s in group)

	# format hit sequence with variants applied
	# only for for snps right now
	def modified_seq(seq, start, end, applied_variants, variants):
		ret = list(seq[start:end].lower())
		for i in applied_variants:
			# print(ret)
			# print(variants[i].start_pos - start)
			ret[variants[i].start_pos - start] = (variants[i].modified_seq).upper()
		return "".join(ret)


	# print results
	def print_results(header, seq, matrices, matrix_names, results, results_snps=[], variants=[]):
		if args.no_rc:
			fr = results
			frs = results_snps
			rr = [[] for m in matrices]
			rrs = [[] for m in matrices]
		else:
			fr = results[:len(matrix_names)]
			frs = results_snps[:len(matrix_names)]
			rr = results[len(matrix_names):]
			rrs = results_snps[len(matrix_names):]

		# mix the results together, use + and - to indicate strand
		results = [ [(r.pos, r.score, '+', ()) for r in fr[i]] + [(r.pos, r.score, '-', ()) for r in rr[i]] + [(r.pos, r.score, '+', r.variants) for r in frs[i]] + [(r.pos, r.score, '-', r.variants) for r in rrs[i]] for i in xrange(len(matrix_names))]
		for (matrix, matrix_name,result) in zip(matrices, matrix_names, results):
			if args.verbosity >= 2:
				print("{}: {}: {} matches for {}".format(os.path.basename(__file__), header, len(result), matrix_name), file=sys.stderr)
			if len(matrix) == 4:
				l = len(matrix[0])
			if len(matrix) == 16:
				l = len(matrix[0]) + 1
			for r in sorted(result, key=lambda r: r[0]):
				strand = r[2]
				pos = r[0]
				hitseq = seq[pos:pos+l]
				if len(r[3]) >= 1:
					hitseq_actual = modified_seq(seq, pos, pos+l, r[3], variants)
				else:
					hitseq_actual = ""
				print(args.sep.join([header,matrix_name, str(pos), strand, str(r[1]), hitseq, hitseq_actual]), file=output_target)

	# --- Load matrices ---
	if args.verbosity >= 1:
		print("{}: reading matrices ({} matrix files)".format(os.path.basename(__file__), len(args.matrix_files)+len(args.lo_matrix_files)), file=sys.stderr)

	matrix_names = map(os.path.basename, [n for n in args.matrix_files + args.lo_matrix_files])
	matrices = []
	matrices_rc = []
	for filename in args.matrix_files:
		valid = False
		if filename[-4:] != '.adm': # let's see if it's pfm
			valid, matrix = pfm_to_log_odds(filename)
			if valid and args.verbosity >= 3:
				print("{}: read matrix file {} as pfm (converted to PWM)".format(os.path.basename(__file__), filename), file=sys.stderr)
		if not valid: # well maybe it is adm
			valid, matrix = adm_to_log_odds(filename)
			if valid and args.verbosity >= 3:
				print("{}: read matrix file {} as adm (converted to high-order PWM)".format(os.path.basename(__file__), filename), file=sys.stderr)
		if not valid:
			print("{}: error: could not parse matrix file {}".format(os.path.basename(__file__), filename), file=sys.stderr)
			sys.exit(1)
		else:
			matrices.append(matrix)
			matrices_rc.append(MOODS.tools.reverse_complement(matrix,4))
	for filename in args.lo_matrix_files:
		valid = False
		if filename[-4:] != '.adm': # let's see if it's pfm
			valid, matrix = pfm(filename)
			if valid and args.verbosity >= 3:
				print("{}: read matrix file {} as pfm".format(os.path.basename(__file__), filename), file=sys.stderr)
		if not valid: # well maybe it is adm
			valid, matrix = adm(filename)
			if valid and args.verbosity >= 3:
				print("{}: read matrix file {} as adm".format(os.path.basename(__file__), filename), file=sys.stderr)
		if not valid:
			print("{}: error: could not parse matrix file {}".format(os.path.basename(__file__), filename), file=sys.stderr)
			sys.exit(1)
		else:
			matrices.append(matrix)
			matrices_rc.append(MOODS.tools.reverse_complement(matrix,4))

	if args.no_rc:
		matrices_all = matrices
	else:
		matrices_all = matrices + matrices_rc

	# --- Scanning ---
	if args.p_val is not None or args.t is not None:
		scanner = None
		# build scanner if using fixed threshold
		if args.t is not None:
			thresholds = [args.t for m in matrices_all]
			scanner = MOODS.scan.Scanner(7)
			bg = args.bg # used for search optimisation only
			scanner.set_motifs(matrices_all, bg, thresholds)
		# build scanner if using p-value and --batch
		if args.p_val is not None and args.batch:
			bg = args.bg
			if args.verbosity >= 1:
				print("{}: computing thresholds from p-value".format(os.path.basename(__file__)), file=sys.stderr)
			thresholds = [MOODS.tools.threshold_from_p(m,bg,args.p_val,4) for m in matrices_all]
			if args.verbosity >= 3:
				for (m, t) in zip(matrix_names, thresholds):
					print("{}: threshold for {} is {}".format(os.path.basename(__file__), m, t), file=sys.stderr)
			if args.verbosity >= 1:
				print("{}: preprocessing matrices for scanning".format(os.path.basename(__file__)), file=sys.stderr)
			scanner = MOODS.scan.Scanner(7)
			bg = args.bg
			scanner.set_motifs(matrices_all, bg, thresholds)

		for seq_file in args.sequence_files:
			if args.verbosity >= 1:
				print("{}: reading sequence file {}".format(os.path.basename(__file__), seq_file), file=sys.stderr)
			try:
				seq_iterator = read_text_or_fasta(seq_file)
			except Exception:
				print("{}: error: could not parse sequence file {}".format(os.path.basename(__file__), seq_file), file=sys.stderr)
				sys.exit(1)

			for header, seq in seq_iterator:
				# preprocessing for the new sequence if using p-values and not --batch
				if args.p_val is not None and not args.batch:
					bg = MOODS.tools.bg_from_sequence_dna(seq,1)
					if args.verbosity >= 3:
						print("{}: estimated background for {} is {}".format(os.path.basename(__file__), header, bg), file=sys.stderr)
					if args.verbosity >= 1:
						print("{}: computing thresholds from p-value for sequence {}".format(os.path.basename(__file__), header), file=sys.stderr)
					thresholds = [MOODS.tools.threshold_from_p(m,bg,args.p_val,4) for m in matrices_all]
					if args.verbosity >= 3:
						for (m, t) in zip(matrix_names, thresholds):
							print("{}: threshold for {} is {}".format(os.path.basename(__file__), m, t), file=sys.stderr)
					if args.verbosity >= 1:
						print("{}: preprocessing matrices for scanning".format(os.path.basename(__file__)), file=sys.stderr)
					scanner = MOODS.scan.Scanner(7)
					scanner.set_motifs(matrices_all, bg, thresholds)

				if args.verbosity >= 1:
					print("{}: scanning sequence {}".format(os.path.basename(__file__), header), file=sys.stderr)
				results = scanner.scan(seq)
				if args.no_snps:
					results_snps = [[]]*len(matrices_all)
					snps = []
				else:
					snps = MOODS.tools.snp_variants(seq)
					if args.verbosity >= 2:
						 print("{}: {} variants for sequence {}".format(os.path.basename(__file__), len(snps), header), file=sys.stderr)
					results_snps = scanner.variant_matches(seq,snps)
				print_results(header, seq, matrices, matrix_names, results, results_snps, snps)

	# --- Scanning (max hits version) ---
	elif args.max_hits is not None:
		if args.no_rc:
			N = args.max_hits
		else:
			# this is a bit tricky
			N = (args.max_hits * 2)/3
		for seq_file in args.sequence_files:
			if args.verbosity >= 1:
				print("{}: reading sequence file {}".format(os.path.basename(__file__), seq_file), file=sys.stderr)
			try:
				seq_iterator = read_text_or_fasta(seq_file)
			except Exception:
				print("{}: error: could not parse sequence file {}".format(os.path.basename(__file__), seq_file), file=sys.stderr)
				sys.exit(1)

			for header, seq in seq_iterator:
				if args.verbosity >= 1:
					print("{}: scanning sequence {}".format(os.path.basename(__file__), header), file=sys.stderr)
				results = MOODS.scan.scan_best_hits_dna(seq, matrices_all, N)
				if args.no_snps:
					results_snps = [[]]*len(matrices_all)
					snps = []
				else:
					snps = MOODS.tools.snp_variants(seq)
					if args.verbosity >= 2:
						 print("{}: {} variants for sequence {}".format(os.path.basename(__file__), len(snps), header), file=sys.stderr)
					snp_thresholds = [min((r.score for r in rs)) for rs in results]
					if not args.no_rc:
						snp_thresholds = 2 * map(min, zip(snp_thresholds[:len(matrix_names)], snp_thresholds[len(matrix_names):]))
					scanner = MOODS.scan.Scanner(7)
					bg = MOODS.tools.bg_from_sequence_dna(seq,1)
					scanner.set_motifs(matrices_all, bg, snp_thresholds)
					results_snps = scanner.variant_matches(seq,snps)
				print_results(header, seq, matrices, matrix_names, results, results_snps, snps)

	if args.output_file:
		outfile.close()

else:
		print ("Not in the MOODs...")

MOTIF_HITS = args.output_file

####################################------------CHECKPOINT 2-----------------###################################
# From here on the script will process the motif file that was produced from MOODS or provided as the -o option.
################################################################################################################

# The TF file must look something ilke the following with a header column that matches the following column IDs.
#
# Motif_ID		TF_Name
# M00116_1.97d	TFAP2B

#Set the variable to the argument that was given for the file name
tf_name_file = args.tf_name_file

# create a dataframe from the TF_Name list
TF_df = pd.read_csv(tf_name_file, sep="\t", header=0)

# Manipulate the list to get the selected feature
TF_names_df = pd.DataFrame()
TF_names_df["Motif"] = TF_df["Motif_ID"]
TF_names_df["TF"] = TF_df["TF_Name"]
TF_names_df.set_index("Motif")

# Crete a dictionary of motif names and tf names
dicted = dict(zip(TF_names_df.Motif, TF_names_df.TF))

def parse_moods(x):
	df = pd.read_csv(x, header=None, sep="|")

	df[1] = df[1].str.replace('_JASPAR.txt.pfm', '')

	df.drop(6, axis=1, inplace=True)

	df_tmp1 = df[0].str.split(":", expand=True)
	df_tmp2 = df_tmp1[1].str.split("-", expand=True)

	df.drop(columns=0, inplace=True)

	df["chr"] = df_tmp1[0]
	df["start"] = df_tmp2[0]
	df["stop"] = df_tmp2[1]

	df.columns = ["MOTIF_ID", "TF_POS", "STRAND", "MATCH_SCORE", "MOTIF_SEQ", "chr", "start", "stop"]

	df["TF_start"] = df["start"].apply(int) + 1 + df["TF_POS"]
	df["TF_end"] = df["TF_start"] + df["MOTIF_SEQ"].str.len() - 1
	df["PEAK_ID"] = df["chr"] + "_" + df["start"].map(str) + "_" + df["stop"].map(str)
	df["MOTIF_POS"] = df["chr"] + "_" + df["TF_start"].map(str) + "_" + df["TF_end"].map(str)
	df["MOTIF_LEN"] = df["TF_end"] - df["TF_start"] + 1
	df["TF_Name"] = df["MOTIF_ID"].map(dicted)

	return df

y = parse_moods(MOTIF_HITS)

MOTIF_HITS_BED = MOTIF_HITS.replace(".txt", ".bed")

y.to_csv(MOTIF_HITS_BED, sep="\t", index=False)
