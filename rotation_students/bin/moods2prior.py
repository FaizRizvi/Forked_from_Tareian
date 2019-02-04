@@ -14,11 +14,10 @@ import snippets
# Set up the argument parser with all of the options and defaults set up.
parser = argparse.ArgumentParser()

parser.add_argument("-m", "--MOODS", dest='IN_MOODS', help="Input MOODS file: sep = |", required=True)
parser.add_argument("-b", "--BED", dest='IN_BED', help="Input BED file", required=True)
parser.add_argument("-t", "--TF", dest='IN_TF', help="Input TF information file: sep = \t", required=True)
parser.add_argument("-T", "--TSS", dest='TSS', help="Input TSS information file: sep = \t", required=True)
parser.add_argument("-g", "--GENE", dest='GENE', help="Input GENE information file: sep = \t", required=True)
parser.add_argument("-d", "--DOMAIN", dest='DOMAIN', help="Input TAD information file: sep = \t", required=True)
parser.add_argument("-o", "--OUT_DIR", dest='OUT_DIR', help="Output dir", required=True)

# set the arguments from the command line to variables in the args object
@ -28,11 +27,11 @@ args = parser.parse_args()
if not os.path.exists(args.OUT_DIR):
        os.makedirs(args.OUT_DIR)

# Change output directory 
# Change output directory
os.chdir(args.OUT_DIR)

# Get the basename of the files
base = os.path.basename(args.IN_MOODS)
base = os.path.basename(args.IN_BED)
basename = os.path.splitext(base)[0]

# Set up the file names based on the basename
@ -50,22 +49,24 @@ print("Loading and Parsing Gene Information")
dict_gene = snippets.genemeta2dict(args.GENE)

##########################-----------LOAD MOODS DATA------------##############################
print("Loading and Parsing MOODS Information")

bed_moods = snippets.parse_moods(args.IN_MOODS, dict_TF)
print("Loading and Parsing BED Information")
#replace with the actual bedfile
#bed_moods = snippets.parse_moods(args.IN_MOODS, dict_TF)
#bed_moods = snippets.parse_moods(args.IN_BED, dict_TF)
bed_moods = pybedtools.BedTool(args.IN_BED)

##########################-----------LOAD TSS DATA------------##############################
print("Loading TSS file and converting to BED")

bed_TSS = pybedtools.BedTool(args.TSS)

##########################-----------Intersect 10kb of TSS------------##############################
print("Intersecting BEDS")

bed_TSS_and_bed_moods = bed_TSS.window(bed_moods, w=10000)

columns = ["Gene_chr", "Gene_start", "Gene_end", "Gene_name", "Gene_strand", "TF_chr", "TF_start", "TF_end", "TF_name", "PEAK_ID"]
keep = ["Gene_name", "TF_name", "PEAK_ID"]
columns = ["Gene_chr", "Gene_start", "Gene_end", "Gene_name", "Gene_strand", "TF_chr", "TF_start", "TF_end", "TF_name"]

keep = ["Gene_name", "TF_name", "Gene_chr"]

df_intersect = bed_TSS_and_bed_moods.to_dataframe(names=columns, usecols=keep)
