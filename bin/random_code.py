
"""this code will take columns from multiple files, ignore the header, and output them into a file."""
#awk -F '\t' 'BEGIN {OFS=FS}  FNR>1 {print $93, $33, $35, $2, $2, $2, $1}' part* > TCRB_GLIPH_preRA.csv

"""Unix code for adding a header to a file"""
#awk 'BEGIN {print "CDR3b\tTRBV\tTRBJ\tCDR3a\tTRAV\tTRAJ\tPatientCounts"} {print}' TCRB_GLIPH_preRA.csv > TCRB_GLIPH_preRA_header.txt

"""Python code to take columns from a dataframe, manipulate them, and place them back, update them."""
data = pd.read_csv("TCRB_GLIPH_preRA_header.txt", delimiter='\t', header=0)

dfv = data['TRBV'].str.split('*').str.get(0)
dfj = data['TRBJ'].str.split('*').str.get(0)
dfs = data['CDR3b'].str.replace(' ', "")
data.update(dfv)
data.update(dfj)
data.update(dfs)
data.to_csv('TCR_GLIPH_READY.csv', sep='\t', index=False)

"""GLIPH run code"""
#~/bin/gliph/bin/gliph-group-discovery.pl --tcr TCR_GLIPH_READY.csv --simdepth=10000 --lcminp=0.01
#different options
#gliph-group-discovery.pl --tcr mytcrtable.txt --simdepth=10000 --lcminp=0.001

#~/bin/gliph/bin/gliph-group-scoring.pl --convergence_file TCR_GLIPH_READY.csv-convergence-groups.txt \
                         --clone_annotations=TCR_GLIPH_READY.csv \
                         --motif_pval_file=TCR_GLIPH_READY.csv-kmer_resample_10000_minp0.01_ove10.txt

"""plot with seaborn"""
df = pd.read_csv('TCR_isoelectric_tissue.txt', header=0)

plt.figure()
sns.stripplot(x="specimen_tissue", y="cdr3_length.mean.median", data=df);
plt.show()


"""upset"""
import pyupset as pyu
df = pd.read_csv('TCR_isoelectric_tissue.txt', header=0)
pyu.plot(data_dict)
