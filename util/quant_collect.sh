#!/bin/bash
for fn in /Users/caz3so/Dropbox/data/salmon_quants/quants/*_quant;
do
samp=`basename ${fn} _quant`
echo "Processing sample ${samp}"
cp /Users/caz3so/Dropbox/data/salmon_quants/quants/${samp}_quant/quant.sf /Users/caz3so/Dropbox/data/salmon_collected_quants/${samp}_quant.sf
done 