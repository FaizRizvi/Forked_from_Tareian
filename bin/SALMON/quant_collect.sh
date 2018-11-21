#!/bin/bash
for fn in /data/miraldiLab/reference_data/Gaublomme/salmon_quants/*;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
cp /data/miraldiLab/reference_data/Gaublomme/salmon_quants/${samp}/quant.sf /data/miraldiLab/reference_data/Gaublomme/Gobluomme_quants/${samp}_quant.sf
done 