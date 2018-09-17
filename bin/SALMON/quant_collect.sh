#!/bin/bash
for fn in /data/miraldiLab/LCL_repo/reference_datasets/transcriptome/GEUVADIS_salmon_output/*;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
cp /data/miraldiLab/LCL_repo/reference_datasets/transcriptome/GEUVADIS_salmon_output/${samp}/quant.sf /data/miraldiLab/Tareian_workspace/GEUVADIS_salmon/${samp}_quant.sf
done 