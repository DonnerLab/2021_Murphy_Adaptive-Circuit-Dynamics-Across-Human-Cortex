#!/bin/sh

IDs=('DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'EXJ' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ')

for i in {0..17}
#for i in 0
do

echo ${IDs[i]}

docker run -ti --rm -v ~/Desktop/Experiments/Surprise_accumulation/Analysis/MRIs/fs_conv/${IDs[i]}:/input nben/occipital_atlas:latest

done

