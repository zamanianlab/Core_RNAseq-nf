#!/bin/bash

# echo core, thread, and memory
echo "CPU threads: $(grep -c processor /proc/cpuinfo)"
grep 'cpu cores' /proc/cpuinfo | uniq
echo $(free -g)

# set home () and mk dirs
export HOME=$PWD
mkdir data work output

# transfer data from staging
cp -r /staging/mzamanian/data/$1 data

# clone nextflow git repo
git clone https://github.com/zamanianlab/Core_RNAseq-nf.git

# run nextflow
export NXF_OPTS='-Xms1g -Xmx4g'
nextflow run Core_RNAseq-nf/WB-pe.nf -w work -c Core_RNAseq-nf/chtc.config --dir $1 --release "WBPS14" --species "brugia_malayi" --prjn "PRJNA10729" --rlen "150"

# rm files you don't want transferred back to /home/{username}
rm -r work
rm -r data

# mv large output files to staging to avoid their transfer back to /home/{username}
mv output/$1/ /staging/mzamanian/output/
