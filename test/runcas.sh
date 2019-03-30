#!/bin/bash
#runcas
cd /Users/yifan/MiSeq/data/MiSeq20190218/R1/p53/txt/all/
for input in ./*.txt; do
    output=$input".output.txt"
    echo $input
    /Users/yifan/MiSeq/test/cas-offinder $input G $output
    echo The output for above file is generated
    echo $output
done