#!/bin/bash

set -x
set -e

filename=$1
ga=$2
output_folder=$3

if [ -z "$output_folder" ];
then
    name=`basename $filename`
    name=(${name//./ })
    name=${name[0]}

    output_folder=$name
fi

mkdir -p ${output_folder}

tmp_mask=$output_folder/mser_`basename $filename`

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

python $DIR/fetalMask_detection.py $filename $ga $tmp_mask
python $DIR/fetalMask_segmentation.py --img $filename --ga $ga --mask $tmp_mask --output_dir $output_folder

