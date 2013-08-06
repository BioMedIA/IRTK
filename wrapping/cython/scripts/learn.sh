#!/bin/bash

set -x
set -e

# training_patients="/home/kevin/Imperial/PhD/MyPHD/Detection/BOW/PIPELINE/metadata/training.tsv"
# original_folder="/home/kevin/Imperial/PhD/DATASETS/Originals/"
# clean_brainboxes="/home/kevin/Imperial/PhD/MyPHD/Detection/BOW/PIPELINE/metadata/clean_brainboxes.tsv"
# NEW_SAMPLING=1.0
# vocabulary="vocabulary.npy"
# vocabulary_step=2
# mser_detector="mser_detector"
# ga_file="/home/kevin/Imperial/PhD/MyPHD/Detection/BOW/PIPELINE/metadata/ga.csv"
# box_detector="box_detector"

training_patients="/vol/biomedic/users/kpk09/pipeline2/metadata/training.tsv"
original_folder="/vol/biomedic/users/kpk09/DATASETS/Originals/"
clean_brainboxes="/vol/biomedic/users/kpk09/pipeline2/metadata/list_boxes.tsv"
NEW_SAMPLING=0.8
vocabulary="vocabulary.npy"
vocabulary_step=2
mser_detector="mser_detector.linear"
ga_file="/vol/biomedic/users/kpk09/pipeline2/metadata/ga.csv"

# python ./create_bow.py \
#     --training_patients $training_patients \
#     --new_sampling $NEW_SAMPLING \
#     --original_folder $original_folder \
#     --step $vocabulary_step \
#     --output $vocabulary

python ./learn_mser.py \
    --training_patients $training_patients \
    --new_sampling $NEW_SAMPLING \
    --original_folder $original_folder \
    --clean_brainboxes $clean_brainboxes \
    --ga_file $ga_file \
    --vocabulary $vocabulary \
    --output $mser_detector


