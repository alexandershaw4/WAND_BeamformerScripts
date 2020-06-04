#!/bin/bash
# fix segmentation violation
# usage:
# fixsve.sh "$(pwd)/314_04037_VisualCut.ds/"

cp ~/Downloads/BadChannels $1
echo "Copied bad chan file"

NAME=`echo "$1" | cut -d'.' -f1`
EXT=`echo "$1" | cut -d'.' -f2`
NEW="${NAME}Fix.${EXT}"
echo "${NEW}"
newDs $1 $NEW
echo "Made new (fixed) dataset"

