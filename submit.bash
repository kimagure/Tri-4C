#!/bin/bash

s3fs -o nonempty tri4c /media/s3fs

source align_function.bash
source postalign.bash

fin(){
    path=$1
    indfile=$2
    s3fs=$3
execute $indfile $path
analysis $path
cp -r $path /media/s3fs
mv /media/s3fs/$path /media/s3fs/$path_`date |cut -d ' ' -f 1-3|sed 's/ /_/g'`
}
