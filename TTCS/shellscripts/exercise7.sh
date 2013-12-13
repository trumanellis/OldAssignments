#!/bin/bash
cd ~
FILES=`ls -a`
nfiles=0
nexec=0
for f in $FILES
do
    nfiles=`expr $nfiles + 1`
    if [ -f $f ]; then
        if [ -x $f ]; then
            nexec=`expr $nexec + 1`
        fi
    fi
    if [ -d $f ]; then
        echo `du -sh $f`
    fi
done
echo "$nfiles total files"
echo "$nexec are executable"
