#!/bin/bash

s1dir=$1

cd $2


site=${1#$2}


for d in */; do
    echo "preparing $d directory for filtering"
    for i in {1..25}; do
	if [[ $i = 26 ]]; then
	    continue
	else
	    cp ${d}chr$i/final_bir_locs.txt ${d}chr$i/non_commons_fromAll.txt
	fi
    done
done


for d in */; do
    echo $d
    echo $site
    if [ "$d" == "$site" ]; then
	echo "omitting"
	continue
    else
        echo "filtering from dir $d"
        for i in {1..25}; do
            if [[ $i = 26 ]]; then
                continue
            else
                python3 ../tools/findCommonsV6_all.py ${site}chr$i/non_commons_fromAll.txt ${d}chr$i/consensus_reads.txt
		mv commons_fromAlltemp.txt commons_fromAll.txt
		mv non_commons_fromAlltemp.txt non_commons_fromAll.txt 
                mv commons_fromAll.txt ${site}chr$i/
                mv non_commons_fromAll.txt ${site}chr$i/
            fi
        done
    #./getfilteredbirlocs.sh $1        
    fi
done
