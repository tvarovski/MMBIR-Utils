#!/usr/bin/sh

FILTER="False"
OUT_SUFFIX="_raw.txt"


while getopts "f" option; do
   case $option in
      f) #if set, filter out commons
         FILTER="True"
         OUT_SUFFIX="_allchrmmbir_filtered_all.txt"
         #exit;;
   esac
done



for sample in *-*/; do
    echo "inside $sample"
    for site in $sample*/; do
        echo "    analyzing ${site} ..."
        if [ $FILTER == "True" ]; then
          echo "using filtered setting"
          tools/filterAll ${site} ${sample}
        fi
        out1=${site%/}
        output_file=${out1#$sample}$OUT_SUFFIX

        python tools/annotate4.py $site $output_file $FILTER
    done
done

