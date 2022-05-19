FILTER=False
OUT_SUFFIX="_raw.txt"

while getopts ":f" option; do
   case $option in
      f) #if set, filter out commons
         $FILTER=True
         $OUT_SUFFIX="_allchrmmbir_filtered_all.txt"
         exit;;
   esac
done



for sample in ind*/; do
    echo "inside $sample"
    for site in $sample*/; do
        echo "    analyzing ${site} ..."

        tools/filterAll ${site} ${sample}

        out1=${site%/}
        output_file=${out1#$sample}$OUT_SUFFIX

        python tools/annotate3.py $site $output_file $FILTER
    done
done

