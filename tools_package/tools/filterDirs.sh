for sample in ind*/; do
    echo "inside $sample"
    for site in $sample*/; do
        echo "    analyzing ${site} ..."

        tools/filterAll ${site} ${sample}

        out1=${site%/}
        #_raw.txt
        #_allchrmmbir_filtered_all.txt
        output_file=${out1#$sample}_allchrmmbir_filtered_all.txt

        python tools/annotate2.py $site $output_file
    done
done

