for n in {1..25..1}; do
	for f in config1404-chr${n}.txt ; do
		head -92 $f > temp1.txt
		cat new_config.txt >> temp1.txt
		tail -79 $f >> temp1.txt
		mv temp1.txt config2404-chr${n}.txt
	done
done
