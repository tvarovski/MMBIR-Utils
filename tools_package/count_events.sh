cd $1

for sample in *.txt; do
    lines=$(wc -l < "$sample")
    echo "$sample	$[ lines - 1]"
done

