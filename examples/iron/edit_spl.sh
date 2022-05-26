files=aspl_*
for i in $files
do
	awk 'NR > 10 {print}' $i > e${i}
done
