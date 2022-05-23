## Wrapper for analysis of ufo mutant microarray data with Limma
Limma_script=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/scripts/microarray_limma.r

ufo_microarrays=/home/312.6-Flo_Re/312.6.1-Commun/data/UFO_microarrays
ufo_inflor_dir=/home/312.6-Flo_Re/312.6.1-Commun/data/UFO_microarrays/ufo_inflor


## Rename files: ATGE_29 -> Col-0_inflor; ATGE_52 -> ufo_inflor
old_names=("ATGE_29_A2.cel" "ATGE_29_B2.cel" "ATGE_29_C2.cel"
		   "ATGE_52_A.cel" "ATGE_52_B.cel" "ATGE_52_C.cel")
new_names=("Col-0_inflor_A.cel" "Col-0_inflor_B.cel" "Col-0_inflor_C.cel"
		   "ufo_inflor_A.cel" "ufo_inflor_B.cel" "ufo_inflor_C.cel")
for i in {0..5}; do
	if [[ ! -f $ufo_inflor_dir/${new_names[$i]} ]]; then
		mv $ufo_inflor_dir/${old_names[$i]} $ufo_inflor_dir/${new_names[$i]}
	else
		echo "okay!"
	fi
done


############# Microarray analysis using Limma & volcano plots #############
## Limma microarray analysis and volcano plots for each:
# ufo_flow directory and then ufo_inflor:
if [[ ! -f $ufo_inflor_dir/ggvolcano.png ]]; then
	Rscript $Limma_script -n wt_ufo_inflor -i $ufo_inflor_dir -o $ufo_inflor_dir -r 3 -c Col -t ufo -s 3 &> $ufo_inflor_dir/microarray_processing.log
fi


## create txt files containing only list of ATXG of DEG genes in each exp:
## for FC+-0.5
awk -F ";" '{print $2}' $ufo_inflor_dir/wt_ufo_inflor_DEG.csv | sed 's/"//g' | grep -v "ATM" | grep -v "ATC" > $ufo_inflor_dir/wt_ufo_inflor_DEG_ATXG.txt
wc -l <$ufo_inflor_dir/wt_ufo_inflor_DEG_ATXG.txt

echo 'done!'


