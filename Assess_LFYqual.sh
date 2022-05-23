main_dir=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY-UFO/results

scores_prog=/home/312.6-Flo_Re/312.6.1-Commun/scripts/TFgenomicsAnalysis/bin/scores.py

############################### Assess LFYBS quality in dLUBS
name1=LFYamp
name2=LFYUFO_cf
# table=$main_dir/comparisons/${name1}_${name2}/table_${name1}_${name2}.tsv

table=$main_dir/comparisons/${name1}_${name2}/table_${name1}_${name2}_RIP_CFC.tsv



echo "old dlubs matrix requested"
dlubs_matrice=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY-UFO/results/motifs/dLUBS/LFYUFO_cf/LFYUFO_cf.pfm
dirout=$main_dir/Assess_LFYqual/LFY_in_dLUBS
LFY_matrice=/home/312.6-Flo_Re/312.6.1-Commun/data/LFY.pfm



mkdir -p $dirout


## take 50 bp around peak maximum and sort
awk -v OFS="\t" '{print $1,$2+int(($3-$2)/2)-25,$2+int(($3-$2)/2)+25}' $table > $dirout/LU_peak.bed
# head $dirout/LU_peak.bed
tair10=/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas
sort -k1,1 -k2,3n $table > $dirout/table_sorted.bed



## dLUBS scores calculation on central 50 bp
if [[ ! -f $dirout/LU_peak.fas.scores ]]; then
	bedtools getfasta -fi $tair10 -fo $dirout/LU_peak.fas -bed $dirout/LU_peak.bed
	python $scores_prog -f $dirout/LU_peak.fas -m $dLUBS_matrice -o $dirout/
fi


## select best dLUBS score per peak and cut around dLUBS matrix
if [[ ! -f $dirout/dLUBS_bestscore_info.bed ]]; then 
	sort -k1,1 -k8,8nr $dirout/LU_peak.fas.scores \
			| awk -v OFS="\t" '{print $2,$4,$8,$1}' \
			| uniq -f3 \
			| awk -v OFS="\t" '{print $4,$1,$2,$3}' > $dirout/LU_peak_dLUBS_bestscore.bed
	
	
# 	
	if [[ $new == "yes" ]]; then
		sort -u -k1,1 $dirout/LU_peak_dLUBS_bestscore.bed | sed 's/-/\t/' | sed 's/:/\t/g' | awk -v FS="\t" -v OFS="\t" '{print $1,$2+$4-1,$2+$4+25-1,$5}' | sort -k1,1 -k2,3n > $dirout/dLUBS_bestscore_info.bed ## new
	else
		sort -u -k1,1 $dirout/LU_peak_dLUBS_bestscore.bed | sed 's/-/\t/' | sed 's/:/\t/g' | awk -v FS="\t" -v OFS="\t" '{print $1,$2+$4-1,$2+$4+26-1,$5}' | sort -k1,1 -k2,3n > $dirout/dLUBS_bestscore_info.bed ## old
	fi
	
fi


## retrieve FASTA sequences and and double-check dLUBS coordinates 
if [[ ! -f $dirout/dLUBS_bestscore_info.fas.scores ]]; then
	bedtools getfasta -fi $tair10 -fo $dirout/dLUBS_bestscore_info.fas -bed $dirout/dLUBS_bestscore_info.bed
	python $scores_prog -f $dirout/dLUBS_bestscore_info.fas -m $dLUBS_matrice -o $dirout/
fi


## shift by 11 bp on the left and extend by 4 bp on the right --> LFY matrix
if [[ $new == "yes" ]]; then
	awk -v OFS="\t" '$4=="+"{print $1,$2+9,$3+3}$4=="-"{print $1,$2-3,$3-9}' $dirout/dLUBS_bestscore_info.bed > $dirout/LFY_coord_on_dLUBS_bestscore.bed ## new
else
	awk -v OFS="\t" '$4=="+"{print $1,$2+10,$3+3}$4=="-"{print $1,$2-3,$3-10}' $dirout/dLUBS_bestscore_info.bed > $dirout/LFY_coord_on_dLUBS_bestscore.bed ## old
fi


## retrieve FASTA sequences and calculate LFY scores within dLUBS
if [[ ! -f $dirout/LFY_on_dLUBS.fas.scores ]]; then
	bedtools getfasta -fi $tair10 -fo $dirout/LFY_on_dLUBS.fas -bed $dirout/LFY_coord_on_dLUBS_bestscore.bed
	python $scores_prog -f $dirout/LFY_on_dLUBS.fas -m $LFY_matrice -o $dirout/
fi


## associate each bestscore with its 50bp peak coordinates
sed 's/-/:/' $dirout/LU_peak_dLUBS_bestscore.bed | awk -v FS="[:\t]" -v OFS="\t" '{print $1,$2,$3,$6}' | sort -k1,1 -k2,3n > $dirout/LU_peak_dLUBS_bestscore_positions.bed



## create file with 50bp peak coordinates + dLUBS and LFY scores
paste $dirout/LU_peak_dLUBS_bestscore_positions.bed <(sed 's/-/\t/' $dirout/LFY_on_dLUBS.fas.scores | sed 's/:/\t/g' | sort -k1,1 -k2,3n | awk '{print $10}' ) | awk -v OFS="\t" '{print $0}' > $dirout/dLUBS_LFY_scores_50bppeaks.tsv



## LFY scores calculation on central 50 bp
if [[ ! -f $dirout/LU_peak_LFY.fas.scores ]]; then
	bedtools getfasta -fi $tair10 -fo $dirout/LU_peak_LFY.fas -bed $dirout/LU_peak.bed ## use same bed coordinates used for dLUBS
	python $scores_prog -f $dirout/LU_peak_LFY.fas -m $LFY_matrice -o $dirout/
fi



## select best LFY score per peak and cut around LFY matrix
if [[ ! -f $dirout/LFY_bestscore_info.bed ]]; then 
	sort -k1,1 -k8,8nr $dirout/LU_peak_LFY.fas.scores \
			| awk -v OFS="\t" '{print $2,$4,$8,$1}' \
			| uniq -f3 \
			| awk -v OFS="\t" '{print $4,$1,$2,$3}' > $dirout/LU_peak_LFY_bestscore.bed
	
	sort -u -k1,1 $dirout/LU_peak_LFY_bestscore.bed | sed 's/-/\t/' | sed 's/:/\t/g' | awk -v FS="\t" -v OFS="\t" '{print $1,$2+$4-1,$2+$4+19-1,$5}' | sort -k1,1 -k2,3n > $dirout/LFY_bestscore_info.bed
fi



# 	rm $dirout/table_all_scores.tsv
if [[ ! -f $dirout/table_all_scores.tsv ]]; then
	echo -e "chrom\tstart\tstop\tfile\tLFYamp\tLFY-UFO\tCFC\tdLUBS\tLFY_on_dLUBS\tLFY" > $dirout/table_all_scores.tsv

	paste $dirout/table_sorted.bed <( sort -k1,1 -k2,3n $dirout/dLUBS_LFY_scores_50bppeaks.tsv | awk '{print $4,$5}') <(sed 's/-/\t/' $dirout/LU_peak_LFY_bestscore.bed | sed 's/:/\t/g' | sort -k1,1 -k2,3n | awk -v OFS="\t" '{print $6}' ) | awk -v OFS="\t" '{print $0}' >> $dirout/table_all_scores.tsv
fi



## produce paper fig plot from Rstudio:
if [[ ! -f $dirout/boxplots_dLUBS_vs_LFY.png ]]; then
	conda run -n LFYUFO_figs Rscript scores_dLUBS_LFY_clean.R $dirout
fi


## create files with coordinates of best site for dLUBS, LFY in dLUBS and LFY for IGB
## dLUBS bestscore
sort -u -k1,1 $dirout/LU_peak_dLUBS_bestscore.bed | sed 's/-/\t/' | sed 's/:/\t/g' | awk -v FS="\t" -v OFS="\t" '{print $1,$2+$4-1,$2+$4+25-1,$6}' | sort -k1,1 -k2,3n > $dirout/dLUBS_bestscore_coord_IGB.bed


## LFY in dLUBS
awk -v OFS="\t" '{print $1,$8}' $dirout/LFY_on_dLUBS.fas.scores | sed 's/:/\t/g' | sed 's/-/\t/' | awk -v OFS="\t" '{print $0}' > $dirout/LFY_on_dLUBS_bestscore_coord_IGB.bed


## LFY bestscore
sort -u -k1,1 $dirout/LU_peak_LFY_bestscore.bed | sed 's/-/\t/' | sed 's/:/\t/g' | awk -v FS="\t" -v OFS="\t" '{print $1,$2+$4-1,$2+$4+19-1,$6}' | sort -k1,1 -k2,3n > $dirout/LFY_bestscore_coord_IGB.bed


exit 0

