#!/bin/bash
source /home/312.6-Flo_Re/312.6.1-Commun/scripts/TFgenomicsAnalysis/compil_functions.sh

LFY_targets=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets
LFY_UFO=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY-UFO

dir_comparisons=$LFY_UFO/results/comparisons
dir_peakcalling=$LFY_UFO/results/PeakCalling

names=("Sayou" "Moyroud" "Wellmer" "callus")
output=$LFY_UFO/results/comparisons/LFY_ChIP_comp

# comparison -n names[@] -id $LFY_targets/results/ChIP-seq/PeakCalling -od $output

## Produce corresponding .bed file:
# awk -F "\t" -v OFS="\t" '{print $1,$2,$3}' $output/table_peaks.csv > $output/ChIP_union_peaks.bed



## Hybrid method: CFC to select LFY-UFO-spe peaks compared to LFYamp, but bedtools to cross with LFY ChIP union peaks
## comparison LFY ChIP union vs LFYUFO: only keep common peaks
outdir=$dir_comparisons/ChIP_ampDAP_CFC_bt_2305
mkdir -p $outdir

percentage=20

## Intersect between LFY ChIP union and top 20% of CFC (LFYUFO/LFYamp)
LFYUFO_spe_peaks=$LFY_UFO/results/NSG_ts/LFYUFO_cf_top20CFC_pos.bed
ChIP_union_peaks=$dir_peakcalling/LFY_ChIP_union/LFY_ChIP_union_narrow.bed

bedtools intersect -a $ChIP_union_peaks -b $LFYUFO_spe_peaks -wa -f 0.8 -F 0.8 -e > $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt.bed
echo "Hybrid comparison --- peaks bound in LFY ChIP and LFYUFO ampDAP (top $percentage% CFC): $(wc -l <$outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt.bed)"



## comparison [ LFY ChIP union && LFYUFO ] vs LFYamp: exclude LFYamp
LFYamp_peaks=$dir_peakcalling/LFYamp/LFYamp_narrow.bed
bedtools intersect -a $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt.bed -b $LFYamp_peaks -v -wa -f 0.6 -F 0.6 -e > $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_nonLFYamp.bed



## Retrieve corresponding bound genes
## Retrieve genes corresponding to LFY ChIP_and_LFYUFO DAP peaks
echo "Retrieving genes..."
if [[ ! -f $LFY_UFO/data/A_thaliana_extendedgenes.bed ]]; then
	genes=/home/312.6-Flo_Re/312.6.1-Commun/data/tair10_prep.gff3
	awk -v FS="[\t=;]" -v OFS="\t" '$3=="gene" && $7 =="+" {print $1,$4-3000,$5+1000,$10;next} $3=="gene" && $7 =="-" {print $1,$4-1000,$5+3000,$10;next}' $genes | sed 's/C/c/' | awk -v FS="[\t\.]" -v OFS="\t" '$2<0{print $1,"0",$3,$4;next}{print $1,$2,$3,$4}' > $LFY_UFO/data/A_thaliana_extendedgenes.bed
fi

extendedAT_genes=$LFY_UFO/data/A_thaliana_extendedgenes.bed

## intersect gene regions with peaks
bedtools intersect -a $extendedAT_genes -b $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_nonLFYamp.bed -wa -wb -f 0.8 -F 0.8 -e | awk -v OFS="\t" '{print $4}' | sort > $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_bound_genes.bed
wc -l <$outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_bound_genes.bed

## same but keeping peaks coordinates
bedtools intersect -a $extendedAT_genes -b $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_nonLFYamp.bed -wa -wb -f 0.8 -F 0.8 -e | awk -v OFS="\t" '{print $5,$6,$7,$8,$4}' | sort -k1,1 -k2,3 | uniq -u > $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_peaks_genes.bed


echo "Hybrid comparison --- genes bound in LFY ChIP and LFYUFO ampDAP but not in LFYamp: $(wc -l <$outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_bound_genes.bed)"




## Compare with lists of DEGs
pval=0.05
th=("0.5")

DEG_result=/home/312.6-Flo_Re/312.6.1-Commun/data/UFO_microarrays/ufo_inflor/wt_ufo_inflor_volcanodata.csv


for ((i=0;i<${#th[@]};i++)); do
	if [[ ! -f /home/312.6-Flo_Re/312.6.1-Commun/data/UFO_microarrays/ufo_inflor/DEG_wt_ufo_inflor_FC${th[i]}_updown.txt ]]; then
		awk -F ";" -v pval=$pval -v th=${th[i]} '$4<pval && $3>=th {print $2} $4<pval && $3<=(-th) {print $2}' $DEG_result | sed 's/"//g' | grep -v "ATMG" | grep -v "ATCG" | awk 'NR>1{print $0}' > /home/312.6-Flo_Re/312.6.1-Commun/data/UFO_microarrays/ufo_inflor/DEG_wt_ufo_inflor_FC${th[i]}.txt
		
		awk -F ";" -v pval=$pval -v th=${th[i]} '$4<pval && $3>=th {print $2,"up"} $4<pval && $3<=(-th) {print $2,"down"}' $DEG_result | sed 's/"//g' | grep -v "ATMG" | grep -v "ATCG" | awk 'NR>1{print $0}'> /home/312.6-Flo_Re/312.6.1-Commun/data/UFO_microarrays/ufo_inflor/DEG_wt_ufo_inflor_FC${th[i]}_updown.txt

	fi

	DEGs=/home/312.6-Flo_Re/312.6.1-Commun/data/UFO_microarrays/ufo_inflor/DEG_wt_ufo_inflor_FC${th[i]}.txt
	
# echo "Total number of DEGs at FC +/- ${th[i]}: $(wc -l </home/312.6-Flo_Re/312.6.1-Commun/data/UFO_microarrays/ufo_inflor/DEG_wt_ufo_inflor_FC${th[i]}.txt)"
	
	
	

# 		if [[ ! -f $dir_comparisons/ChIP_ampDAP_comp/DEG_comp/${name1}_${name2}_DEG_FC${th[i]}.txt ]]; then
	echo "Crossing with DEGs..."
	
	
		## common genes between DEGs set and LFYUFO-spe bound genes
		LFYUFO_spe_genes=$outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_bound_genes.bed
		
		grep -f $DEGs $LFYUFO_spe_genes > $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_DEG_FC${th[i]}.txt
		
		grep -f $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_DEG_FC${th[i]}.txt /home/312.6-Flo_Re/312.6.1-Commun/data/UFO_microarrays/ufo_inflor/DEG_wt_ufo_inflor_FC${th[i]}_updown.txt | sort -k1,1 > $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_DEG_FC${th[i]}_updown.txt
		
		
		echo "genes bound in LFY ChIP (ChIP union), LFY-UFO DAP and also DEG in inflorescence experiments at FC +/- ${th[i]}: $(wc -l <$outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_DEG_FC${th[i]}.txt)"
		
		## Retrieve corresponding gene positions to scroll IGB
# 			grep -f $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_DEG_FC${th[i]}.txt $extendedAT_genes > $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_DEG_FC${th[i]}.bed
		
# 		
		
		### Get gene descriptions of these genes:
		bash $LFY_targets/scripts/get_gene_descriptions.sh -f $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_DEG_FC${th[i]}.txt -o $outdir -a /home/312.6-Flo_Re/312.6.1-Commun/data/A_thaliana_phytozome_v12/Phytozome/PhytozomeV12/Athaliana/annotation/annotation_AT_TAIR10.csv
	
# 		fi
	
	## cross this list of genes with DEG in lfy mutant inflo
	lfy_inflo_results=$LFY_targets/data/lfyinflor/limma/wtlfyinfl_volcanodata.csv
	awk -F";" -v OFS="\t" '($3>=0.5&&$4<0.05){print $2};($3<=(-0.5)&&$4<0.05){print $2}' $lfy_inflo_results | grep -v ATM | grep -v ATC | sort > $outdir/wt_lfy_infl_DEG05.txt
	
	comm -12 $outdir/ChIP_and_LFYUFO-spe_top${percentage}CFC_bt_DEG_FC${th[i]}.txt $outdir/wt_lfy_infl_DEG05.txt | sort > $outdir/ChIP_and_LFYUFO-spe_top20CFC_bt_DEG_FC0.5_ufo_lfy.txt
	



