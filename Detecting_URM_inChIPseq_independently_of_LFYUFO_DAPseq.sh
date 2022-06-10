#version="DAP_global_analysis_p3.7"
version="TFgenomicsAnalysis"
PATH_TO_COMPIL=/home/312.6-Flo_Re/312.6.1-Commun/scripts/$version
source $PATH_TO_COMPIL/compil_functions.sh
source $PATH_TO_COMPIL/compil_usages.sh
scores_prog=$PATH_TO_COMPIL/bin/scores.py
genome=/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas
annotation=/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.bed
meme_prog=/home/prog/meme/meme_4.12.0/bin/meme-chip # v4.12.0
meme2meme=/home/prog/meme/meme_4.12.0/bin/meme2meme # v4.12.0
meme2pfm=$PATH_TO_COMPIL/bin/meme2pfm.sh

PFM=/home/312.6-Flo_Re/312.6.1-Commun/data/LFY.pfm # this is the LFY pfm with position interdependancies from Moiroud et al 2011 (PMID: 21515819)

#let's start from the original ChIP peaks
Peaks_dir_ChIP=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/results/Wellmer_q01/PeakCalling/Wellmer
Peaks_dir_DAP=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/results/DAP_q00001/PeakCalling/FL

main_dir=/home/312.6-Flo_Re/312.6.1-Commun/LFY/ChIPnonDAP_motifs/resuts
Peaks_dir=$main_dir/Peakcalling
work_dir=$main_dir/ChIPoriWellmer_may2022


#bedgraphs
covout=$work_dir/coverage_onChIP
mkdir -p $covout
LFY_ChIP=$Peaks_dir_ChIP/Wellmer_cpm.bdg
LFY_dap=$Peaks_dir_DAP/FL_cov.bdg
LFYUFO_dap=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY-UFO/results/PeakCalling/LFYUFO_cf/LFYUFO_cf_cov.bdg

ChIPpeaks=$Peaks_dir_ChIP/Wellmer_narrow.bed
cat $ChIPpeaks | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $1,$2,$3,"ChIP"}' | sed '1 i\chr\tbegin\tend\tname'  > $covout/tmp_peaks

bdg_array=("$LFY_ChIP" "$LFY_dap")
names_array=("ChIP" "LFYdap")
if [ ! -f $covout/table.tsv ]; then add_coverage -t $covout/tmp_peaks -b bdg_array[@] -n names_array[@] -od $covout; fi

cat $covout/table.tsv |  tail -n +2 | awk -v OFS='\t' '{print $1,$2,$3,$5/$6}' | sort -k4,4rn > $covout/table_chiptodap.tsv

cat $covout/table_chiptodap.tsv | awk '$4>2' | cut -f -3 > $covout/ChIPori_sup2ChIPtoLFYDAP.bed


# now we have a list of 298 peaks 
# those peaks have a relative binding affinity twice greater in ChIP relative to DAP 
# let's look for motif...
array_names=("ChIPori_sup2ChIPtoLFYDAP")

for name in "${array_names[@]}"
do
dirout=$work_dir/$name
mkdir -p $dirout
bedtools getfasta -fi $genome -bed  $covout/${name}.bed > $dirout/${name}.fasta

#get peaks summit:
size=55
dir25=$dirout/summit${size}bp
if [ ! -f $dir25/summit.fasta ]; then
mkdir -p $dir25
bedtools intersect -a $covout/${name}.bed -b $Peaks_dir_ChIP/Wellmer_cpm.bdg -wa -loj > $dir25/tmp.bed
cat $dir25/tmp.bed | awk -v OFS="" '{print $1,"_",$2,"_",$3"\t",$5,"\t",$6,"\t",$7}' | sort -k1,1 -k4,4rn |  sort -uk1,1 | awk -v s=$size '{a=int(($2+$3)/2); $5=a-s; $6=a+s; print}' | tr '_' '\t' | tr " " "\t" | cut -f 1,7,8 | awk '!seen[$0]++' > $dir25/summit.bed
bedtools getfasta -fi $genome -bed $dir25/summit.bed > $dir25/summit.fasta
fi


# scan LFY canonical BS (PWM with dependencies, Moiroud et al 2011) in the resized peaks
SC=-23 # PWM score cutoff as derived in Sayou et al 2016 Nat Comm
mkdir -p $dir25/scores_LFYdep
res_scores=$dir25/scores_LFYdep/summit.fasta.scores
best_scores=$dir25/scores_LFYdep/summit.fasta.scores.best
if [ ! -f $res_scores ]; then $Python_TFFM $scores_prog -m $PFM -f $dir25/summit.fasta -o $dir25/scores_LFYdep; fi
#cat $res_scores | sort -k1,1 -k8,8rn | sort -uk1,1 > $best_scores #we keep the best scoring binding site per resized peaks
cat $res_scores | sort -k1,1 -k8,8rn | awk -v Sco=$SC '$8>Sco'> $best_scores # we keep all binding sites with score above threshold per peak


# Loop through matrix size (i) and distance to TFBS (j) to reconstruct motif with MEME
echo "LET S LOOP!"
for i in {4..10}; do 
	for j in {1..20}; do
		echo $i
		echo $j
		
		target="size"${i}"_pos"${j}
		dirTarget=$dir25/LOOP_${SC}/${target}
		mkdir -p $dirTarget
		cat $best_scores | tr ':' '\t' | tr '-' '\t' |awk -v OFS='\t' -v I=$i -v J=$j '{print $1,$2+$4-(J+I),$2+$4-J}' > $dirTarget/left.bed
		cat $best_scores | tr ':' '\t' | tr '-' '\t' |awk -v OFS='\t' -v I=$i -v J=$j '{print $1,$2+$5+(J-2),$2+$5+(I+J-2)}' > $dirTarget/right.bed
		bedtools getfasta -fi $genome -bed $dirTarget/left.bed > $dirTarget/left.fasta
		bedtools getfasta -fi $genome -bed $dirTarget/right.bed > $dirTarget/right.fasta
		# rev comp (works only if each seq fits in one line)
		cat $dirTarget/right.fasta | while read L; do  echo $L; read L; echo "$L" | rev | tr "ATGC" "TACG" ; done > $dirTarget/right_rc.fasta
		cat $dirTarget/left.fasta $dirTarget/right_rc.fasta > $dirTarget/left_and_right.fasta
		
		SEQ=$dirTarget/left_and_right.fasta
		nSEQ=$(grep -c ">" $SEQ)
		mi=$i
		ma=$i
		if [ ! -f $dirTarget/meme_${mi}_${ma}/meme-chip.html ]
		then 
		mkdir -p $dirTarget/meme_${mi}_${ma}
		$meme_prog -oc $dirTarget/meme_${mi}_${ma} -nmeme $nSEQ -meme-maxsize $(calc $nSEQ*1000) -meme-minw $mi -meme-maxw $ma -meme-nmotifs 5 -dreme-m 0 -noecho $SEQ -seed 165 -db /home/312.6-Flo_Re/312.6.1-Commun/data/meme_db/motif_databases/JASPAR/JASPAR2018_CORE_plants_non-redundant.meme
		fi		
		
		done
		
done
done # end  the for name in "${array_names[@]}" loop



 