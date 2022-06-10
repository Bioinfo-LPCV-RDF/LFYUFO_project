source /home/312.6-Flo_Re/312.6.1-Commun/scripts/TFgenomicsAnalysis/compil_functions.sh
source /home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY-UFO/scripts/LFYUFO_github_scripts/get_bdg_from_peaks.sh

main_dir=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY-UFO/results
echo 'hello'


fastq_dir=$main_dir/Fastq
mapping_dir=$main_dir/Mapping
dir_peakcalling=$main_dir/PeakCalling
out_motif=$main_dir/motifs
mkdir -p $main_dir/comparisons
dir_comparisons=$main_dir/comparisons

data=/home/312.6-Flo_Re/312.6.1-Commun/data
genome=$data/tair10.fas
Arabidopsis_genome_bed=$data/tair10_whole_chrom.bed
gff3=$data/A_thaliana_phytozome_v12/Phytozome/PhytozomeV12/Athaliana/annotation/Athaliana_167_TAIR10.gene.gff3
annotation=$data/tair10.bed

scores_prog=/home/312.6-Flo_Re/312.6.1-Commun/scripts/TFgenomicsAnalysis/bin/scores.py


############# Mapping
seed=2785
threads=20


names_mapping=(	"LFYUFO_cfa" "LFYUFO_cfb" "LFYUFO_cfc" 
				"LFYm_a" "LFYm_b"
				"LFYmUFO_cfa" "LFYmUFO_cfb")
IDs=("1" "2" "3" "10" "11" "12" "13")

i=0


for name in "${names_mapping[@]}"; do
	if [[ ! -f ${i}.minimal.stats ]]; then
		Fastq1=/home/312.3-StrucDev/312.3.1-Commun/DAP_LFYraw/${IDs[$i]}_R1_001.fastq.gz
		Fastq2=/home/312.3-StrucDev/312.3.1-Commun/DAP_LFYraw/${IDs[$i]}_R2_001.fastq.gz
		
#  		main_mapping_Fastq -f1 $Fastq1 -f2 $Fastq2 -fd $fastq_dir -md $mapping_dir -n $name -s $seed -pr $threads
		i=$((i+1))
	fi
done

if [[ ! -f $mapping_dir/md5sum_summary.txt ]]; then
	for n in "${IDs[@]}"; do
		fastq1=/home/312.3-StrucDev/312.3.1-Commun/DAP_LFYraw/${n}_R1_001.fastq.gz
		echo $fastq1
		md5sum $fastq1 >> $mapping_dir/md5sum_summary.txt
		
		fastq2=/home/312.3-StrucDev/312.3.1-Commun/DAP_LFYraw/${n}_R2_001.fastq.gz
		md5sum $fastq2 >> $mapping_dir/md5sum_summary.txt
		echo $fastq2
	done
	
fi


############# PeakCalling
Genome_length=120000000
additional_args="-q 0.0001 --call-summits"

names_peakcalling=("LFYUFO_cf" "LFYm" "LFYmUFO_cf")

for name in "${names_peakcalling[@]}"; do
	if [[ ! -d ${dir_peakcalling}/${name} ]]; then
		echo "proceeding with peakcalling for ${name}"
		
		if [[ -d "$mapping_dir/${name}a" ]]; then
			if [[ ! -d "$mapping_dir/${name}c" ]]; then
				bam_dir=("$mapping_dir/${name}a" "$mapping_dir/${name}b")
# 				echo ${bam_dir[@]}
			else 
				bam_dir=("$mapping_dir/${name}a" "$mapping_dir/${name}b" "$mapping_dir/${name}c")
# 				echo ${bam_dir[@]}
			fi
		elif [[ -d "$mapping_dir/${name}_a" ]]; then
			bam_dir=("$mapping_dir/${name}_a" "$mapping_dir/${name}_b")
# 			echo ${bam_dir[@]}
		fi
		
		input=("$mapping_dir/control")
#  		main_peakcalling -id bam_dir[@] -cd input -od $dir_peakcalling -nc $name -g $Genome_length -add "$additional_args" -ps 400
		echo -e "peakcalling done\n"
	else
		echo -e "peakcalling was successful; proceeding to motif search\n"
		echo -e "===========================================\n"
	fi
done



## rep/rep plots
bash rep_rep_plots.sh
# plot them on Rstudio with plot_pretty_reprep.r



############# Comparisons
name1=LFYamp
LFY_peakcalling_dir=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY-pioneer/PeakCalling # directory containing LFY ampDAP peakcalling
if [[ ! -d $dir_peakcalling/${name1}3 ]]; then
	echo 'copying files...'
# 	mkdir -p $dir_peakcalling/$name1
# 	mkdir -p $dir_peakcalling/${name1}1
# 	mkdir -p $dir_peakcalling/${name1}2
# 	mkdir -p $dir_peakcalling/${name1}3
# 	cp $LFY_peakcalling_dir/$name1/* $dir_peakcalling/$name1
# 	cp $LFY_peakcalling_dir/${name1}1/* $dir_peakcalling/${name1}1
# 	cp $LFY_peakcalling_dir/${name1}2/* $dir_peakcalling/${name1}2
# 	cp $LFY_peakcalling_dir/${name1}3/* $dir_peakcalling/${name1}3
	echo 'done!'
fi


if [[ ! -f $dir_peakcalling/$name1/${name1}_narrow.bed ]]; then ## if LFY DAP data is not present in the current PeakCalling directory, then make a copy of it
	echo "copying file in right directory"
	mkdir -p $dir_peakcalling/$name1
	
	cp $LFY_peakcalling_dir/${name1}_narrow.bed $dir_peakcalling/$name1/${name1}_narrow.bed
	cp $LFY_peakcalling_dir/${name1}_cov.bdg $dir_peakcalling/$name1/${name1}_cov.bdg
fi


name2=LFYUFO_cf

reps1=("${name1}1" "${name1}2" "${name1}3")
reps2=("${name2}a" "${name2}b" "${name2}c")

## First comparison: LFY ampDAP vs LFY+UFO ampDAP
if [[ ! -f $dir_comparisons/${name1}_${name2}/table_${name1}_${name2}.csv ]]; then
	initial_comparison -n1 ${name1} -n2 ${name2} -od $dir_comparisons -id $dir_peakcalling -bd $mapping_dir -rep1 reps1[@] -rep2 reps2[@]
fi

tsvtable=$dir_comparisons/${name1}_${name2}/table_${name1}_${name2}_RiL_RiP.tsv ## table produced by initial_comparison and containing the union list of all peaks in the two input files, with coverage corresponding to each peak in each file and normalized by Reads in Library or Reads in Peaks



## Produce a tsv table (starting from the csv table obtained with initial_comparison) containing cov LFY | cov LFY+UFO | calculated CFR (LFY+UFO/LFY) and sorted based on the CFR value (7th column, from greatest to smallest)
if [[ ! -f $dir_comparisons/${name1}_${name2}/table_${name1}_${name2}_RIP_CFC.tsv ]]; then
	awk -F "\t" -v OFS="\t" 'NR >1 {print $1,$2,$3,$4,$7,$8,$8/$7}' $tsvtable | sort -k7,7nr > $dir_comparisons/${name1}_${name2}/table_${name1}_${name2}_RIP_CFC.tsv # Produce tsv table containing CFC value in 7th column, sorted by smallest to biggest CFC
fi

tsvtable=$dir_comparisons/${name1}_${name2}/table_${name1}_${name2}_RIP_CFC.tsv




### Produce bdg files for LFYamp and LFYUFO
name1="LFYamp"
name2="LFYUFO_cf"
if [[ ! -d $dir_comparisons/${name1}_${name2}/inPeaks/${name2} ]]; then 
	peaks=$dir_comparisons/${name1}_${name2}/table_${name1}_${name2}.csv
	output=$dir_comparisons/${name1}_${name2}
	
	## LFYamp
	sample_names=("LFYamp1" "LFYamp2" "LFYamp3")
	get_RIP_RIL_bdgs -n ${name1} -p $peaks -bam $mapping_dir -sn sample_names[@] -od $output -pd $dir_peakcalling -m "inPeaks" -s 123
	## remove unnecessary bdg files
	rm $output/inPeaks/${name1}/${sample_names[0]}_RIP_cov.bdg
	rm $output/inPeaks/${name1}/${sample_names[1]}_RIP_cov.bdg
	rm $output/inPeaks/${name1}/${sample_names[2]}_RIP_cov.bdg
	
	
	## LFYUFO_cf
	sample_names=("LFYUFO_cfa" "LFYUFO_cfb" "LFYUFO_cfc")
	get_RIP_RIL_bdgs -n ${name2} -p $peaks -bam $mapping_dir -sn sample_names[@] -od $output -pd $dir_peakcalling -m "inPeaks" -s 123
	## remove unnecessary bdg files
	rm $output/inPeaks/${name2}/${sample_names[0]}_RIP_cov.bdg
	rm $output/inPeaks/${name2}/${sample_names[1]}_RIP_cov.bdg
	rm $output/inPeaks/${name2}/${sample_names[2]}_RIP_cov.bdg
fi











############# Motif generation
## LFY-spe motif:
if [[ ! -f $out_motif/${name1}/${name1}.pfm ]]; then
	echo -e "\nCompute motif for LFY-specific peaks\n"
	maxmean_file=$dir_peakcalling/${name1}/${name1}_maxMean.bed
	tail -n 600 $maxmean_file | awk -v OFS="\t" '{print $1,$2,$3}' > $dir_peakcalling/${name1}/${name1}_spePeaks.bed
	bottom_600_peaks=$dir_peakcalling/${name1}/${name1}_spePeaks.bed ## use this peak file to compute motif 
	head $bottom_600_peaks
	echo "LFY preparation..."
	
	name=LFYamp
	compute_motif -p $bottom_600_peaks -n $name -g $genome -od $out_motif -ls 600 -nm 1 -mim 19 -mam 19 -s 1234 -pal # compute motif for LFY+UFO-specific peaks
fi



## LFY+UFO-spe:
if [[ ! -f $out_motif/mLUBS/${name2}/${name2}.pfm ]]; then
	## Select 600 peaks that are most specific for LFY+UFO (highest CFC)
	head -n 600 $tsvtable | awk -v OFS="\t" '{print $1,$2,$3}' > $dir_comparisons/${name1}_${name2}/${name2}_spePeaks_RIP.bed
	best_600_peaks=$dir_comparisons/${name1}_${name2}/${name2}_spePeaks_RIP.bed ## use this peak file to compute motif 
	head $best_600_peaks
	echo "dLUBS preparation..."
	
	## dLUBS: min size 20, max 30
	add_motif="-ls 600 -nm 1 -mim 20 -mam 30 -s 1234" 
	name=LFYUFO_cf
	compute_motif -p $best_600_peaks -n $name -g $genome -od $out_motif/dLUBS -ls 600 -nm 1 -mim 25 -mam 25 -s 1234 # compute motif for LFY+UFO-specific peaks
	
	
	echo "mLUBS preparation..."
	## mLUBS: min size 16, max 19
	add_motif="-ls 600 -nm 1 -mim 16 -mam 19 -s 1234" 
	name=LFYUFO_cf
	mkdir -p $out_motif/mLUBS_0802 ## create mLUBS folder for results
	compute_motif -p $best_600_peaks -n $name -g $genome -od $out_motif/mLUBS -ls 600 -nm 1 -mim 16 -mam 19 -s 1234 # compute motif for LFY+UFO-specific peaks
	
else
	echo "dLUBS and mLUBS matrices have already been produced, moving to NSG"
fi


############# NSG and ROCS
############# NSG and ROCs on top 20% CFC, excluding those already used to generate the motifs:
if [[ ! -f ${main_dir}/NSG_ts/${name2}_top20CFC_testingset.bed ]]; then
	echo "NGS generation"
	mkdir -p ${main_dir}/NSG_ts
	## calculate top 20% of the peak total
	top20_CFC=$(calc $(wc -l <$tsvtable)*0.2 | awk '{print int($1+0.5)}')
	echo $top20_CFC
	
	## Exclude peaks used for motif generation from set -> testing set
	head -n $top20_CFC $tsvtable | tail -n+601 > ${main_dir}/NSG_ts/${name2}_top20CFC_testingset_RIP.bed
	top20_CFC_ts=${main_dir}/NSG_ts/${name2}_top20CFC_testingset_RIP.bed ## produce testing set by taking only peaks that were not used to create motifs
	
	
	echo 'computing NS...'
	name=${name2}_top20CFC
	if [[ ! -f ${main_dir}/NSG_ts/${name}_1_neg.bed ]]; then
		compute_NS -p $top20_CFC_ts -n $name -g $genome -od ${main_dir}/NSG_ts -nb 1 -anf $annotation -s 1234
	else
		echo -e "\nNegative testing set already generated"
	fi
	echo 'here'
fi




############# Compute ROC curves & AUC
### dLUBS matrix vs mLUBS vs LFY (SYM-T) vs UFO trace on best LFY-UFO-spe peaks:
peakname=${name2}_top20CFC
peaks=("${main_dir}/NSG_ts/${peakname}_pos.bed")

matrix1=dLUBS; matrix2=mLUBS; matrix3=LFY;


names=("$matrix1" "$matrix2" "$matrix3")

neg_set=("${main_dir}/NSG_ts/${peakname}_1_neg.bed")


matrices=("$out_motif/$matrix1/${name2}/${name2}.pfm"
		  "$out_motif/$matrix2/${name2}/${name2}.pfm"
		  "/home/312.6-Flo_Re/312.6.1-Commun/data/${matrix3}.pfm")


colors=('#006B00' '#9ACD32' '#007FFF') 

name=${peakname}
echo "compute ROC"
if [[ ! -f ${main_dir}/ROCS/${name}/ROC.png ]]; then
	mkdir -p ${main_dir}/ROCS/${name}
	compute_ROCS -p peaks[@] -ns neg_set[@] -m matrices[@] -n names[@] -g $genome -od ${main_dir}/ROCS/${name} -color colors[@]
fi





## produce paper figure 3 with covplot, motifs and ROC curve on Rstudio
Rscript ./fig_3_Rstudio.r




## assess LFY site quality in dLUBS sites:
bash Assess_LFYqual.sh
## plot LFY quality supp figure with Rstudio and script (will be run from within Assess_LFYqual.sh directly)
Rscript scores_dLUBS_LFY_Rstudio.R




###### ChIP vs ampDAP comparison
bash microarray_ufo.sh # to analyse ufo microarray data
bash LFYChIP_CFC2_comparisons.sh # to merge all ChIPs, compare to ampDAP LFYUFO and cross with results from microarray analysis



### LFY K249R mutation
## comparison LFYm vs LFYmUFO
name1="LFYm"
name2="LFYmUFO_cf"
reps1=("${name1}_a" "${name1}_b")
reps2=("${name2}a" "${name2}b")

if [[ ! -f $dir_comparisons/${name1}_${name2}/table_${name1}_${name2}.csv ]]; then
	echo "comparison ${name1} vs ${name2}"
	initial_comparison -n1 ${name1} -n2 ${name2} -od $dir_comparisons -id $dir_peakcalling -bd $mapping_dir -rep1 reps1[@] -rep2 reps2[@]
fi


## comparison LFYamp vs LFYm
name1="LFYamp"
name2="LFYm"
reps1=("${name1}1" "${name1}2" "${name1}3")
reps2=("${name2}_a" "${name2}_b")

if [[ ! -f $dir_comparisons/${name1}_${name2}/table_${name1}_${name2}.csv ]]; then
	echo "comparison ${name1} vs ${name2}"
	mkdir -p $dir_comparisons/${name1}_${name2}
# 	exit 0
	initial_comparison -n1 ${name1} -n2 ${name2} -od $dir_comparisons -id $dir_peakcalling -bd $mapping_dir -rep1 reps1[@] -rep2 reps2[@]
fi






##### calculate dLUBS scores genome-wide
## retrieve whole genome sequence
# rm -r $out_motif/dLUBS_0802/scores
mkdir -p $out_motif/dLUBS/scores
out_scores=$out_motif/dLUBS/scores

if [[ ! -f $out_scores/scores_dLUBS_th15.bdg ]]; then
	echo 'here'
	name2="LFYUFO_cf"
	dLUBS=$out_motif/dLUBS_0802/${name2}/${name2}.pfm
	
	th=15
	echo $th; echo -$th
	matrixlen=$(wc -l $dLUBS | awk '{print $1-2}')
	
	echo 'calculating dLUBS genome-wide scores'
	python $scores_prog -m $dLUBS -f $genome -o $out_scores -th -${th}
	echo 'scores done!'
	
	## create bdg for IGB
	awk -v OFS="\t" -v matrixlen=$matrixlen -v th=$th '{print $1,$2,$2+matrixlen,$8+th}' $out_scores/tair10.fas.scores > $out_scores/scores_dLUBS_th15.bdg
	
	
	## create bed
	awk -v OFS="\t" -v matrixlen=$matrixlen -v th=$th '$3==1{print $1,$2,$2+matrixlen,"peak_"$1"_"$2 ,$8+th,"+";next}{print $1,$2,$2+matrixlen,"peak_"$1"_"$2 ,$8+th,"-"}' $out_scores/tair10.fas.scores > $out_scores/scores_dLUBS_th15.bed
fi




##### start from LFYamp-LFYUFO comparison peaks and calculate RIP
##### for LFYm and LFYmUFO
echo "use LFYUFO/LFY peaks for LFYm and LFYmUFO"
awk -v OFS="\t" '{print $1,$2,$3}' $tsvtable | sort -k1,1 -k2,3n > $dir_comparisons/LFYamp_LFYUFO_LFYm_LFYmUFO/LFYamp_LFYUFO_cf_peaks.bed


## retrieve RIP normalization
mode="inPeaks"
sample_names=(	"LFYamp1" "LFYamp2" "LFYamp3"
				"LFYUFO_cfa" "LFYUFO_cfb" "LFYUFO_cfc"
				"LFYm_a" "LFYm_b"
				"LFYmUFO_cfa" "LFYmUFO_cfb")
peaks=$dir_comparisons/LFYamp_LFYUFO_LFYm_LFYmUFO/LFYamp_LFYUFO_cf_peaks.bed

compute_rpkmrip_rpkmril -p $peaks -bd $main_dir/Mapping -pd $main_dir/PeakCalling -sn sample_names[@] -m $mode -o $dir_comparisons/LFYamp_LFYUFO_LFYm_LFYmUFO


## get bdgs for mutant LFY experiments
if [[ ! -d $dir_comparisons/LFYamp_LFYUFO_LFYm_LFYmUFO/inPeaks/LFYmUFO ]]; then 
	sample_names3=("LFYm_a" "LFYm_b")
	sample_names4=("LFYmUFO_cfa" "LFYmUFO_cfb")
	
	name3="LFYm"
	name4="LFYmUFO"
	
	peaks=$dir_comparisons/LFYamp_LFYUFO_LFYm_LFYmUFO/LFYamp_LFYUFO_cf_peaks.bed
	output=$dir_comparisons/${name1}_${name2}_${name3}_${name4}
	
	## LFYm
	get_RIP_RIL_bdgs -n ${name3} -p $peaks -bam $mapping_dir -sn sample_names3[@] -od $output -pd $dir_peakcalling -m "inPeaks" -s 123
	rm $output/inPeaks/${name3}/${sample_names3[0]}_RIP_cov.bdg
	rm $output/inPeaks/${name3}/${sample_names3[1]}_RIP_cov.bdg
	
	## LFYmUFO
	get_RIP_RIL_bdgs -n ${name4} -p $peaks -bam $mapping_dir -sn sample_names4[@] -od $output -pd $dir_peakcalling -m "inPeaks" -s 123
	rm $output/inPeaks/${name4}/${sample_names4[0]}_RIP_cov.bdg
	rm $output/inPeaks/${name4}/${sample_names4[1]}_RIP_cov.bdg
fi



if [[ ! -f $dir_peakcalling/LFYmUFO_cf/LFYK249RUFO_narrow.bed ]]; then
	submitted_data=$main_dir/GEO_submitted_data
	cp $dir_comparisons/LFYamp_LFYUFO_LFYm_LFYmUFO/inPeaks/LFYmUFO/LFYmUFO_RIP_cov.bdg $submitted_data/LFYK249RUFO_RIP_cov.bdg
	cp $dir_comparisons/LFYamp_LFYUFO_LFYm_LFYmUFO/inPeaks/LFYm/LFYm_RIP_cov.bdg $submitted_data/LFYK249R_RIP_cov.bdg
	cp $dir_comparisons/LFYamp_LFYUFO_cf/inPeaks/LFYUFO_cf/LFYUFO_cf_RIP_cov.bdg $submitted_data/LFYUFO_RIP_cov.bdg
	cp $dir_peakcalling/LFYUFO_cf/LFYUFO_cf_narrow.bed $submitted_data/LFYUFO_narrow.bed
	cp $dir_peakcalling/LFYm/LFYm_narrow.bed $submitted_data/LFYK249R_narrow.bed
	cp $dir_peakcalling/LFYmUFO_cf/LFYmUFO_cf_narrow.bed $submitted_data/LFYK249RUFO_narrow.bed
	echo 'done!'
fi




## make all LFYK249R-related figures
if [[ ! -f $dir_comparisons/LFYamp_LFYUFO_LFYm_LFYmUFO/violinplot_top20perc_CFC_pub.png ]]; then
	conda run -n LFYUFO_figs Rscript LFY_mutant_plots.r "$dir_comparisons/LFYm_LFYmUFO_cf" "$dir_comparisons/LFYamp_LFYm" "$dir_comparisons/LFYamp_LFYUFO_LFYm_LFYmUFO"
fi


### retrieve RIL normalization for the same set
# echo 'retrieve RIL'
compute_rpkmrip_rpkmril -p $peaks -bd $main_dir/Mapping -pd $main_dir/PeakCalling -sn sample_names[@] -m "inLibs" -o $dir_comparisons/LFYamp_LFYUFO_LFYm_LFYmUFO/RIL


##### calculate dLUBS and LFY scores on AP3 promoter
## retrieve AP3 promoter sequence (~1 kb upstream)
name1=LFYamp
name2=LFYUFO_cf

out_AP3=$dir_comparisons/${name1}_${name2}/AP3
mkdir -p $out_AP3
if [[ ! -f $out_AP3/AP3_mLUBS_light.scores ]]; then 
	echo 'AP3 scores calculation'
	grep "AT3G54340" $gff3 | grep "five_prime_UTR" | awk -v OFS="\t" '{print $1,$4-1,$4+929}' | sed 's/Chr/chr/g' > $out_AP3/AP3_prom.bed ## print coordinates upstream of ATG
	
	bedtools getfasta -fi $genome -bed $out_AP3/AP3_prom.bed > $out_AP3/AP3_prom.fa
	
	# dLUBS score calculation
	dLUBS=$out_motif/dLUBS/${name2}/${name2}.pfm
	python $scores_prog -m $dLUBS -f $out_AP3/AP3_prom.fa -o $out_AP3 ## dLUBS
	mv $out_AP3/AP3_prom.fa.scores $out_AP3/AP3_dLUBS.scores
	awk -v OFS="\t" '{print $2,$4,$8}' $out_AP3/AP3_dLUBS.scores > $out_AP3/AP3_dLUBS_light.scores  #+929
	
	
	# mLUBS score calculation
	mLUBS=$out_motif/mLUBS/${name2}/${name2}.pfm
	python $scores_prog -m $mLUBS -f $out_AP3/AP3_prom.fa -o $out_AP3 ## mLUBS
	mv $out_AP3/AP3_prom.fa.scores $out_AP3/AP3_mLUBS.scores
	awk -v OFS="\t" '{print $2,$4,$8}' $out_AP3/AP3_mLUBS.scores > $out_AP3/AP3_mLUBS_light.scores
	
	
	# LFY scores
	LFY_matrix=/home/312.6-Flo_Re/312.6.1-Commun/data/LFY.pfm
	python $scores_prog -m $LFY_matrix -f $out_AP3/AP3_prom.fa -o $out_AP3
	mv $out_AP3/AP3_prom.fa.scores $out_AP3/AP3_LFY.scores
	awk -v OFS="\t" '{print $2,$4,$8}' $out_AP3/AP3_LFY.scores > $out_AP3/AP3_LFY_light.scores
	
fi


## plot scores with RStudio
if [[ ! -f $out_AP3/AP3_prom_scores_LFY.png ]]; then
	conda run -n LFYUFO_figs Rscript plot_AP3_scores.r ${out_AP3} AP3
fi


## dLUBS and LFY scores on RBE promoter
mkdir -p $dir_comparisons/${name1}_${name2}/RBE
out_RBE=$dir_comparisons/${name1}_${name2}/RBE
if [[ ! -f $out_RBE/RBE_mLUBS_light.scores ]]; then 
	echo 'RBE scores calculation'
	grep "AT5G06070" $gff3 | grep "five_prime_UTR" | awk -v OFS="\t" '{print $1,$4-1,$4+1563}' | sed 's/Chr/chr/g' > $out_RBE/RBE_prom.bed
	bedtools getfasta -fi $genome -bed $out_RBE/RBE_prom.bed > $out_RBE/RBE_prom.fa
	
	# dLUBS score calculation
	dLUBS=$out_motif/dLUBS/${name2}/${name2}.pfm
	python $scores_prog -m $dLUBS -f $out_RBE/RBE_prom.fa -o $out_RBE
	mv $out_RBE/RBE_prom.fa.scores $out_RBE/RBE_dLUBS.scores
	awk -v OFS="\t" '{print $2,$4,$8}' $out_RBE/RBE_dLUBS.scores > $out_RBE/RBE_dLUBS_light.scores  #+929
	
	
	# mLUBS scores
	mLUBS=$out_motif/mLUBS/${name2}/${name2}.pfm
	python $scores_prog -m $mLUBS -f $out_RBE/RBE_prom.fa -o $out_RBE
	mv $out_RBE/RBE_prom.fa.scores $out_RBE/RBE_mLUBS.scores
	awk -v OFS="\t" '{print $2,$4,$8}' $out_RBE/RBE_mLUBS.scores > $out_RBE/RBE_mLUBS_light.scores
	
	
	# LFY scores
	LFY_matrix=/home/312.6-Flo_Re/312.6.1-Commun/data/LFY.pfm
	python $scores_prog -m $LFY_matrix -f $out_RBE/RBE_prom.fa -o $out_RBE
	mv $out_RBE/RBE_prom.fa.scores $out_RBE/RBE_LFY.scores
	awk -v OFS="\t" '{print $2,$4,$8}' $out_RBE/RBE_LFY.scores > $out_RBE/RBE_LFY_light.scores
	
	
fi


if [[ ! -f $out_RBE/RBE_prom_scores_LFY.png ]]; then
	conda run -n LFYUFO_figs Rscript plot_RBE_scores.r $out_RBE RBE
fi




##### calculate dLUBS and mLUBS scores on PI promoter
## retrieve PI promoter sequence (1 kb upstream)
echo 'PI scores calculation'
mkdir -p $dir_comparisons/${name1}_${name2}/PI
out_PI=$dir_comparisons/${name1}_${name2}/PI
if [[ ! -f $out_PI/PI_LFY_light.scores ]]; then 
	grep "AT5G20240" $gff3 | grep "five_prime_UTR" | awk -v OFS="\t" '{print $1,$5-1500,$5}' | sed 's/Chr/chr/g' > $out_PI/PI_prom.bed ## print coordinates 1,5kb upstream of ATG



	bedtools getfasta -fi $genome -bed $out_PI/PI_prom.bed > $out_PI/PI_prom.fa

	# score calculation
	dLUBS=$out_motif/dLUBS/${name2}/${name2}.pfm
	mLUBS=$out_motif/mLUBS/${name2}/${name2}.pfm
	LFY_matrix=/home/312.6-Flo_Re/312.6.1-Commun/data/${matrix3}.pfm

	python $scores_prog -m $dLUBS -f $out_PI/PI_prom.fa -o $out_PI ## dLUBS
	mv $out_PI/PI_prom.fa.scores $out_PI/PI_dLUBS.scores
	awk -v OFS="\t" '{print $2-1501,$4,$8}' $out_PI/PI_dLUBS.scores > $out_PI/PI_dLUBS_light.scores 
	
	python $scores_prog -m $mLUBS -f $out_PI/PI_prom.fa -o $out_PI ## mLUBS
	mv $out_PI/PI_prom.fa.scores $out_PI/PI_mLUBS.scores
	awk -v OFS="\t" '{print $2-1501,$4,$8}' $out_PI/PI_mLUBS.scores > $out_PI/PI_mLUBS_light.scores 


	python $scores_prog -m $LFY_matrix -f $out_PI/PI_prom.fa -o $out_PI ## LFY
	mv $out_PI/PI_prom.fa.scores $out_PI/PI_LFY.scores
	awk -v OFS="\t" '{print $2-1501,$4,$8}' $out_PI/PI_LFY.scores > $out_PI/PI_LFY_light.scores

	## plot scores with RStudio
	plot_PI_scores.r
fi

if [[ ! -f $out_PI/PI_prom_scores_LFY.png ]]; then
	conda run -n LFYUFO_figs Rscript plot_PI_scores.r $out_PI PI 
fi

exit 0

