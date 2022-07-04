# list passed as arg
# list=("a" "b")
# test(){ testVar=("${!1}"); echo "${testVar[@]}";  }
# test list[@]
IFS=$'\n\t'
###### Needed in path
# fastqc #v0.11.7 mapping_Fastq_bowtie2
# NGmerge #v0.2_dev mapping_Fastq_bowtie2
# bedtools #v2.27.1 peakcalling_MACS2
# mspc #v4.0.0 main_peakcalling
# inkscape #v0.92.2 compute_space
# python 2.7.9  #Python 3 version in progress
# R 3.5.0

PATH_TO_COMPIL=/home/312.6-Flo_Re/312.6.1-Commun/scripts/TFgenomicsAnalysis
source $PATH_TO_COMPIL/compil_usages.sh
PATH_TO_DNAshapedTFBS=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DNAshapedTFBS-master

# Python Envs (made with conda, need manual installation of GHMM and TFFM packages, both are stuck in python 2.7)
Python_TFFM=/home/312.6-Flo_Re/312.6.1-Commun/Programs/Anaconda3/envs/py27_tffm/bin/python
R_36=/home/312.6-Flo_Re/312.6.1-Commun/Programs/Anaconda3/envs/py38_R36/bin/Rscript
#
# PATH_to_PWMscan=/home/312.6-Flo_Re/312.6.1-Commun/lib/pwmscan
# download_SRA
SRAtoolkit=/home/312.3-StrucDev/312.3.1-Commun/SRAtoolkit/sratoolkit.3.0.0-ubuntu64/bin/
# mapping_Fastq_bowtie2
PATH_TO_BOWTIE_INDEX="/home/312.3-StrucDev/312.3.1-Commun/bowtie/bowtie2-2.3.4.1-linux-x86_64/indexes/at"
PATH_TO_BOWTIE2=$(dirname $(which bowtie2)) # v2.3.4.1
PATH_TO_SAMTOOLS=$(dirname $(which samtools)) #v1.8 (using htslib v1.8)
PATH_TO_BEDTOOLS2=$(dirname $(readlink -f $(which bedtools))) #v2.27.1
# peakcalling_MACS2
PATH_TO_MACS=$(dirname $(readlink -f $(which macs2))) # v2.1.1.20160309
mspc_mk=/home/312.6-Flo_Re/312.6.1-Commun/lib/MSPC/mspc/mspc
# initial_comparison
export TMPDIR2=/nobackup/ # use a disk with some space (~20Go for safety) # Jérémy j'ai ajouté "2" à "TMPDIR" car ça créé un conflit avec le TMPDIR du wrapper TF_benchmark... mais le mieux serait de définir cette variable en local si possible ???. Romain
do_plot_quant=$PATH_TO_COMPIL/bin/Hist_cov_gen.r
do_plot_quantRIP=$PATH_TO_COMPIL/bin/Hist_cov_genRIP_RIL.r
#do_plot_quant=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY-UFO/results/test_initial_comp_14012022/Hist_cov_gen_modifROmain.r

merge_peaks=$PATH_TO_COMPIL/bin/merge_all_peaks.py
compute_coverage=$PATH_TO_COMPIL/bin/compute_coverage.py
# comparison
merge_peaks_Nsets=$PATH_TO_COMPIL/bin/merge_all_peaks_Nsets.py
# compute_rpkmrip_rpkmril_WIP (this tow compare to conditions, need to be checked/finished...
norma_cov_plots=$PATH_TO_COMPIL/bin/rpkmrip_rpkmril_perRep_per_Condition_covPlots.R
# compute_rpkmrip_rpkmril
rpkmrip_rpkmril=$PATH_TO_COMPIL/bin/rpkmrip_rpkmril.R
# Compute motif
meme_prog=/home/prog/meme/meme_4.12.0/bin/meme-chip # v4.12.0
meme2meme=/home/prog/meme/meme_4.12.0/bin/meme2meme # v4.12.0
meme2pfm=$PATH_TO_COMPIL/bin/meme2pfm.sh
prepMEMEforPalTFFM=$PATH_TO_COMPIL/bin/prepMEMEforPalTFFM.py
pfmTOtffm=$PATH_TO_COMPIL/bin/get_tffm.py
transpose=$PATH_TO_COMPIL/bin/transpose_matrix.r
# Define where are located the DNA shape data from GBshape
araTha=/home/312.6-Flo_Re/312.6.1-Commun/data/DNAshape/A_thaliana
helt=$araTha/HeIT.bw;
mgw=$araTha/MGW.bw;
prot=$araTha/ProT.bw;
roll=$araTha/Roll.bw;
helt2=$araTha/HeIT2o.bw;
mgw2=$araTha/MGW2o.bw;
prot2=$araTha/ProT2o.bw;
roll2=$araTha/Roll2o.bw;
ComputeDNAshaped=$PATH_TO_DNAshapedTFBS/DNAshapedTFBS.py
HeatmapDNAshape=$PATH_TO_DNAshapedTFBS/feature_importance_heatmap.py
# compute kmer
GEM=/home/312.6-Flo_Re/312.6.1-Commun/Romain/gem/gem.jar
# prep_annotation
prepare_gff=$PATH_TO_COMPIL/bin/prepare_gff.py # DEPRECATED
generate_BedFromGff=$PATH_TO_COMPIL/bin/generate_bed_gff.py # DEPRECATED
bedops=/home/312.3-StrucDev/312.3.1-Commun/bedops/bin/bedops # v2.4.38 # DEPRECATED
bedmap=/home/312.3-StrucDev/312.3.1-Commun/bedops/bin/bedmap # v2.4.38 # DEPRECATED
# compute_NS
negative_set_script=$PATH_TO_COMPIL/bin/neg_set_gen.py
# compute_ROCS
parse_KSM_scores=$PATH_TO_COMPIL/bin/parse_KSM_scores.py
pocc_pfm=$PATH_TO_COMPIL/bin/compute_POcc.py
tffmscores=$PATH_TO_COMPIL/bin/get_best_score_tffm.py
scores_prog=$PATH_TO_COMPIL/bin/scores.py
plot_ROCS_prog=$PATH_TO_COMPIL/bin/plots_ROCS_multiple.py
F_score=$PATH_TO_COMPIL/bin/fscore.R
dimer_builder=$PATH_TO_COMPIL/bin/build_dimer_matrix.sh
# compute_space
spacing_mk=$PATH_TO_COMPIL/bin/get_interdistances.py
tffm_all_scores=/home/312.6-Flo_Re/312.6.1-Commun/scripts/get_all_score_tffm.py
Zscore_spacing=$PATH_TO_COMPIL/bin/Zacing.R
get_interdistances_2TF=$PATH_TO_COMPIL/bin/get_interdistances_2TF.py
Zacing_2TFs=$PATH_TO_COMPIL/bin/Zacing_2TF.R

# heatmap_reads
heatmap_mk=$PATH_TO_COMPIL/bin/show_reads.py
bdg_to_bw=$PATH_TO_COMPIL/bin/bedGraphToBigWig
impactspacing_mk=$PATH_TO_COMPIL/bin/compute_spacingScorev2.r
#methylation
methMap=/home/312.6-Flo_Re/312.6.1-Commun/data/Methylation_maps/Zhu_lab_PNAS_2016/GSM1876327_Col.mC.bed
full_methylation=$PATH_TO_COMPIL/bin/full_methylation_NEW.py
plot_meth_full=$PATH_TO_COMPIL/bin/plot_meth_full_NEW.R
figs_meth_violin=$PATH_TO_COMPIL/bin/figs_meth_violin_NEW.R

# Function to compute simple math with bash script
# usage : calc $M1/$M2 or calc $M1/$M2*M3
calc(){
awk "BEGIN { print "$*" }"; 
}

# used to create delimited string from array
join_by(){
local IFS="$1"; shift; echo "$*"; 
}

# join multiple files recursively
multijoin() {
	out=$1
	shift 1
	cat $1 | awk '{print $1}' > $out
	for f in $*; do join $out $f > tmp; mv tmp $out; done
}

# Function to write to log. Note that if "verbose" is set,
# the log messages will be on the std also
logme(){
    if [ ! -z $verbose ]
       then
           echo "${*}" >> $log 2>&1;
    else
           echo "${*}" >> $log
    fi
}


download_SRA(){ #download_SRA -f <LIST> -o <PATH> -t [INT]

while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-f)
			local list_file=("${!2}")
			echo "list of SRA ID set to: ${list_file[@]}";shift 2;;
        -o)
			local dir_out=$2
			echo "name of output directory set to: ${2}";shift 2;;
        -t)
			local threads=$2
			echo "number of threads set to: ${2}";shift 2;;
		-n)
			local name=("${!2}")
			echo "name set to: ${name[@]}";shift 2;;
		-h)
			usage ; return;;
		--help)
			usage ; return;;
		*)
			echo "Error in arguments"
			echo $1; usage download_SRA ; return;;
	esac
done

local Errors=0
if [ -z $list_file ]; then echo "ERROR: -f argument needed";Errors+=1;fi
if [ -z $dir_out ]; then echo "ERROR: -o argument needed";Errors+=1;fi
if [ -z $threads ]; then echo "-t argument not used, using 1 thread by default";local threads=1 ;fi
if [ $Errors -gt 0 ]; then usage download_SRA ; return 1; fi


mkdir -p -m 774 $dir_out
for ((i=0;i<${#list_file[@]};i++))
do
	if [ ! -s $dir_out/${list_file[i]}.fastq.gz ] && [ ! -f $dir_out/${list_file[i]}_R1.fastq.gz ] ; then
		$SRAtoolkit/prefetch `echo ${list_file[i]}` --max-size 30g -O $dir_out
		$SRAtoolkit/fasterq-dump `echo ${list_file[i]}` -O $dir_out -t $TMPDIR2 --split-3 -p -e $threads
		if [ -s $dir_out/${list_file[i]}.fastq ]; then
			gzip -c $dir_out/${list_file[i]}.fastq > $dir_out/${name[i]}.fastq.gz
		else
			gzip -c $dir_out/${list_file[i]}_1.fastq > $dir_out/${name[i]}_R1.fastq.gz
			rm $dir_out/${list_file[i]}_1.fastq
			gzip -c $dir_out/${list_file[i]}_2.fastq > $dir_out/${name[i]}_R2.fastq.gz
			rm $dir_out/${list_file[i]}_2.fastq
			rm -Rf $dir_out/${list_file[i]}
		fi
	fi
done
}


mapping_Fastq_bowtie2 (){
#mapping_Fastq_bowtie2 -id <PATH> -od <PATH> -s <INT> -pr <INT>
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-id)
			local in_dir=$2
			 echo "data directory set to: ${2}";shift 2;;
		-od)
			local out_dir=$2
			 echo "output directory set to: ${2}";shift 2;;
		-s)
			local seed=$2
			 echo "seed set to: ${2}";shift 2;;
		-pr)
			local proc=$2
			 echo "thread(s) number set to: ${2}";shift 2;;
		-in)
			local PATH_TO_BOWTIE_INDEX=$2
			echo "bowtie index set to: ${2}"; shift 2;;
		-h)
			usage mapping_Fastq_bowtie2; return;;
		--help)
			usage mapping_Fastq_bowtie2; return;;
		*)
			echo "Error in arguments"
			echo $1; usage mapping_Fastq_bowtie2; return;;
	esac
done
local Errors=0
if [ -z $proc ]; then echo "-pr argument not used, using 1 processor"; local proc=1; fi
if [ -z $seed ]; then echo "-s argument not used, using default seed for random (1254)"; local seed=1254; fi
if [ -z $in_dir ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ -z $PATH_TO_BOWTIE_INDEX ]; then echo "no bowtie2 index using default"; PATH_TO_BOWTIE_INDEX="/home/312.3-StrucDev/312.3.1-Commun/bowtie/bowtie2-2.3.4.1-linux-x86_64/indexes/at" ; fi
if [ $Errors -gt 0 ]; then usage mapping_Fastq_bowtie2; return 1; fi


# PATH_TO_BOWTIE_INDEX="/home/312.3-StrucDev/312.3.1-Commun/bowtie/bowtie2-2.3.4.1-linux-x86_64/indexes/at"
# PATH_TO_BOWTIE2=$(dirname $(which bowtie2))
# PATH_TO_SAMTOOLS=$(dirname $(which samtools))
# PATH_TO_BEDTOOLS2=$(dirname $(readlink -f $(which bedtools)))

echo "Processing dataset $in_dir";
echo $(find $in_dir -name "*.fastq.gz")

for fastq in $(find $in_dir -name "*.fastq.gz")
do 
	if [ -e $fastq ]
	then
		echo "fastq found"
		local fastq_file=${fastq##*/}
		local fastq_filename=${fastq_file%.*.gz}
# 		local fileIN=(${fastq_filename//_/ })
 		local filename=${fastq_filename} #get the file name
		local order=0
		
		if [[ $fastq_filename == *"_R1"* ]]; then local order=1; local filename=${fastq_filename%_R1*}; fi
		if [[ $fastq_filename == *"_R2"* ]]; then local order=2; local filename=${fastq_filename%_R2*}; fi
# 		local paired_file=${fileIN[2]%.*} #get last element which is the name of the paired file
		echo "Step 1: alignment using bowtie2";
		#no paired file
		local out=$out_dir/${out_dir##*/}  #$fastq_filename
		echo "out tmp: " $out
		if [ $order -eq 0 ] && [ ! -e "$out.sam" ] #[ "$paired_file" == "$(dirname $fastq)" ]
		then
			echo "Single end mapping"
# 			local out=$out_dir/$filename #$fastq_filename
			if [ ! -e "$out.sam" ]; 
			then
				fastqc -t $proc -o $in_dir $fastq
				$PATH_TO_BOWTIE2/bowtie2 -x $PATH_TO_BOWTIE_INDEX -U $fastq -S $out.sam -p $proc -3 16 2> ${out}_log.txt
			fi
		else
		#figure out the order of the pair-ended files
			local paired_file_order=1
			local fastqIsFirst=false
			if [ "$order" == "$paired_file_order" ]
			then
				local paired_file_order=2
				local fastqIsFirst=true
			fi
		#build up the file name of the pair
			local paired_file=${fastq/"R"$order/"R"$paired_file_order}
		#${paired_file%.*}_${paired_file_order}_$(basename $filename).fastq.gz
			echo "The paired file for $fastq is $paired_file"
# 			if [[ ${fileIN[2]} != "."* ]]; then
# 				local out=$out_dir/$(basename $filename)_${fileIN[2]%.*} #$(dirname $fastq)/$(basename $filename)_${fileIN[2]%.*}
# 				echo "here"
# 			else
				local out=$out_dir/${out_dir##*/} 
# 			fi
			echo "out: " $out
# 			mkdir -p $out
			local R_NGmerge=$out
 			fastqc -t $proc -o $in_dir $fastq $paired_file 
			# check if the paired file exists in the folder or exit with error
			if [ -e ${paired_file} ] && [ ! -e "$out.sam" ]
			then
# 				local out=$out_dir/$(basename $filename)_${fileIN[2]%.*} #$(dirname $fastq)/$(basename $filename)_${fileIN[2]%.*}
				if $fastqIsFirst
				then 
  					NGmerge -a -1 $fastq -2 $paired_file -q 34 -o $R_NGmerge -n $proc
 					$PATH_TO_BOWTIE2/bowtie2 --seed $seed -x $PATH_TO_BOWTIE_INDEX -1 ${out}_1.fastq.gz -2 ${out}_2.fastq.gz -S $out.sam --dovetail -p $proc 2> ${out}_log.txt
					#echo "$PATH_TO_BOWTIE2/bowtie2 -x $PATH_TO_BOWTIE_INDEX -1 $fastq -2 $paired_file -S $fastq.sam"
					echo "Processed file ${out}.sam"
				else
  					NGmerge -a -2 $fastq -1 $paired_file -q 34 -o $R_NGmerge -n $proc
  					$PATH_TO_BOWTIE2/bowtie2 --seed $seed -x $PATH_TO_BOWTIE_INDEX -1 ${out}_1.fastq.gz -2 ${out}_2.fastq.gz -S $out.sam --dovetail -p $proc 2> ${out}_log.txt
					echo "$PATH_TO_BOWTIE2/bowtie2 -x $PATH_TO_BOWTIE_INDEX -1 $paired_file -2 $fastq -S $fastq.sam"
					echo "Processed file ${out}.sam" 
				fi
			else
				echo "Paired file for $fastq does not exist for data set $(basename $in_dir)"
				continue
			fi
			#remove the fastq files
			rm $fastq
			rm $paired_file
		fi
		echo "Step 2: filtering SAM";
		# filter out reads with suboptimal scores
		# filter out reads with more than 2 mismatches, mapping qual <30 and with subalignement (romain 23june2022' interpretation of command below)
		$PATH_TO_SAMTOOLS/samtools view -Sh $out.sam | \
		grep -e "^@" -e 'XM:i:[012][^0-9]' | awk '$1~/@/ || $5>30 {print $0}' | grep -v "XS:i:" > $out.filtered.sam
# 		$PATH_TO_SAMTOOLS/samtools view -Sh $out.sam | \
# 		grep -v "XS:i:" > $out.filtered.sam 
		echo "Step 3: SAM to BAM conversion"
		$PATH_TO_SAMTOOLS/samtools view -Sh -b $out.filtered.sam \
		> $out.filtered.bam
		$PATH_TO_SAMTOOLS/samtools sort -o $out.filtered.sorted.bam \
		$out.filtered.bam
		echo "Step 4: remove PCR duplicates";
		$PATH_TO_SAMTOOLS/samtools rmdup -s $out.filtered.sorted.bam \
		$out.filtered.sorted.nodup.bam
		echo "Step 5: BAM indexing";
		$PATH_TO_SAMTOOLS/samtools index $out.filtered.sorted.nodup.bam
		echo "Step 6: printing mapping statistics";
		if [ $order -eq 0 ]; then
			local lines=$(zcat ${out}.fastq.gz | wc -l)
			local reads=$(expr $lines / 4)
			samtools stats $out.sam > $out.sam_stats
			samtools stats $out.filtered.sorted.nodup.bam > $out.filtered.sorted.nodup.stats
			local R1=$reads
			local RMAPPED=$(grep "reads mapped:" $out.filtered.sorted.nodup.stats | cut -f 3)
			local R2="NA"
			local MEAN1="NA"
			local SD1="NA"
			local MEAN2="NA"
			local SD2="NA"
		else
			local lines=$(zcat ${out}_1.fastq.gz | wc -l)
			local reads=$(expr $lines / 4)
			local R1=$reads
			local lines=$(zcat ${out}_2.fastq.gz | wc -l)
			local reads=$(expr $lines / 4)
			local R2=$reads
			
			samtools stats $out.sam > $out.sam_stats
			samtools stats $out.filtered.sorted.nodup.bam > $out.filtered.sorted.nodup.stats
			local RMAPPED=$(grep "reads mapped:" $out.filtered.sorted.nodup.stats | cut -f 3)
			echo "...printing insert size (PE only)"
			samtools view -f66 ${out}.filtered.sorted.nodup.bam  | cut -f9 | awk '{print sqrt($0^2)}' > ${out}.tmpIS
			local MEAN2=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' ${out}.tmpIS)
			local SD2=$(awk '{x+=$0;y+=$0^2}END{print sqrt(y/NR-(x/NR)^2)}' ${out}.tmpIS)
			rm ${out}.tmpIS
			local MEAN1=$(grep "insert size average" $out.filtered.sorted.nodup.stats | cut -f 3)
			local SD1=$(grep "insert size standard deviation" $out.filtered.sorted.nodup.stats | cut -f 3)
 		fi
 		echo "R1 R2 RMAPPED MEAN_IS SD_IS MEAN_IS_f66 SD_IS_f66" > $out.minimal.stats
		echo $R1 $R2 $RMAPPED $MEAN1 $SD1 $MEAN2 $SD2 >> $out.minimal.stats		

		echo "Step 7: removing temporary sam";
		rm $out.filtered.sam &
		rm $out.sam &
		rm ${out}.filtered.bam &
		rm ${out}.filtered.sorted.bam &
		wait
		if [ $order -eq 0 ]; then
			rm ${out}.fastq.gz
		else
			rm ${out}_1.fastq.gz
			rm ${out}_2.fastq.gz
		fi
		
		date
	fi
done

}

main_mapping_Fastq(){
# main_mapping_Fastq -fd <PATH> -md <PATH> -f1 <FILE> -f2 [FILE] -n <STRING> -s [INT] -pr [INT]
local Fastq2="NA"; local seed=1254; local proc=1; PATH_TO_BOWTIE_INDEX="/home/312.3-StrucDev/312.3.1-Commun/bowtie/bowtie2-2.3.4.1-linux-x86_64/indexes/at";
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-fd)
			local Fastq_dir=$2
			 echo "data directory set to: ${2}";shift 2;;
		-md)
			local Mapping_dir=$2
			 echo "Mapping directory set to: ${2}";shift 2;;
		-f1)
			local Fastq1=$2
			 echo "Fastq pair 1 / Fastq non paired end set to: ${2}";shift 2;;
		-f2)
			local Fastq2=$2
			 echo "Fastq pair 2 set to: ${2}";shift 2;;
		-n)
			local Name=$2
			 echo "name of Fastq directory set to: ${2}";shift 2;;
		-s)
			local seed=$2
			 echo "seed for random set to: ${2}";shift 2;;
		-pr)
			local proc=$2
			 echo "thread(s) number set to: ${2}";shift 2;;
		-in)
			local PATH_TO_BOWTIE_INDEX=$2
			echo "index set to: ${2}";shift 2;;
		-h)
			usage main_mapping_Fastq; return;;
		--help)
			usage main_mapping_Fastq; return;;
		*)
			echo "Error in arguments"
			echo $1; usage main_mapping_Fastq; return;;
	esac
done
local Errors=0
if [ -z $proc ]; then echo "-pr argument not used, using 1 processor"; fi
if [ -z $seed ]; then echo "-s argument not used, using default seed for random (1254)"; fi
if [ -z $Fastq2 ]; then echo "-f2 argument not used, assuming non paired-end analysis"; fi
if [ -z $Fastq_dir ]; then echo "ERROR: -fd argument needed"; Errors+=1; fi
if [ -z $Mapping_dir ]; then echo "ERROR: -md argument needed"; Errors+=1; fi
if [ -z $Fastq1 ]; then echo "ERROR: -f1 argument needed"; Errors+=1; fi
if [ -z $Name ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
if [ $Errors -gt 0 ]; then usage main_mapping_Fastq; return 1; fi
echo "==================="
if [[ $Fastq1 == *".done"* ]]; then
	local tmpName=$Fastq2
	echo $tmpName
	if [[ -f $Fastq_dir/${tmpName}_R1.fastq.gz ]]; then
		local Fastq1=$Fastq_dir/${tmpName}_R1.fastq.gz
	if [[ -f $Fastq_dir/${tmpName}_R2.fastq.gz ]]; then
		local Fastq2=$Fastq_dir/${tmpName}_R2.fastq.gz
	fi
	fi
	if [[ -f $Fastq_dir/${tmpName}.fastq.gz ]]; then
		local Fastq1=$Fastq_dir/${tmpName}.fastq.gz
	fi
fi

echo $Fastq1
echo $Fastq2
echo "==================="

if [ -d $Mapping_dir/$Name ]; then
  rm -Rf $Mapping_dir/$Name
fi
mkdir -p $Mapping_dir/$Name
if [ $(dirname $Fastq1) != $Mapping_dir/$Name ];then
	if [[ $Fastq1 == *"R1"* ]]; then
		cp $Fastq1 $Mapping_dir/$Name/${Name}_R1.fastq.gz
	else
		cp $Fastq1 $Mapping_dir/$Name/${Name}.fastq.gz
	fi
fi
if [ $Fastq2 != "NA" ]; then
	if [ $(dirname $Fastq2) != $Mapping_dir/$Name ];then
		if [[ $Fastq1 == *"R1"* ]]; then
			cp $Fastq2 $Mapping_dir/$Name/${Name}_R2.fastq.gz
		fi
	fi
fi



#mapping_Fastq_bowtie2 -id <PATH> -od <PATH> -s <INT> -pr <INT>
mapping_Fastq_bowtie2 -id $Mapping_dir/$Name -od $Mapping_dir/$Name -s $seed -pr $proc -in ${PATH_TO_BOWTIE_INDEX}

}


peakcalling_MACS2 (){
# peakcalling_MACS2 -id <PATH> -od <PATH> -g <FILE> -add <STRING>
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-id)
			local in_dir=$2
			 echo "data directory set to: ${2}";shift 2;;
		-od)
			local out_dir=$2
			 echo "output directory set to: ${2}";shift 2;;
		-g)
			local Genome_length=$2
			 echo "Genome length (mappable) set to: ${2}";shift 2;;
		-s)
			local seedrandom=$2
			echo "random seed set to: ${2}"; shift 2;;
		-add)
			local additionnal_args="$2"
			 echo "additional arguments for MACS2: ${2}";shift 2;;
		-h)
			usage peakcalling_MACS2; return;;
		--help)
			usage peakcalling_MACS2; return;;
		*)
			echo "Error in arguments"
			echo $1; usage peakcalling_MACS2; return;;
	esac
done
local Errors=0
if [ -z $Genome_length ]; then echo "-g argument not used, assuming A.thaliana is used: 120000000"; local Genome_length=120000000; fi
if [ -z "$additionnal_args" ]; then echo "-add argument not used, no additional argument passed to MACS2"; local additionnal_args=""; fi
if [ -z $seedrandom ]; then echo "-s missing, default random seed used"; local seedrandom=168159; fi
if [ -z $in_dir ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ $Errors -gt 0 ]; then usage peakcalling_MACS2; return 1; fi

# PATH_TO_MACS=$(dirname $(readlink -f $(which macs2)))
# PATH_TO_SAMTOOLS=$(dirname $(which samtools))
# PATH_TO_BEDTOOLS2=$(dirname $(readlink -f $(which bedtools)))

# START PROCESSING
echo "Processing folder: ${in_dir##*/}"
# add current folder to the output path
local out_dir=${out_dir}/${in_dir##*/}
#check if the output folder is not already there 
if [ ! -d "$out_dir" ]
     then
		echo "creating folder for sample" ## LT
        #create folder structure
        mkdir -m 774 -p $out_dir
        mkdir -m 774 $out_dir/controls
        mkdir -m 774 $out_dir/replicates     
    
        #split controls and replicates separately (ENCODE datasets have control in the name of the file)     
        for bam in $(find $in_dir -name "*.filtered.sorted.nodup.bam")
           do
             if [[ ! -L $bam && -f $bam ]]
             then
                 if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" == *"control"* ]]               
                 then
                    cp $bam $out_dir/controls
                 else
                    cp $bam $out_dir/replicates
# 				local libsize=$(samtools view -c $bam)
# 				local scale=$(calc 1000000/$libsize)
# 	# 		    bedtools genomecov -bga -scale $scale -ibam $bam -g /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.size > $out_dir/${in_dir##*/}_cpm.bdg
# 				bedtools genomecov -bga -scale $scale -ibam $bam > $out_dir/${in_dir##*/}_cpmold.bdg
                 fi
             else
                  echo "link detected as bam files. Those are not accepted. File is ignored"
             fi
        done
        local log=$out_dir/${in_dir##*/}.log;


        #check for controls and unzip them or log "no controls"
        local all_ok=true;
        cd $out_dir/controls
        local controls=(*)
        echo "controls: ${controls[@]}" ## LT 07/01

        if [ ${#controls[@]} -eq 0 ]
            then
               echo "No control files present"
               local all_ok=false;
            else
#               gunzip *.gz
				local controls=(*)
				#prefix with the absolute path
				local controls=("${controls[@]/#/$out_dir/controls/}")
				local controls=`join_by " " "${controls[@]}"`
				#merge all control files
				$PATH_TO_SAMTOOLS/samtools merge control.bam $controls
				#convert to BED
				$PATH_TO_BEDTOOLS2/bamToBed -i control.bam > control.bed

				#clean up
	#                rm *.bam
				$PATH_TO_BEDTOOLS2/sortBed -i control.bed > control.bed.sorted
				mv control.bed.sorted control.bed
        fi

        #check for replicates and unzip them or log "no replicates" 
        cd $out_dir/replicates
        local replicates=(*)

        if [ ${#replicates[@]} -eq 0 ]
            then
               echo "No replicates present"
               local all_ok=false;
            else
              echo "Replicates"
#               gunzip *.gz;
               local replicates=(*)
               #prefix with the absolute path       
               echo ${replicates[@]}
               local replicates=("${replicates[@]/#/$out_dir/replicates/}")
               echo ${replicates[@]}
               local replicates=`join_by " " "${replicates[@]}"`
               echo ${replicates[@]}
               #merge all replicates files
               $PATH_TO_SAMTOOLS/samtools merge replicate.bam $replicates #this step is not needed anymore because we run macs2 we treate sample by sample and then merge with mspc
               #convert to BED
               $PATH_TO_BEDTOOLS2/bamToBed -i replicate.bam > replicate.bed
               #clean up
               rm *.bam
               $PATH_TO_BEDTOOLS2/sortBed -i replicate.bed > replicate.bed.sorted
               mv replicate.bed.sorted replicate.bed
        fi

        # launch the peak calling by peak caller type
        if [ $all_ok ]
            then
		cd $out_dir
		R2_value="NA"
		for bam in $(find $in_dir -name "*.filtered.sorted.nodup.bam") # checking if bam is PE or not
		do
			if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
				R2_value=$(awk -v FS=" " 'NR>1{print $2}' ${bam%.filtered.sorted.nodup.bam}.minimal.stats)
			fi
		done
		local formatBed="BED"
		if [ ${R2_value} != "NA" ]; then
			local formatBed="BEDPE"
		fi
		
		echo "number of controls: ${#controls[@]}" ## LT
		echo "Started MACS..."
		time $PATH_TO_MACS/macs2 callpeak -t $out_dir/replicates/replicate.bed -c $out_dir/controls/control.bed -B -f $formatBed -n ${in_dir##*/} -g $Genome_length --call-summits $additionnal_args >> $log 2>&1;
       elif [ ${#controls[@]} -eq 0 ]; then
			time $PATH_TO_MACS/macs2 callpeak -t $out_dir/replicates/replicate.bed -B -f $formatBed -n ${in_dir##*/} -g $Genome_length --call-summits --seed $seedrandom $additionnal_args >> $log 2>&1;
		else
           echo "Not all needed input present."
           return 1
       fi
       
       if [ $formatBed == "BEDPE" ]; then
		fragment_length=$(cat $log | grep "fragment size = " | awk -v OFS="\t" '{print $12}')
       else
		fragment_length=$(cat $log | grep "predicted fragment length is" | awk -v OFS="\t" '{print $13}') # SE
       fi
       
#        
       
#        R2_value=$(awk -v FS=" " 'NR>1{print $2}' $in_dir/${in_dir##*/}.minimal.stats)
# 		R2_value=$(awk -v FS=" " 'NR>1{print $2}' $in_dir/*.minimal.stats)
		R2_value="NA"
		for bam in $(find $in_dir -name "*.filtered.sorted.nodup.bam")
		do
			if [[ ! -L $bam && -f $bam ]]; then
				if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
					R2_value=$(awk -v FS=" " 'NR>1{print $2}' ${bam%.filtered.sorted.nodup.bam}.minimal.stats)
				fi
				if [ ${R2_value} == "NA" ]; then #TODO correct this IF as it send an error as R2_value is not set
					if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
						bedtools bamtobed -i $bam > $in_dir/${in_dir##*/}.bamtobed.bed
						awk -v fraglen=${fragment_length} -v OFS="\t" '$6=="+"{print $1,$2,$3+fraglen,$4,$5,$6}$6=="-"{print $1,$2-fraglen,$3,$4,$5,$6}' $in_dir/${in_dir##*/}.bamtobed.bed | awk -v OFS="\t" '$2<=0{print $1,1,$3,$4,$5,$6;next}{print $0}' > $in_dir/${in_dir##*/}.ext.bed
						local libsize=$(samtools view -c $bam)
						local scale=$(calc 1000000/$libsize)
						
						sed -i 's/chr/Chr/g' $in_dir/${in_dir##*/}.ext.bed  ## Laura 27/05/2021
						bedtools genomecov -bga -scale $scale -i $in_dir/${in_dir##*/}.ext.bed -g /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.size | sed 's/Chr/chr/g' > $out_dir/${in_dir##*/}_cpm.bdg ## Laura 27/05/2021
						
						
						cp $out_dir/${in_dir##*/}_cpm.bdg $out_dir/${in_dir##*/}_cov.bdg
					fi
				else
					if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
						local libsize=$(samtools view -c $bam)
						local scale=$(calc 1000000/$libsize)
						bedtools genomecov -bga -scale $scale -ibam $bam > $out_dir/${in_dir##*/}_cpm.bdg
						cp $out_dir/${in_dir##*/}_cpm.bdg $out_dir/${in_dir##*/}_cov.bdg
					fi
				fi
			fi
		done
       
       
       
       # filter out peaks that are present in Input
       i=0
       for bdg in $(find $out_dir -name "*_treat_pileup.bdg")
         do
		   i+=1
           echo ${in_dir##*/}
           local basename=${bdg%_treat_pileup.bdg}
           local mean_cov=$(awk '{S+=$4}END{print (S/NR)*2}' ${basename}_control_lambda.bdg)
           echo ${basename#${out_dir}/} $mean_cov
           sort -k1,1 -k2,2n ${basename}_peaks.narrowPeak > ${basename}_sorted.narrowPeak
           bedtools intersect -a ${basename}_sorted.narrowPeak -b ${basename}_treat_pileup.bdg -wb | awk -v OFS="\t" '{print $11,$12,$13,$14}' > ${basename}_treatment.bdg
           bedtools intersect -a ${basename}_sorted.narrowPeak -b ${basename}_control_lambda.bdg -wb | awk -v OFS="\t" '{print $11,$12,$13,$14}' > ${basename}_control.bdg
           bedtools unionbedg -i ${basename}_treatment.bdg ${basename}_control.bdg > ${basename}_tmp.bdg
            awk -v OFS="\t" -v Mean=${mean_cov} '$5 > Mean && $2 < $3 && $4 < $5 && $5 >= 1 {print $0}' ${basename}_tmp.bdg > ${basename}_forbidden_regions.bed
            rm ${basename}_tmp.bdg
           bedtools intersect -v -a ${basename}_peaks.narrowPeak -b ${basename}_forbidden_regions.bed -wa > ${basename}_filtered.narrowPeak
           bedtools intersect -a ${basename}_peaks.narrowPeak -b ${basename}_forbidden_regions.bed -u  > ${basename}_removed.narrowPeak
         done
		if [ $i -eq 0 ]; then
			echo "ERROR no BEDGRAPH created, MACS2 must have FAIL to compute peaks"
		fi
#       clean up
       rm -R $out_dir/controls;
       rm -R $out_dir/replicates;
       awk -v OFS="\t" '{print $1,$2,$3,$4,$8}' ${basename}_filtered.narrowPeak > ${basename}_peaks.bed
	   awk -v OFS="\t" '{print $1, $2+$10-"'$phs'", $2+$10+"'$phs'"}' ${basename}_filtered.narrowPeak > ${basename}_narrow.bed
else
    echo "Folder $out_dir already exists in results. Skipping...";
fi
}


main_peakcalling (){
# main_peakcalling -id -cd -od -nc -g -add
local top=0; local mspcpc=100
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-id)
			local in_dir=("${!2}")
			 echo "data directory set to: ${in_dir[@]}";shift 2;;
		-cd)
			local input_dir=("${!2}")
			 echo "control directory set to: ${input_dir[@]}";shift 2;;
		-od)
			local out_dir=$2
			 echo "output directory set to: ${2}";shift 2;;
		-nc)
			local name_cons=$2
			 echo "name for consensus directory set to: ${2}";shift 2;;
		-g)
			local Genome_length=$2
			 echo "genome length (mappable) set to: ${2}";shift 2;;
		-add)
			local additionnal_args="$2"
			 echo "additionnal argument set to MACS2: ${2}";shift 2;;
		-top)
			local top=$2
			echo "maximum number of peaks set to: ${2}";shift 2;;
		-ps)
			local phs=$(calc $2/2.0) # $2 is the peak size! but we take only half, hence peak half size...
			echo "peak size: ${2}";shift 2;;
		-mspc)
			local mspcpc=$2
			echo "percentage for mspc set to: $2"; shift 2;;
		-s)
			local seedrandom=$2
			echo "random seed set to: ${2}"; shift 2;;
		-h)
			usage main_peakcalling; return;;
		--help)
			usage main_peakcalling; return;;
		*)
			echo "Error in arguments"
			echo $1; usage main_peakcalling; return;;
	esac
done
local Errors=0
if [ -z $Genome_length ]; then echo "-g argument not used, assuming A.thaliana is used: 120000000"; local Genome_length=120000000; fi
if [ -z $phs ]; then echo "-ps argument not used; by default the peak regions will be resized 200 bp at each side of the max height of the peak ( which correspond to ps=400)"; local phs=200; fi
if [ -z "$additionnal_args" ]; then echo "-add argument not used, no additional argument passed to MACS2"; local additionnal_args=""; fi
if [ -z $seedrandom ]; then echo "-s missing, default random seed used"; local seedrandom=168159; fi
if [ -z $in_dir ]; then echo "ERROR: -id argument needed or need to be a LIST"; Errors+=1; fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ -z $input_dir ]; then echo "ERROR: -cd argument needed or need to be a LIST"; Errors+=1; fi
if [ -z $name_cons ]; then echo "ERROR: -nc argument needed"; Errors+=1; fi
if [ $Errors -gt 0 ]; then usage main_peakcalling; return 1; fi
local list_peaks_mspc=()
local list_bdg=()

for dir in ${in_dir[@]};
do
	bam_ech=$(find $dir -name "*.filtered.sorted.nodup.bam")
  local i=1
  for bam in $(find $input_dir -name "*.filtered.sorted.nodup.bam") # TODO JEREMY: implement a loop on control dirs
  do
	if [ ! -f $dir/control_$i.filtered.sorted.nodup.bam ] || [[ $bam -nt $dir/control_$i.filtered.sorted.nodup.bam ]]; then
		cp $bam $dir/control_$i.filtered.sorted.nodup.bam
		local i=$(($i+1))
    fi
  done
  local name=${dir##*/}
  list_peaks_mspc+=("$out_dir/$name/${name}_peaks.bed")
  list_bdg+=("$out_dir/$name/${name}_cpm.bdg")
  peakcalling_MACS2 -id $dir -od $out_dir -g $Genome_length -s $seedrandom -add "$additionnal_args"

for elt in ${bam_ech[@]};
do
 	if [ ! -f $out_dir/$name/${name}_stats.txt ] || [[ $out_dir/$name/${name}_peaks.bed -nt $out_dir/$name/${name}_stats.txt ]]; then
		if [[ $elt != *"control"* ]]; then # filter out all bam that contain control in his name.
			echo $elt
			total=$(samtools view -c $elt)
			inpeak=$(bedtools sort -i $out_dir/$name/${name}_peaks.bed | bedtools merge -i stdin | bedtools intersect -u -a $elt -b stdin -ubam | samtools view -c)

			FRIP=$(calc $inpeak/$total*100)
			echo "$name $FRIP% $inpeak / $total" > $out_dir/$name/${name}_FreqReadInPeak.txt
			
			R2_value=$(awk -v FS=" " 'NR>1{print $2}' ${elt%.filtered.sorted.nodup.bam}.minimal.stats)
			if [ $R2_value == "NA" ]; then
				totTags=$(grep "total tags in treatment" $out_dir/$name/${name}.log | cut -d " " -f 15)
				filtTags=$(grep "tags after filtering in treatment" $out_dir/$name/${name}.log | cut -d " " -f 16)
			fi
			if [ $R2_value != "NA" ]; then
				totTags=$(grep "total fragments in treatment" $out_dir/$name/${name}.log | cut -d " " -f 15)
				filtTags=$(grep "fragments after filtering in treatment" $out_dir/$name/${name}.log | cut -d " " -f 16)
			fi
			
			filtPeaks=$(wc -l $out_dir/$name/${name}_filtered.narrowPeak | cut -d " " -f 1)
			echo "Sample totalTags filtTags filtPeaks FRIP" > $out_dir/$name/${name}_stats.txt
			echo $name $totTags $filtTags $filtPeaks $FRIP >> $out_dir/$name/${name}_stats.txt
		fi
	fi
done

done
if [[ -d $out_dir/$name_cons ]]; then rm -Rf $out_dir/$name_cons; fi
mkdir -p -m 774 $out_dir/$name_cons/temp
if [ "${#in_dir[@]}" -gt 1 ]; then

	if [[ $(uname -a) =~ "el7" ]]; then
		echo "${mspc_mk} -i ${list_peaks_mspc[@]} -r Tec -w 1e-4 -s 1e-8 -c 100% -o $out_dir/$name_cons -d 1"
		${mspc_mk} -i ${list_peaks_mspc[@]} -r Tec -w 1e-4 -s 1e-8 -c ${mspcpc}% -o $out_dir/$name_cons -d 1
	else
		echo "Error SL7 node is needed for mspc" && return 1
	fi


	

	if [ $top -eq 0 ]; then  
		sed '1d' $out_dir/$name_cons/ConsensusPeaks.bed | awk -v OFS="\t" '{print $1,$2,$3}' | sort -k1,1 -k2,2n > $out_dir/$name_cons/${name_cons}.bed
	else
		sed '1d' $out_dir/$name_cons/ConsensusPeaks.bed | sort -k5,5nr | awk -v OFS="\t" '{print $1,$2,$3}' > $out_dir/$name_cons/${name_cons}_comp.bed
		head -${top} $out_dir/$name_cons/${name_cons}_comp.bed | sort -k1,1 -k2,2n > $out_dir/$name_cons/${name_cons}.bed
	fi
	


	local i=0
	local files=()
	for dir in ${in_dir[@]};
	do
	local name=${dir##*/}
	if [ $i -eq 0 ];
	then # why $5+$13 ? $5 start=filtered narropeak & $13 position of  maximum
		bedtools intersect -a $out_dir/$name_cons/$name_cons.bed -b $out_dir/$name/${name}_filtered.narrowPeak -loj | awk -v OFS="\t" '{print $1,$2,$3,$5+$13}' | awk -v OFS="\t" -v chr="" -v start="" -v stop="" -v save="" 'start!=$2 && stop !=$3{if(save!=""){print save};save=$0;chr=$1;start=$2;stop=$3;next} chr==$1 && start==$2 && stop ==$3 {save=save" "$4}END{print save}' > $out_dir/$name_cons/temp/tmp_peaks_$i.bed 
	else
		bedtools intersect -a $out_dir/$name_cons/$name_cons.bed -b $out_dir/$name/${name}_filtered.narrowPeak -loj | awk -v OFS="\t" '{print $1,$2,$3,$5+$13}' | awk -v OFS="\t" -v chr="" -v start="" -v stop="" -v save="" 'start!=$2 && stop !=$3{if(save!=""){print save};save=$0;chr=$1;start=$2;stop=$3;next} chr==$1 && start==$2 && stop ==$3 {save=save" "$4}END{print save}' | awk -v OFS=" " '{$1=$2=$3="";print $0}' | sed 's/   //' > $out_dir/$name_cons/temp/tmp_peaks_$i.bed
	fi
	files+=("$out_dir/$name_cons/temp/tmp_peaks_$i.bed")
	local i=$(($i+1))
	done

	paste ${files[@]} > $out_dir/$name_cons/temp/${name_cons}_narrow.bed
	echo "nb of replicates :${#files[@]}"
	thpc=$(calc ${#files[@]}*${mspcpc}/100 )
	echo $thpc
	awk -v SIZE=$phs -v th=${thpc} -v FS="[ \t]" -v OFS="\t" 'function abs(v) {return v < 0 ? -v : v} {nbmax=0; 
	for(i=4;i<=NF;i++) { 
		if($i!= -1){ 
			if(nbmax!=0){
				ok=1;
				for(lmax in max){
					localmax=max[lmax];
					if(abs(localmax-$i)<=SIZE){ 
						moy[localmax]+=$i; nb[localmax]+=1; ok=0; break 
					}
				} 
				if(ok==1){ 
					max[i-3]=$i; moy[$i]=$i;nb[$i]=1
				} 
			} 
			if(nbmax==0){
				max[i-3]=$i; moy[$i]=$i;nb[$i]=1;nbmax=1
			} 
		}
	}; 
	for(kmax in max){
		localmax=max[kmax]; 
		if(nb[localmax]>=th){
			print($1, int(moy[localmax]/nb[localmax])-SIZE, int(moy[localmax]/nb[localmax])+SIZE)
		}; 
		delete max[kmax];delete nb[localmax]; delete moy[localmax]
	}
	}' $out_dir/$name_cons/temp/${name_cons}_narrow.bed | sed 's/\t\+/\t/g;s/^\t//' > $out_dir/$name_cons/${name_cons}_narrow.bed
	
# 	local list_bdg=`join_by " " "${list_bdg[@]}"`
	bedtools unionbedg -i ${list_bdg[@]} > $out_dir/$name_cons/${name_cons}_cov_sep.bdg
	awk -v OFS="\t" '{moy=0;nb=0;for(i=4;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/$name_cons/${name_cons}_cov_sep.bdg > $out_dir/$name_cons/${name_cons}_cov.bdg
else

	local name=${in_dir[0]##*/}
	#awk -v OFS="\t" '{print $1,$2,$3}' $out_dir/$name/${name}_filtered.narrowPeak > $out_dir/$name_cons/${name_cons}_narrow.bed
	#awk -v OFS="\t" '{print $1, $2+$10-"'$phs'", $2+$10+"'$phs'"}' $out_dir/$name/${name}_filtered.narrowPeak > $out_dir/$name_cons/${name_cons}_narrow.bed
	cp $out_dir/$name/${name}_narrow.bed $out_dir/$name_cons/${name_cons}_narrow.bed
	cp ${list_bdg[0]} $out_dir/$name_cons/${name_cons}_cov.bdg

fi

local i=0
local files=()
for dir in ${in_dir[@]};
do
  local name=${dir##*/}
  if [ $i -eq 0 ];
  then
	bedtools intersect -a $out_dir/$name_cons/${name_cons}_narrow.bed -b $out_dir/$name/${name}_cpm.bdg -wao | awk -v OFS="\t" '{print $1":"$2"-"$3,$7}' | sort -k1,1 -k2,2nr | sort -u -k1,1 > $out_dir/$name_cons/temp/tmp_max_${i}.txt
  else
    bedtools intersect -a $out_dir/$name_cons/${name_cons}_narrow.bed -b $out_dir/$name/${name}_cpm.bdg -wao | awk -v OFS="\t" '{print $1":"$2"-"$3,$7}' | sort -k1,1 -k2,2nr | sort -u -k1,1 | awk '{print $2}' > $out_dir/$name_cons/temp/tmp_max_${i}.txt
  fi 
  files+=("$out_dir/$name_cons/temp/tmp_max_${i}.txt")
  local i=$(($i+1))
done
paste ${files[@]} | sed "s/[-:]/\t/g" > $out_dir/$name_cons/${name_cons}_max.bed
awk -v OFS="\t" '{mean=0;for(i=4;i<=NF;i++) {mean+=$i};print $1,$2,$3,mean/(i-3)}' $out_dir/$name_cons/${name_cons}_max.bed | sort -k4,4n > $out_dir/$name_cons/${name_cons}_maxMean.bed
rm -Rf $out_dir/$name_cons/*/
}


replicates_comparisons (){
# main_peakcalling -id -cd -od -nc -g -add
local top=0
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=$2
			 echo "consensus peaks file set to: ${2}";shift 2;;
		-nb)
			local nbrep=$2
			 echo "number of replicates in the sample set to: ${2}";shift 2;;
		-od)
			local result=$2
			 echo "output directory set to: ${2}";shift 2;;
		-gd)
			local general_dir=$2
			 echo "general result directory set to: ${2}";shift 2;;
		-n)
			local name=$2
			 echo "name of sample set to: ${2}";shift 2;;
		-h)
			usage replicates_comparisons ; return;;
		--help)
			usage replicates_comparisons ; return;;
		*)
			echo "Error in arguments"
			echo $1; usage replicates_comparisons ; return;;
	esac
done
local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument is needed"; Errors+=1; fi
if [ -z $nbrep ]; then echo "ERROR: -nb argument is needed"; Errors+=1; fi
if [ -z $result ]; then echo "ERROR: -od argument is needed"; Errors+=1; fi
if [ -z $name ]; then echo "ERROR: -n argument is needed"; Errors+=1; fi
if [ -z $general_dir ]; then echo "ERROR: -gd argument is needed"; Errors+=1; fi
if [ $Errors -gt 0 ]; then usage replicates_comparisons ; return 1; fi

Pairwise_comps=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/bin/Pairwise_comps.R

local rep_list=()
local bdg_list=()
local remove_list=()
mkdir -p -m 774 $result/$name/IN $result/$name/NIN $result/$name/IRN
for rep in $(seq $nbrep)
do
	rep_list+=("${name}rep${rep}")
	bdg_list+=("$general_dir/Peakcalling/${name}rep${rep}/${name}rep${rep}_cpm.bdg")
	remove_list+=("$result/$name/IN/${name}rep${rep}.inter" "$result/$name/IN/${name}rep${rep}.inter.cov" "$result/$name/IN/${name}rep${rep}.inter.cov.normed" "$result/$name/NIN/${name}rep${rep}_filt.cov.bed" "$result/$name/IRN/${name}rep${rep}_filt.cov.bed")
# 	cp $result/Peakcalling/${name}rep${rep}/${name}rep${rep}_cpm.bdg $result/NIN/${name}rep${rep}_filt.cov.bed
done

# Input Normalized method # Romain 12/01/2022: in fact this is not input normalized! the difference with compute_rpkmrip_rpkmril being only the "way we count reads" i.e. Arnaud's way based on bedgraph and romain's way with bedtools coverage in bam files. the later count reads within the peak coordinate.
awk -v OFS='\t' '{print $1,$2,$3,"fake"}' $peaks | sort -k1,1 -k2,2n | sed '1ichr\tstart\tend\tname' > $result/$name/IN/tmp_peaks
add_coverage -b bdg_list[@] -n rep_list[@] -t $result/$name/IN/tmp_peaks -od $result/$name/IN
cat $result/$name/IN/table.tsv | cut -f 1,2,3,5,6,7,8,9,10,11,12,13 > $result/$name/IN/peaks_perSample_INRPM.txt
$R_36 $Pairwise_comps $result/$name/IN/ $result/$name/IN/peaks_perSample_INRPM.txt

# Not Input Normalized method
list_data=("$general_dir/Peakcalling")
list_bamdir=("$general_dir/Mapping")
compute_rpkmrip_rpkmril -p $peaks -bd list_bamdir[@] -pd list_data[@] -sn rep_list[@] -o $result/$name/NIN/ -m "inPeaks"
$R_36 $Pairwise_comps $result/$name/NIN/ $result/$name/NIN/peaks_perSample_rpkminPeaks.txt

# Input ReNormalized method
compute_rpkmrip_rpkmril -p $peaks -bd list_bamdir[@] -pd list_data[@] -sn rep_list[@] -o $result/$name/IRN/ -m "inLibs"
$R_36 $Pairwise_comps $result/$name/IRN/ $result/$name/IRN/peaks_perSample_rpkminLibs.txt

rm ${remove_list[@]}
rm $result/$name/*/tmp_peaks
rm $result/$name/NIN/tmpFiltTags.txt || echo "No $result/$name/NIN/tmpFiltTags.txt to remove !"
rm $result/$name/IRN/tmpFiltTags.txt || echo "No $result/$name/IRN/tmpFiltTags.txt to remove !"
}


initial_comparison(){
# initial_comparison -n1 -n2 -od -id -r1 -r2
local get_bedtools_cov="yes"; local filterCov=0; local filterHeight=0; local ratio1=0; local ratio2=0; local bamdir2="NA"; local data2="NA";
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-n1)
			local name1=$2
			echo "name of directory for dataset 1 is: ${2}";shift 2;;
		-n2)
			local name2=$2
			echo "name of directory for dataset 2 is: ${2}";shift 2;;
		-od)
			local result=$2
			echo "output directory set to: ${2}";shift 2;;
		-id)
			local data=$2
			echo "general data directory (i.e. peakcalling directory) set to: ${2}";shift 2;;
		-id2)
			local data2=$2
			echo "general data directory (i.e. peakcalling directory) for second dataset set to: ${2}";shift 2;;
		-r1)
			local ratio1=$2
			echo "ratio of coverage for specific peaks of dataset 1 is: ${2}";shift 2;;
		-r2)
			local ratio2=$2
			echo "ratio of coverage for specific peaks of dataset 2 is: ${2}";shift 2;;
		-f)
			local filterCov=$2
			echo "coverage must be at least >${2} in both sample to be considered as peaks";shift 2;;
		-he)
			local filterHeight=$2
			echo "height must be at least >${2} in both sample to be considered as peaks";shift 2;;
		-bd)
			local bamdir=$2
			echo "general bam directory set to: ${2}";shift 2;;
		-bd2)
			local bamdir2=$2
			echo "general bam directory for second dataset set to: ${2}";shift 2;;
		-gcov)
			local get_bedtools_cov=$2
			echo "'yes' or 'no' to compute reads count from bam at each peak: ${2}";shift 2;;
		-rep1)
			local list_rep1=("${!2}")
			echo "list of replicates names for dataset 1 in bam directory set to: ${list_rep1[@]}";shift 2;;
		-rep2)
			local list_rep2=("${!2}")
			echo "list of replicates names for dataset 2 in bam directory set to: ${list_rep2[@]}";shift 2;;
		-h)
			usage initial_comparison; return;;
		--help)
			usage initial_comparison; return;;
		*)
			echo "Error in arguments"
			echo $1; usage initial_comparison; return;;
	esac
done
local Errors=0
if [ -z $name1 ]; then echo "ERROR: -n1 argument needed"; Errors+=1; fi
if [ -z $name2 ]; then echo "ERROR: -n2 argument needed"; Errors+=1; fi
if [ -z $result ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ -z $data ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
if [ -z $bamdir ]; then echo "-bd not included, not RPKMrip computation"; fi
if [ -z $get_bedtools_cov ]; then echo "will compute reads count for each bam if -bamdir used"; fi
if [ -z $filterCov ]; then echo "no filter on coverage applied"; fi
if [ -z $filterHeight ]; then echo "no filter on peak height applied"; fi
if [ $Errors -gt 0 ]; then usage initial_comparison; return 1; fi
if [ $ratio1 -eq 0 ]; then 
	if [ $ratio2 -eq 0 ]; then 
		echo "no ratios specified, using 2 fold for both"
		local ratio1=2
		local ratio2=0.5
	else
		echo "WARNING: only ratio for $name2 specified, assuming equivalent fold ratio for $name1"
		local ratio1=$(calc 1/${ratio2})
	fi
fi
if [ -z $ratio2 ]; then 
	echo "WARNING: only ratio for $name1 specified, assuming equivalent fold ratio for $name2"
	local ratio2=$(calc 1/${ratio1})
fi
if [ $data2 == "NA" ]; then
	echo "id2 not used, assuming same Peakcalling directory for both datasets"
	local data2=$data
fi
if [ $bamdir2 == "NA" ]; then
	echo "bd2 not used, assuming same bam directory for both datasets"
	local bamdir2=$bamdir
fi


# 	export TMPDIR=/nobackup
# 	local do_plot_quant=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/Hist_cov_gen.r
# 	local merge_peaks=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/merge_all_peaks.py
# 	local compute_coverage=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/compute_coverage.py

	local out_dir=$result/${name1}_${name2}
	mkdir -p -m 774 $out_dir
	local cov1=$data/$name1/${name1}_cov.bdg
	local cov2=$data2/$name2/${name2}_cov.bdg

	###### few lines added for selecting peaks on their height:
	if [ $filterHeight -eq 0 ];then
		local peaks1=$data/$name1/${name1}_narrow.bed
    		local peaks2=$data2/$name2/${name2}_narrow.bed
	else
# 			awk -v filt=$filterHeight -v OFS="\t" '$4>filt && $5>filt && $6>filt' $data/$name1/${name1}_max.bed > $data/$name1/${name1}_narrow.heightFiltered.bed
			awk -v filt=$filterHeight -v OFS="\t" -v keep=1 '{for(i=4;i<=NF;i++) {if($i<=filt){keep=0}};if(keep==1){print $0};keep=1}' 	$data/$name1/${name1}_max.bed > $data/$name1/${name1}_narrow.heightFiltered.bed
			local peaks1=$data/$name1/${name1}_narrow.heightFiltered.bed
# 			awk -v filt=$filterHeight -v OFS="\t" '$4>filt && $5>filt && $6>filt' $data/$name2/${name2}_max.bed > 	$data/$name2/${name2}_narrow.heightFiltered.bed
			awk -v filt=$filterHeight -v OFS="\t" -v keep=1 '{for(i=4;i<=NF;i++) {if($i<=filt){keep=0}};if(keep==1){print $0};keep=1}' 	$data2/$name2/${name2}_max.bed > $data2/$name2/${name2}_narrow.heightFiltered.bed
			local peaks2=$data2/$name2/${name2}_narrow.heightFiltered.bed
	fi

    awk -v OFS="\t" -v name=$name1 '{print $1,$2,$3,name}' $peaks1  > $out_dir/${name1}_peaks.bed
    awk -v OFS="\t" -v name=$name2 '{print $1,$2,$3,name}' $peaks2  > $out_dir/${name2}_peaks.bed
    
    cat $out_dir/${name1}_peaks.bed $out_dir/${name2}_peaks.bed | sort -k1,1 -k2,2n > $out_dir/${name1}_${name2}_peaks.bed
    python $merge_peaks -f1 $name1 -f2 $name2 -o $result
    echo "merged"
    local peak_file=$out_dir/${name1}_${name2}_peaks_processed.bed
    cat $out_dir/${name1}_${name2}_peaks_merged.bed $out_dir/${name1}_${name2}_peaks_uniques.bed | sort -k1,1 -k2,2n | uniq > $peak_file
    echo "sorted"
    bedtools intersect -a $peak_file -b $cov1 -wa -wb -sorted -loj | awk  -v OFS="\t" '$6 != "-1" {print $0} $6=="-1" {print $1,$2,$3,$4,$1,$2,$3,0}' > $out_dir/$name1.inter &
    bedtools intersect -a $peak_file -b $cov2 -wa -wb -sorted -loj | awk  -v OFS="\t" '$6 != "-1" {print $0} $6=="-1" {print $1,$2,$3,$4,$1,$2,$3,0}' > $out_dir/$name2.inter &
    wait
    echo "intersected"
    python $compute_coverage -i $out_dir/$name1.inter -m 
    python $compute_coverage -i $out_dir/$name2.inter -m &
    wait 
    echo "cov computed"
    awk  '{print (1000*$5)/($3-$2)}' $out_dir/$name1.inter.cov  > $out_dir/$name1.inter.cov.normed # adjust cov by len peak
    awk  '{print (1000*$5)/($3-$2)}' $out_dir/$name2.inter.cov  > $out_dir/$name2.inter.cov.normed
    echo "norm adjusted"
	if [ $filterCov -eq 0 ];then
		paste $peak_file $out_dir/$name1.inter.cov.normed $out_dir/$name2.inter.cov.normed  |  sed "1ichr\tbegin\tend\tname\t$name1\t$name2" >  $out_dir/table_${name1}_${name2}.csv
	else
		paste $peak_file $out_dir/$name1.inter.cov.normed $out_dir/$name2.inter.cov.normed | awk -v filter=$filterCov -v OFS="\t" '$5>=filter || $6>=filter {print $0}' |  sed "1ichr\tbegin\tend\tname\t$name1\t$name2" >  $out_dir/table_${name1}_${name2}.csv
	fi
	
	if [ ! -z $bamdir ] && [ -z $list_rep1 ]; then # Bamdir defined but not list_rep
		local list_rep1=()
		local list_rep2=()
		for rep in {1..10}; do ## LT
# 		for rep in 1 2 3 4 5 6 7 8 9 10; do #$(ls $bamdir)
			echo $rep ## LT
			
# 			echo $bamdir/${name1}rep${rep} ## LT
			if [[ -d $bamdir/${name1}rep${rep} ]]; then
				list_rep1+=("${name1}rep${rep}")
			fi
			
# 			echo $bamdir/${name2}rep${rep} ## LT
			if [[ -d $bamdir2/${name2}rep${rep} ]]; then
				list_rep2+=("${name2}rep${rep}")
			fi
			
			if [[ ! -d $bamdir/${name1}rep${rep} ]] && [[ ! -d $bamdir2/${name2}rep${rep} ]]; then
# 				break 
				continue ## LT
			fi
			
		done
# 	exit 0
# 	
	fi
	if [ -z $bamdir ] && [ -z $list_rep1 ]; then # Bamdir & list_repnot defined
	/home/prog/R/R-3-5-0/bin/Rscript $do_plot_quant $result $name1 $name2 $out_dir/table_${name1}_${name2}.csv RiL
	else
	# if [ ${#list_rep1[@]} -gt 0 ]; then
		# compute_rpkmrip_rpkmril -p $out_dir/table_${name1}_${name2}.csv -bd ${bamdir} -pd ${data} -sn list_rep1[@] -m "inLibs" -o $out_dir/ -g $get_bedtools_cov
		# awk -v OFS="\t" '{moy=0;nb=0;for(i=4;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/peaks_perSample_rpkminLibs.txt > $out_dir/${name1}_RPKMrip.txt
	# fi
	
	# if [ ${#list_rep2[@]} -gt 0 ]; then
		# compute_rpkmrip_rpkmril -p $out_dir/table_${name1}_${name2}.csv -bd ${bamdir} -pd ${data} -sn list_rep2[@] -m "inLibs" -o $out_dir/ -g $get_bedtools_cov
		# awk -v OFS="\t" '{moy=0;nb=0;for(i=4;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/peaks_perSample_rpkminLibs.txt > $out_dir/${name2}_RPKMrip.txt
	# fi
	
	local list_all=("${list_rep1[@]}" "${list_rep2[@]}")
	local list_bamdir=("$bamdir" "$bamdir2")
	local list_data=("$data" "$data2")
	compute_rpkmrip_rpkmril -p $out_dir/table_${name1}_${name2}.csv -bd list_bamdir[@] -pd list_data[@] -sn list_all[@] -m "inLibs" -o $out_dir/ -g $get_bedtools_cov
	last_col_n1="$((${#list_rep1[@]}-1+4))"
	echo $last_col_n1
	first_col_n2="$(($last_col_n1+1))"
	echo $first_col_n2
	awk -v OFS="\t" -v c=$last_col_n1 '{moy=0;nb=0;for(i=4;i<=c;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/peaks_perSample_rpkminLibs.txt > $out_dir/${name1}_RPKMril.txt
	awk -v OFS="\t" -v c=$first_col_n2 '{moy=0;nb=0;for(i=c;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/peaks_perSample_rpkminLibs.txt > $out_dir/${name2}_RPKMril.txt
	
	
	compute_rpkmrip_rpkmril -p $out_dir/table_${name1}_${name2}.csv -bd list_bamdir[@] -pd list_data[@] -sn list_all[@] -m "inPeaks" -o $out_dir/ -g "no"
	awk -v OFS="\t" -v c=$last_col_n1 '{moy=0;nb=0;for(i=4;i<=c;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/peaks_perSample_rpkminPeaks.txt > $out_dir/${name1}_RPKMrip.txt
	awk -v OFS="\t" -v c=$first_col_n2 '{moy=0;nb=0;for(i=c;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/peaks_perSample_rpkminPeaks.txt > $out_dir/${name2}_RPKMrip.txt

# 	exit 0
	if [ -f $out_dir/${name1}_RPKMrip.txt ] && [ -f $out_dir/${name2}_RPKMrip.txt ] && [ -f $out_dir/${name1}_RPKMril.txt ] && [ -f $out_dir/${name2}_RPKMril.txt ]; then
		echo "name1 : ${name1}"
		echo "name2 : ${name2}"
		paste <(awk 'NR==1 && $2!="begin"{print $0}NR>1{print $0}' $out_dir/table_${name1}_${name2}.csv) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name1}_RPKMrip.txt) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name2}_RPKMrip.txt) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name1}_RPKMril.txt) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name2}_RPKMril.txt) | sed "1ichr\tbegin\tend\tname\t${name1}_RiLgenomecov\t${name2}_RiLgenomecov\t${name1}_RiP\t${name2}_RiP\t${name1}_RiL\t${name2}_RiL" > $out_dir/table_${name1}_${name2}_RiL_RiP.tsv
	fi
	echo "ça plante ici JL"
	/home/prog/R/R-3-5-0/bin/Rscript $do_plot_quantRIP $result $name1 $name2 $out_dir/table_${name1}_${name2}_RiL_RiP.tsv RiP
	#Rscript $do_plot_quant $result $name1 $name2 $out_dir/table_${name1}_${name2}_RiL_RiP.tsv RiP
	fi
	
}


analyze_decile(){
	local typeofnorm="RiP"
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-n1)
			local name1=$2
			echo "name of directory for dataset 1 is: ${2}";shift 2;;
		-n1)
			local name2=$2
			echo "name of directory for dataset 1 is: ${2}";shift 2;;
		-o)
			local out_dir=$2
			echo "ouput directory set to: ${2}"; shift 2;;
		-t)
			local typeofnorm=$2
			echo "type of normalization to use set to: ${2}"; shift 2;;
		-f)
			local table_comp=$2
			echo "comparison table set to: ${2}"; shift 2;;
		-h)
			usage initial_comparison; return;;
		--help)
			usage initial_comparison; return;;
		*)
			echo "Error in arguments"
			echo $1; usage initial_comparison; return;;
	esac
done
local Errors=0
if [ -z $name1 ]; then echo "ERROR: -n1 argument needed"; Errors+=1; fi
if [ -z $name2 ]; then echo "ERROR: -n2 argument needed"; Errors+=1; fi
if [ -z $out_dir ]; then echo "ERROR: -o argument needed"; Errors+=1; fi
if [ -z $table_comp ]; then echo "ERROR: -f argument needed"; Errors+=1; fi

if [ $Errors -gt 0 ]; then usage initial_comparison; return 1; fi




	local $name1=$1
	local $name2=$2
	
	local dir_comp=$3
	local typeofnorm=$4 # either RIP or RIL
	
	
	
	mkdir -p -m 774 ${dir_comp}/${name1}_${name2}/
	
	local total=$(wc -l ${dir_comp}/${name1}_${name2}/table_${name1}_${name2}_RiL_RiP.tsv | awk '{print $1-1}')
	local decile=$(calc $total/10 | awk '{print int($1+0.5)}')
# 	local decile=$(calc $total/10)
	echo "$total $decile"
	if [[ $typeofnorm == *"RiP"* ]] || [[ $typeofnorm == *"RIP"* ]] || [[ $typeofnorm == *"rip"* ]]; then
		awk 'NR!=1{print $1,$2,$3,$7,$8,$7/$8}' ${dir_comp}/${name1}_${name2}/table_${name1}_${name2}_RiL_RiP.tsv | sort -k6,6n | awk -v decile=$decile -v OFS="\t" '{printf "%s\t%d\t%d\t%f\t%f\t%f\t%d\n", $1,$2,$3,$4,$5,$6,((NR-1)/decile)+1}' | awk -v OFS="\t" '$7==11{print $1,$2,$3,$4,$5,$6,10; next}{print $0}' > ${dir_comp}/${name1}_${name2}/tmp_table_${name1}_${name2}.tsv
	else
		awk 'NR!=1{print $1,$2,$3,$9,$10,$9/$10}' ${dir_comp}/${name1}_${name2}/table_${name1}_${name2}_RiL_RiP.tsv | sort -k6,6n | awk -v decile=$decile -v OFS="\t" '{printf "%s\t%d\t%d\t%f\t%f\t%f\t%d\n", $1,$2,$3,$4,$5,$6,((NR-1)/decile)+1}' | awk -v OFS="\t" '$7==11{print $1,$2,$3,$4,$5,$6,10; next}{print $0}' > ${dir_comp}/${name1}_${name2}/tmp_table_${name1}_${name2}.tsv
	fi
# 	echo -e "Spacing\tDecile\tER\tIR\tDR" > ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Recap_F2.tsv
	echo -e "Spacing\tDecile\tIR\tER\tDR" > ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Recap_F2.tsv # just for old nomenclature ARFs
	for ((l=1;l<=10;l++)); do
			awk -v OFS="\t" -v th=$l '$7==th{print $0}' ${dir_comp}/${name1}_${name2}/tmp_table_${name1}_${name2}.tsv > ${dir_comp}/${name1}_${name2}/tmp_peaks.tsv
# 			compute_motif -p ${dir_comp}/${name1}_${name2}/tmp_peaks.tsv -n decile${l} -g $genome -od ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Motif -ls 600 -s 257486 -nm 5 -mim 8 -mam 10
			
			list_peaks=("${dir_comp}/${name1}_${name2}/tmp_peaks.tsv")
			list_matrices=("$ARF5_PWM")
			th=("-8" "-9" "-10")
			list_name=("decile${l}")
# 			compute_space -p list_peaks[@] -m list_matrices[@] -n list_name[@] -od ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Spacing -maxs $max_spacing -mins $min_spacing -ol ${offset_left} -or ${offset_right} -g ${genome} -th th[@]
			awk -v OFS="\t" -v decile=$l 'NR!=1{print $1,decile,$3,$7,$11}' ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Spacing/decile${l}/Zscore_stats_F2.tsv >> ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Recap_F2.tsv
			
			
	done
	$R_36 /home/312.6-Flo_Re/312.6.1-Commun/ARF-anr/DAP_052022/Heatmap.R ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Recap_F2.tsv ${dir_comp}/${name1}_${name2}/${name1}vs${name2} $name1 $name2 ${dir_comp}/${name1}_${name2}/tmp_table_${name1}_${name2}.tsv
}


comparison(){
# comparison -n -od -id -f -he
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-n)
			local names=("${!2}")
			echo "name of directory for dataset 1 is: ${names[@]}";shift 2;;
		-od)
			local result=$2
			echo "output directory set to: ${2}";shift 2;;
		-id)
			local data=$2
			echo "general data directory set to: ${2}";shift 2;;
		-f)
			local filterCov=$2
			echo "coverage must be at least >${2} in all samples to be considered as peaks";shift 2;;
		-he)
			local filterHeight=$2
			echo "height must be at least >${2} in all samples to be considered as peaks";shift 2;;
		-h)
			usage comparison; return;;
		--help)
			usage comparison; return;;
		*)
			echo "Error in arguments"
			echo $1; usage comparison; return;;
	esac
done
local Errors=0
if [ -z $names ]; then echo "ERROR: -n argument needed or need to be a LIST"; Errors+=1; fi
if [ -z $result ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ -z $data ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
if [ -z $filterCov ]; then echo "no filter on coverage applied"; local filterCov=0; fi
if [ -z $filterHeight ]; then echo "no filter on peak height applied"; local filterHeight=0; fi
if [ $Errors -gt 0 ]; then usage comparison; return 1; fi


local out_dir=$result ## LT 09/01; was $result/ before
 
mkdir -p ${out_dir}
local peaks=()
echo "here are the names: ${names[@]}" ## LT 09/01
for ((i=0;i<${#names[@]};i++)) # get all bdgs & peaks files in two lists
do
	covs+=("$data/${names[i]}/${names[i]}_cov.bdg")
	peaks+=("${out_dir}/${names[i]}_peaks.bed")
	if [ $filterHeight -eq 0 ];then
		echo ${names[i]} ## LT 09/01
		awk -v OFS="\t" -v name=${names[i]} '{print $1,$2,$3,name}' $data/${names[i]}/${names[i]}_narrow.bed  > ${out_dir}/${names[i]}_peaks.bed
	else # we filter by height of maximum of peaks
		awk -v filt=$filterHeight -v OFS="\t" -v keep=1 '{for(i=4;i<=NF;i++) {if($i<=filt){keep=0}};if(keep==1){print $0};keep=1}' 	$data/${names[i]}/${names[i]}_max.bed > $data/${names[i]}/${names[i]}_narrow.heightFiltered.bed
		awk -v OFS="\t" -v name=${names[i]} '{print $1,$2,$3,name}' $data/${names[i]}/${names[i]}_narrow.heightFiltered.bed  > ${out_dir}/${names[i]}_peaks.bed
	fi
done

local str_peaks=`join_by " " "${peaks[@]}"`


awk '{print}' "${peaks[@]}" | sort -k1,1 -k2,2n > ${out_dir}/all_peaks.bed ## LT 09/01
# cat ${str_peaks} | sort -k1,1 -k2,2n > ${out_dir}/all_peaks.bed ## commented LT 09/01

python ${merge_peaks_Nsets} -b ${out_dir}/all_peaks.bed -o ${out_dir}/merged.bed
echo "peaks merged"

local inter_cov=("${out_dir}/merged.bed")
local header=("chr" "begin" "end" "name")
for ((i=0;i<${#names[@]};i++))
do
	bedtools intersect -a ${out_dir}/merged.bed -b ${covs[i]} -wa -wb -sorted -loj | awk  -v OFS="\t" '$6 != "-1" {print $0} $6=="-1" {print $1,$2,$3,$4,$1,$2,$3,0}' > ${out_dir}/${names[i]}.inter
	python $compute_coverage -i ${out_dir}/${names[i]}.inter -m
	Noise=$(awk '{sum+=$4*($3-$2);total+=$3-$2}END{print sum/total}' ${covs[i]})
# 	awk -v noise=$Noise '{print ($5/($3-$2))-noise}' ${out_dir}/${names[i]}.inter.cov | awk '$1>=0.0{print $0;next}{print 0.0}' > ${out_dir}/${names[i]}.inter.cov.normed
	awk -v noise=$Noise '{print ((1000*$5)/($3-$2))}' ${out_dir}/${names[i]}.inter.cov | awk '$1>=0.0{print $0;next}{print 0.0}' > ${out_dir}/${names[i]}.inter.cov.normed
	inter_cov+=("${out_dir}/${names[i]}.inter.cov.normed")
	header+=("${names[i]}")
	echo "${names[i]} treated"
done
# echo "${inter_cov[@]}" ## LT 09/01
local str_inter_cov=`join_by " " "${inter_cov[@]}"`
local str_header=`join_by "	" "${header[@]}"`
echo ${str_header}
if [ $filterCov -eq 0 ];then
	paste ${inter_cov[@]} | sed "1i${str_header}" >  ${out_dir}/table_peaks.csv ## LT 10/01/2022
# 	paste ${str_inter_cov} | sed "1i${str_header}" >  ${out_dir}/table_peaks.csv
else
	paste ${str_inter_cov} | awk -v filter=$filterCov -v OFS="\t" '$5>=filter || $6>=filter {print $0}' |  sed "1ichr\tbegin\tend\tname\t$name1\t$name2" | sed "1i${str_header}" >  ${out_dir}/table_peaks.csv
fi

}


compute_rpkmrip_rpkmril(){
local getCov="yes"; 
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=$2
			echo "file containing peaks set to ${2}";shift 2;;
		-bd)
			local bam_dir=("${!2}")
			echo "path to mapping directory set to ${bam_dir[@]}";shift 2;;
		-pd)
			local peaks_dir=("${!2}")
			echo "path to peaks calling directory set to ${peaks_dir[@]}";shift 2;;
		-sn)
			local samples_names=("${!2}")
			echo "list of sample's name set to: ${samples_names[@]}";shift 2;;
		-m)
			local mode=$2
			echo "normalization mode set to 'inPeaks' or 'inLibs' (number of reads retained by MACS2) ${2}";shift 2;;
		-o)
			local out_dir=$2
			echo "out directory set to : ${2}";shift 2;;
		-g)
			local getCov=$2
			echo "optional, choose 'no' if you have already got the reads count (bedtools step): ${2}";shift 2;;
		-h)
			usage compute_rpkmrip; return;;
		--help)
			usage compute_rpkmrip; return;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_rpkmrip; return;;
	esac	
done
local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $bam_dir ]; then echo "ERROR: -bd argument needed"; Errors+=1; fi
if [ -z $peaks_dir ]; then echo "ERROR: -pd argument needed"; Errors+=1; fi
if [ -z $samples_names ]; then echo "ERROR: -sn argument needed or need to be a LIST"; Errors+=1; fi
if [ -z $out_dir ]; then echo "ERROR: -o argument needed"; Errors+=1; fi
if [ -z $mode ]; then echo "ERROR: -m argument needed"; Errors+=1; fi
if [ -z $getCov ]; then echo "reads count will be computed using bedtools"; fi
if [ $Errors -gt 0 ]; then usage compute_rpkmrip; return 1; fi	

mkdir -p -m 774 $out_dir
#rm $out_dir/tmpFiltTags.txt
touch $out_dir/tmpTotalTags.txt

awk -v OFS="\t" 'NR==1 && $2!="begin" && $2!="start"{print $1,$2,$3}NR>1{print $1,$2,$3}' $peaks > $out_dir/tmp_peaks
printf "pr -mts <(cut -f -3 $out_dir/tmp_peaks) " > $out_dir/tmp_run.sh
# cut -f -3 $peaks > $out_dir/tmp_peaks # DEPRECATED JL 09/11/2021

if [ getCov=="yes" ]; then 
	## if tmpFiltTags file is already present, remove it
	if [[ -f $out_dir/tmpTotalTags.txt ]]; then
		rm $out_dir/tmpTotalTags.txt
		touch $out_dir/tmpTotalTags.txt
	fi
	if [[ -f $out_dir/tmpFiltTags.txt ]]; then
		rm $out_dir/tmpFiltTags.txt
		touch $out_dir/tmpFiltTags.txt
	fi
	## same thing for FRIP file
	if [[ -f ${out_dir}/RIP.txt ]]; then
		rm ${out_dir}/RIP.txt
	fi
fi


for SAMP in  ${samples_names[@]}
do
	printf "\n\ncomputing reads coverage per peak using bedtools coverage...\n\n"
	echo $SAMP
	printf "<(cut -f 4 ${out_dir}/${SAMP}_filt.cov.bed) " >> $out_dir/tmp_run.sh
	
	
	for bams in ${bam_dir[@]}; do
# 		local bam_file=$(find ${bams} -name "$SAMP.filtered.sorted.nodup.bam" -type f )
		local bam_file=$(find ${bams}/$SAMP -name "*.filtered.sorted.nodup.bam" -type f | grep -v control | grep -v fakename) ## LT
		echo $bam_file ## LT
# 		exit 0
		if [ $bam_file != "" ]; then break; fi
	done
	 #this is to remove the control bam # NOTE NEW by JL 02/12/2021
	
	## check that the bam file listed is a file and not a link, and
	## remove the link from the array if present
# 	for bam in ${bam_file[@]}; do
# 		if [[ ! -L $bam && -f $bam ]];then
# 			echo "bam file is the correct one: $bam"
# 		else
# 			echo "bam_file is a link, we do not consider it"
# 			bam_file=("${bam_file[@]/$bam}")
# 			
# 		fi
# 	done
	
	
	# DEPRECATED by JL 02/12/2021
# 	local tmp241=$bam_dir/$SAMP/*.filtered.sorted.nodup.bam
# 	echo ${tmp241}
# 	local bam_file=$(echo $tmp241 | sed -e 's/ .*$//g') #this is to remove the control bam
	echo "hello"
	echo ${bam_file}
	if [ "$getCov" == "yes" ]
	then
	bedtools coverage \
        -a $out_dir/tmp_peaks \
        -b $bam_file -F 1 > ${out_dir}/${SAMP}_filt.cov.bed
		inpeak=$(bedtools sort -i $out_dir/tmp_peaks | bedtools merge -i stdin | bedtools intersect -u -a $bam_file -b stdin -ubam | samtools view -c)
		total=$(samtools view -c $bam_file)
		FRIP=$(calc $inpeak/$total*100)
		echo "${FRIP}% of reads in peaks (${inpeak}/${total})"
		echo "${inpeak}" >> ${out_dir}/RIP.txt
		
		#format filtered tags (for in libs normalization)
		
		for Peaksdir in ${peaks_dir[@]}
		do
			if [ -f $peaks_dir/$SAMP/${SAMP}_stats.txt ]; then
				cut -d " " -f 2 $peaks_dir/$SAMP/${SAMP}_stats.txt | grep -v totalTags >> $out_dir/tmpTotalTags.txt
				cut -d " " -f 3 $peaks_dir/$SAMP/${SAMP}_stats.txt | grep -v filtTags >> $out_dir/tmpFiltTags.txt
				break
			fi
 		done
	else
		printf "\n\nskipping bedtools coverage because reads coverage has already been computed\n\n"
	fi
	
done

printf "> ${out_dir}/allReps_RC.txt" >> $out_dir/tmp_run.sh
bash $out_dir/tmp_run.sh
rm $out_dir/tmp_run.sh
#cp ${out_dir}/allReps_RC.txt .
printf "\n\n\n\nR script to normalize reads count in peaks or in library/mapped running...\n\n"
/home/prog/R/R-3-5-0/bin/Rscript $rpkmrip_rpkmril ${out_dir}/allReps_RC.txt $mode $out_dir ${samples_names[@]}
printf "\n\nend of rpkmrip_rpkmril\n\n"
}


compute_rpkmrip_rpkmril_WIP(){
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=$2
			echo "file containing peaks set to";shift 2;;
		-bd)
			local bam_dir=$2
			echo "path to mapping directory set to";shift 2;;
		-pd)
			local peaks_dir=$2
			echo "path to peaks calling directory set to";shift 2;;
		-n1)
			local name_cond1=$2
			echo "name of first condition experiment set to";shift 2;;
		-n2)
			local name_cond2=$2
			echo "name of second condition experiment set to";shift 2;;
		-rn1)
			local rep_name1=("${!2}")
			echo "name of replicates in first condition set to: ${rep_name1[@]}";shift 2;;
		-rn2)
			local rep_name2=("${!2}")
			echo "name of replicates in second condition set to: ${rep_name2[@]}";shift 2;;
		-r)
			local ratio=$2
			echo "cut off ratio on fold change set to: ${2}";shift 2;;
		-o)
			local out_dir=$2
			echo "out directory set to : ${2}";shift 2;;
		-g)
			local getCov=$2
			echo "optional, choose 'no' if you have already computed the reads coverage: ${2}";shift 2;;
		-h)
			usage compute_rpkmrip_rpkmril_WIP; return;;
		--help)
			usage compute_rpkmrip_rpkmril_WIP; return;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_rpkmrip_rpkmril_WIP; return;;
	esac	
done
local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $bam_dir ]; then echo "ERROR: -bd argument needed"; Errors+=1; fi
if [ -z $peaks_dir ]; then echo "ERROR: -pd argument needed"; Errors+=1; fi
if [ -z $name_cond1 ]; then echo "ERROR: -n1 argument needed"; Errors+=1; fi
if [ -z $name_cond2 ]; then echo "ERROR: -n2 argument needed"; Errors+=1; fi
if [ -z $rep_name1 ]; then echo "ERROR: -rn1 argument needed"; Errors+=1; fi
if [ -z $rep_name2 ]; then echo "ERROR: -rn2 argument needed"; Errors+=1; fi
if [ -z $ratio ]; then echo "ERROR: -r argument needed"; Errors+=1; fi
if [ -z $out_dir ]; then echo "ERROR: -o argument needed"; Errors+=1; fi
if [ -z $getCov ]; then echo "coverage will be computed using bedtools"; local getCov="yes"; fi

# TODO? pour Romain, check si -rn1 & -rn2 ont la meme taille ?

if [ $Errors -gt 0 ]; then usage main_peakcalling; return 1; fi	
ratio1=$ratio
ratio2=$(calc 1/$ratio1)

# get number of replicates per condition:
num_rep1=${#rep_name1[@]}
num_rep2=${#rep_name2[@]}
printf "first condition has ${num_rep1} replicates\n"
printf "second condition has ${num_rep2} replicates\n"

local all_rep=("${rep_name1[@]}" "${rep_name2[@]}")
echo ${all_rep[@]}
mkdir -p $out_dir
#rm $out_dir/tmpFiltTags.txt
touch $out_dir/tmpFiltTags.txt
printf "pr -mts <(cut -f -3 ${peaks}) " > tmp_run.sh
cut -f -3 $peaks > tmp_peaks # this is to ensure we keep only the 3 first col of the bed
for REP in  ${all_rep[@]}
do
	echo $REP
	printf "<(cut -f 4 ${out_dir}/${REP}_filt.cov.bed) " >> tmp_run.sh
	local tmp241=$bam_dir/$REP/*.filtered.sorted.nodup.bam
	echo $tmp241
	local bam_file=$(echo $tmp241 | sed -e 's/^.* //g') #this is to remove the control bam
	echo $bam_file
	if [ "$getCov" == "yes" ]
	then
	printf "\n\ncomputing reads coverage per peak using bedtools coverage...\n\n"
	bedtools coverage \
        -a tmp_peaks \
        -b $bam_file -F 1 > $out_dir/${REP}_filt.cov.bed
		inpeak=$(bedtools sort -i $peaks | bedtools merge -i stdin | bedtools intersect -u -a $bam_file -b stdin -ubam | samtools view -c)
		total=$(samtools view -c $bam_file)
		FRIP=$(calc $inpeak/$total*100)
		echo $inpeak
		echo $total
		echo $FRIP
		
		#format filtered tags (for in libs normalization)
		#cut -d " " -f 3 $peaks_dir/$REP/*_stats.txt | grep -v filtTags >> $out_dir/tmpFiltTags.txt
	else
		printf "\n\nskipping bedtools coverage because reads coverage has already been computed\n\n"
	fi
	
done
printf "> ${out_dir}/allReps_RC.txt" >> tmp_run.sh
bash tmp_run.sh
rm tmp_run.sh
cp ${out_dir}/allReps_RC.txt .
nb_samp=$((num_rep1+num_rep2))
$norma_cov_plots $ratio1 $ratio2 inPeaks $out_dir $num_rep1 $num_rep2 $name_cond1 $name_cond2 ${rep_name1[@]} ${rep_name2[@]}
convert -density 144 $out_dir/pairwise_rep_coverage_inPeaks.pdf $out_dir/pairwise_rep_coverage_inPeaks.png
convert -density 144 $out_dir/conditions_coverage_inPeaks.pdf $out_dir/conditions_coverage_inPeaks.png
rm allReps_RC.txt
printf "\n\nend of compute_rpkmrip_rpkmril_WIP\n\n"
}


compute_motif(){

local pal=false; local motifstffm=0; local learning_size=600; local min_motif=6; local max_motif=15; local seed=1254; local nb_motifs=1; local col_coverage=0; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local background_file="/home/312.6-Flo_Re/312.6.1-Commun/data/A_thaliana_phytozome_v12/PWMscan_root/araTha1/bg_Arath1.txt"; local selection_tffm=0; local transpose=false; local topPeaks=0
# compute_motif -p -n -g -od -ls -nm -mim -mam -pal -s
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=$2
			echo "peak file set to: ${2}";shift 2;;
		-c)
			local col_coverage=$2
			echo "column number, in the input bed file (e.g. output by initial_comp), that gives sequences coverage (RPKM): ${2}";shift 2;;
		-n)
			local name=$2
			echo "name for sub-directory set to: ${2}";shift 2;;
		-g)
			local genome=$2
			echo "genome to use set to: ${2}";shift 2;;
		-od)
			local results_motif=$2
			echo "result directory set to: ${2}";shift 2;;
		-ls)
			local learning_size=$2
			echo "size for learning set: ${2}";shift 2;;
		-nm)
			local nb_motifs=$2
			echo "number of PWM motif generated: ${2}";shift 2;;
		-t)
			local transpose=true
			echo "PWM motif used will be transposed as required";shift 1;;
		-mim)
			local min_motif=$2
			echo "minimum length for motif set to: ${2}";shift 2;;
		-mam)
			local max_motif=$2
			echo "maximum length for motif set to: ${2}";shift 2;;
		-pal)
			local pal=true
			echo "palindromic mode activated";shift 1;;
		-stffm)
			local selection_tffm=$2
			echo "pfm used for tffm computation set to: ${2}"; shift 2;;
		-s)
			local seed=$2
			echo "seed for random set to: ${2}";shift 2;;
		-top)
			local topPeaks=$2
			echo "Maximum number of peaks to consider: $2";shift 2;;
		-bg)
			local background_file=$2
			echo "background file set to: ${2}";shift 2;;
		-h)
			usage compute_motif; return;;
		--help)
			usage compute_motif; return;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_motif; return;;
	esac
done
local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $name ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
if [ -z $results_motif ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ -z $learning_size ]; then echo "-ls argument not used, using 600 peaks"; fi
if [ -z $min_motif ]; then echo "-mim argument not used, min size of motif is 6"; fi
if [ -z $max_motif ]; then echo "-mam argument not used, max size of motif is 15"; fi
if [ -z $nb_motifs ]; then echo "-nm argument not used, searching for 1 motif"; fi
if [ -z $seed ]; then echo "-s argument not used, using default seed for random (1254)"; fi
if [ -z $col_coverage ]; then echo "-c argument not used, no"; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; fi
if [ -z $background_file ]; then echo "-bg not used, assuming you're using Ath"; fi
if ! $pal ; then echo "palindromic mode not activated"; fi
if ! $transpose ; then echo "transpose mode not activated"; fi
if [ $Errors -gt 0 ]; then usage compute_motif; return 1; fi

if [ $max_motif -lt $min_motif ]; then
	echo "WARNING: things got mixed up, max size of motif can't be lower than min size."
	echo -e "\tUsing min size as max size "
	local max_motif=$(calc $max_motif+$min_motif) # max motif = total
	local min_motif=$(calc $max_motif-$min_motif) # new min_motif = old max_motif (total - old min = old max)
	local max_motif=$(calc $max_motif-$min_motif) # new max_motif = old min_motif (total - new min = old min)
fi

	# prog used here
	# meme_prog=/home/prog/meme/meme_4.12.0/bin/meme-chip
	# meme2meme=/home/prog/meme/meme_4.12.0/bin/meme2meme
	# meme2pfm=/home/312.6-Flo_Re/312.6.1-Commun/scripts/meme2pfm.sh
	# # meme2pfm=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis/meme2pfm_modifJM.sh
	# prepMEMEforPalTFFM=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/prepMEMEforPalTFFM.py
	# pfmTOtffm=/home/312.6-Flo_Re/312.6.1-Commun/LFY/scripts/get_tffm.py
	
	
	
	
	local lowlimit=$(calc $learning_size/2)
	mkdir -p $results_motif/$name/sets $results_motif/$name/meme $results_motif/$name/tffm
	if [ $(wc -l $peaks | awk '{print $1}') -le $(calc $learning_size+$lowlimit) ]; then
		echo "not enough sequences for proper training/testing"
		awk -v OFS='\t' '{print $1, $2, $3}' $peaks > $results_motif/$name/sets/${name}_testingset.bed
		awk -v OFS='\t' '{print $1, $2, $3}' $peaks > $results_motif/$name/sets/${name}.bed
	else

		if [ $col_coverage != 0 ]; then
			printf "\n\n BED file will be sorted according to column $col_coverage ! \n\n"
			echo "$learning_size"
			cat $peaks | cut -f 1,2,3,$col_coverage \
			| sort -k4,4rn > $results_motif/$name/sets/tmp${name}.bed
			head -n $learning_size $results_motif/$name/sets/tmp${name}.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}.bed
			cat $peaks | cut -f 1,2,3,$col_coverage \
			| sort -k4,4rn | sed "1,${learning_size}d" > $results_motif/$name/sets/tmp${name}.bed
			if [ $topPeaks != 0 ]; then
				head -n $topPeaks $results_motif/$name/sets/tmp${name}.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}_testingset.bed
			else
				cat $results_motif/$name/sets/tmp${name}.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}_testingset.bed
			fi
			echo "testing set done"
		else
			shuf $peaks | head -n $learning_size | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}.bed
			if [ $topPeaks != 0 ]; then
				sed "1,${learning_size}d" $peaks | head -n $topPeaks | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}_testingset.bed
			else
				sed "1,${learning_size}d" $peaks | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}_testingset.bed
			fi
		fi
	fi
    #head -$learning_size $peaks | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}.bed
	
	
	bedtools getfasta -fi $genome -fo $results_motif/$name/sets/${name}.fas -bed $results_motif/$name/sets/${name}.bed
# 	/home/312.6-Flo_Re/312.6.1-Commun/data/meme_db/motif_databases/JASPAR/JASPAR2018_CORE_plants_non-redundant.meme
	if $pal ; then
		$meme_prog -oc $results_motif/$name/meme -nmeme $learning_size -meme-pal -meme-maxsize $(calc $learning_size*1000) -meme-minw $min_motif -meme-maxw $max_motif -meme-nmotifs $nb_motifs -dreme-m 0 -noecho $results_motif/$name/sets/${name}.fas -seed $seed -db /home/312.6-Flo_Re/312.6.1-Commun/data/meme_db/motif_databases/JASPAR2022_CORE_nonredundant_pfm.meme
	else
		$meme_prog -oc $results_motif/$name/meme -nmeme $learning_size -meme-maxsize $(calc $learning_size*1000) -meme-minw $min_motif -meme-maxw $max_motif -meme-nmotifs $nb_motifs -dreme-m 0 -noecho $results_motif/$name/sets/${name}.fas -seed $seed -db /home/312.6-Flo_Re/312.6.1-Commun/data/meme_db/motif_databases/JASPAR2022_CORE_nonredundant_pfm.meme
	fi
	$meme2meme $results_motif/$name/meme/meme_out/meme.txt > $results_motif/$name/meme/meme_out/meme_mini.txt
	
	local path_to_meme_mini=$results_motif/$name/meme/meme_out/meme_mini.txt

	bash $meme2pfm $path_to_meme_mini $name #> $results_motif/$name/${name}.pfm

	if $pal ; then
		printf "\n\ngenerate TFFM learning set\n\n"
		python $prepMEMEforPalTFFM -m $results_motif/$name/meme/meme_out/meme.txt -f $results_motif/$name/sets/${name}_tffm_learningset.fas -p ${selection_tffm} # THIS SCRIPT NEEDS TO BE CORRECTED 
		local learning_set_tffm=$results_motif/$name/sets/${name}_tffm_learningset.fas
# 		local learning_set_tffm=$results_motif/$name/sets/${name}.fas
# 		cat $learning_set_tffm | grep ">" -v | tr "A" "@" | tr 'T' '@' | tr 'C' '@' | tr 'G' '@'  | sed 's/@//g' | tr '\n' 'U'
	else
		local learning_set_tffm=$results_motif/$name/sets/${name}.fas
	fi
	printf "\n\npfmTotffm\n\n"
	$Python_TFFM $pfmTOtffm -r $results_motif/$name/tffm/ -f $learning_set_tffm -m $results_motif/$name/meme/meme_out/meme.txt -p ${selection_tffm}

	local pfmtokeep=$(calc $selection_tffm+1) # keeping the pfm used for tffm for spacing and ROCs computation
	
	
	cp $results_motif/$name/tffm/tffm_first_order.xml $results_motif/$name/${name}_tffm.xml
	if $transpose ; then
		head -1 $results_motif/$name/meme/meme_out/Motif_MEME_seperateFiles/Motif_${pfmtokeep}.pfm > $results_motif/$name/${name}.pfm
		echo -e 'A\tC\tG\tT' >> $results_motif/$name/${name}.pfm
		tail -n +3 $results_motif/$name/meme/meme_out/Motif_MEME_seperateFiles/Motif_${pfmtokeep}.pfm | tac | awk -v OFS="\t" '{print $4,$3,$2,$1}' >> $results_motif/$name/${name}.pfm
		cp $results_motif/$name/meme/meme_out/logo_rc${pfmtokeep}.png $results_motif/$name/${name}_logo.png
	else
		cp $results_motif/$name/meme/meme_out/logo${pfmtokeep}.png $results_motif/$name/${name}_logo.png
		cp $results_motif/$name/meme/meme_out/Motif_MEME_seperateFiles/Motif_${pfmtokeep}.pfm $results_motif/$name/${name}.pfm
	fi
}


prep_annotation(){
local promoterLength=1000; local Size="NA"; local Genome="NA"
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-n)
			local gff=$2
			echo "gff file of the genome set to: ${2}";shift 2;;
        -p)
			local promoterLength=$2
			echo "size for the promoters set to: ${2}";shift 2;;
        -g)
			local Genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
		-s)
			local Size=$2
			echo "size file of the genome (fai) set to: ${2}";shift 2;;
		-o)
			local out_file=$2
			echo "output file & directory set to: ${2}";shift 2;;
		-h)
			usage prep_annotation; return;;
		--help)
			usage prep_annotation; return;;
		*)
			echo "Error in arguments"
			echo $1; usage prep_annotation; return;;
	esac
done

local Errors=0
if [ -z $gff ]; then echo "ERROR: -n argument needed";Errors+=1;fi
if [ -z $out_file ]; then echo "ERROR: -o argument needed";Errors+=1;fi
if [ $Genome == "NA" ] && [ $Size == "NA" ]; then 
	echo "ERROR: -g or -s argument needed";Errors+=1;
else
	if [ $Size == "NA" ]; then echo "WARNING: -s argument not used, if fai file is not in same directory as genome file, will be generated"; local Size=${Genome}.fai ;fi
fi
if [ -z $promoterLength ]; then echo "-p argument not used, using 1000bp as promoter length" ;fi
if [ $Errors -gt 0 ]; then usage prep_annotation; return 1; fi

# prepare_gff=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/prepare_gff.py
# generate_BedFromGff=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/generate_bed_gff.py
# bedops=/home/312.3-StrucDev/312.3.1-Commun/bedops/bin/bedops
# bedmap=/home/312.3-StrucDev/312.3.1-Commun/bedops/bin/bedmap

## OLD / DEPRECATED
# python $prepare_gff -gi ${bed_genome}.gff3 -go ${bed_genome}_prep.gff3  -f ${Genome}
# python $generate_BedFromGff -g ${bed_genome}_prep.gff3 -p $promoterLength -o ${bed_genome}.bed -s $Size
# bedtools sort -i ${bed_genome}_prep.bed | uniq -u | $bedops --partition - | $bedmap --echo --echo-map-id --delim '\t' - ${bed_genome}_prep.bed | awk -v OFS="\t" '$4==""{print $1,$2,$3,"intergenic";next} {print $0}' > ${bed_genome}.bed

# TODO Adapt to every genome 
## NEW
local out_dir=$(dirname $out_file)
mkdir -p -m 774 $out_dir

# awk -v OFS="\t" 'split($9,a,";") a[1] ~ ".1.TAIR10$" || a[2] ~ ".1.TAIR10$" {print $0;next} $3=="gene"{print $0;next}' $gff | uniq > $out_dir/tmpPrimaryRNA.gff3 # DEPRECATED

# this awk should do the work to adapt to 'most' gff3, tested with A.thaliana, O. sativa and B. oleracea 
awk -v OFS="\t" '{split($9,a,";"); split(a[2],b,"."); if(b[2]~"1"){print $0}else{if($3=="gene"){print $0}}}' $gff | uniq > $out_dir/tmpPrimaryRNA.gff3

if [ $Size == ${Genome}.fai ]; then
	samtools faidx $Genome
	awk 'OFS="\t" {print $1, $2}' ${Genome}.fai | sort -k1,1 -k2,2n > $out_dir/chromSizes.bed
else
	awk 'OFS="\t" {print $1, $2}' $Size | sort -k1,1 -k2,2n > $out_dir/chromSizes.bed
fi

if [ $promoterLength -gt 0 ]; then
	awk -v OFS="\t" -vFS="[=\t]" -v lengthprom=$promoterLength '$1 ~ /^#/ {next} $7=="+" && $3=="gene" {print $1,$2,$3,$4,$5,$6,$7,$8,$9"="$10"="$11; print $1,".","promoter",$4-lengthprom-1,$4-1,".","+",".","ID="$11".1.TAIR10.PROMOTER;Parent="$11".TAIR10" } $7=="-" && $3=="gene" {print $1,$2,$3,$4,$5,$6,$7,$8,$9"="$10"="$11; print $1,".","promoter",$5+1,$5+1+lengthprom,".","-",".","ID="$11".1.TAIR10.PROMOTER;Parent="$11".TAIR10" }' $out_dir/tmpPrimaryRNA.gff3 | awk -v OFS="\t" '$3 != "gene" {print $0;next}' | awk -v OFS="\t" '$4<0{print $1,".","promoter","1",$5,$6,$7,$8,$9;next}{print $0}' > $out_dir/promoters.gff3
	
	cat $out_dir/tmpPrimaryRNA.gff3 $out_dir/promoters.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > $out_dir/in_sorted.gff
	bedtools complement -i <(cat $out_dir/tmpPrimaryRNA.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}') -g $out_dir/chromSizes.bed > $out_dir/intergenicTMP.bed
else
	# Option exist in case you don't want to use promoters in NSG. 
	cat $out_dir/tmpPrimaryRNA.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > $out_dir/in_sorted.gff
fi



bedtools complement -i $out_dir/in_sorted.gff -g $out_dir/chromSizes.bed > $out_dir/intergenic.bed



awk -v OFS="\t" '$1 ~ /^#/ {next} $3!= "gene" && $3!="mRNA" && $3!="promoter" {print $1, $4-1, $5+1;next}' $out_dir/in_sorted.gff > $out_dir/exon.bed

if [ $promoterLength -eq 0 ]; then
	bedtools complement -i <(cat $out_dir/exon.bed $out_dir/intergenic.bed | sort -k1,1 -k2,2n) -g $out_dir/chromSizes.bed > $out_dir/intron.bed
else
	bedtools complement -i <(cat $out_dir/exon.bed $out_dir/intergenicTMP.bed | sort -k1,1 -k2,2n) -g $out_dir/chromSizes.bed > $out_dir/intron.bed
fi
cat <(awk -v OFS="\t" '{print $0,"intron"}' $out_dir/intron.bed) <(awk -v OFS="\t" '{print $0,"intergenic"}' $out_dir/intergenic.bed) <(awk -v OFS="\t" '$1 ~ /^#/ {next} $3!= "gene" && $3!="mRNA"{print $1, $4, $5,$3;next}' $out_dir/in_sorted.gff) | awk -v OFS="\t" '$2<0{print $1,"1",$3,$4;next}{print $0}' | awk -v OFS="\t" '$4=="CDS" {print $1,$2,$3,"exon";next}{print $0}' | awk -v OFS="\t" '$2!="0" && $4=="intergenic"{print $1,$2+1,$3,$4;next}{print $0}' | uniq | sort -k1,1 -k2,2n > $out_file

[ -f $out_dir/tmpPrimaryRNA.gff3 ] && rm $out_dir/tmpPrimaryRNA.gff3
[ -f $out_dir/promoters.gff3 ] && rm $out_dir/promoters.gff3
[ -f $out_dir/in_sorted.gff ] && rm $out_dir/in_sorted.gff
[ -f $out_dir/chromSizes.bed ] && rm $out_dir/chromSizes.bed
[ -f $out_dir/intergenicTMP.bed ] && rm $out_dir/intergenicTMP.bed
[ -f $out_dir/intergenic.bed ] && rm $out_dir/intergenic.bed
[ -f $out_dir/exon.bed ] && rm $out_dir/exon.bed
[ -f $out_dir/intron.bed ] && rm $out_dir/intron.bed

}


compute_NS(){
# compute_NS -p -n -g -od -anf -nb -ws -gc -lt -s
local window_size=300; local number_of_NS=1; local seed=1254; local deltaGC=0.03; local limit_type=1000;
local timeout=false; local simple_type=false; 
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=$2
			echo "peaks file set to: ${2}";shift 2;;
		-n)
			local name=$2
			echo "prefix for results set to: ${2}";shift 2;;
		-g)
			local genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
		-od)
			local results=$2
			echo "output directory set to: ${2}";shift 2;;
		-anf)
			local annotation_file=$2
			echo "Annotation file set to: ${2}";shift 2;;
		-nb)
			local number_of_NS=$2
			echo "number of negative sets to create set to: ${2}";shift 2;;
		-ws)
			local window_size=$2
			echo "window size set to: ${2}";shift 2;;
		-gc)
			local deltaGC=$2
			echo "delta GC allowed set to: ${2}";shift 2;;
		-lt)
			local limit_type=$2
			echo "number of region by origin set to: ${2}";shift 2;;
		-s)
			local seed=$2
			echo "seed for random set to: ${2}";shift 2;;
		--timeout)
			local timeout=true
			echo "Timeout option active: NSG will start as default, but if it takes too long (>15 min per negative set demanded) it will be launched again with simple_type option";shift 1;;
		-st)
			local simple_type=true
			echo "simple_type option active: compute_NS will start with -st option";shift 1;;
		-h)
			usage compute_NS; return;;
		--help)
			usage compute_NS; return;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_NS; return;;
	esac
done
local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $name ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas" ; fi
if [ -z $results ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ -z $annotation_file ]; then echo "ERROR: -anf argument needed"; Errors+=1; fi
if [ -z $number_of_NS ]; then echo "-nb argument not used, searching for 1 negative set"; fi
if [ -z $window_size ]; then echo "-ws argument not used, using 250 bp as window size"; fi
if [ -z $deltaGC ]; then echo "-gc argument not used, using 0.03 as maximum of GC% divergence"; fi
if [ -z $limit_type ]; then echo "-lt argument not used, using 1000 as minimum group size for origin (use -h for more information)"; fi
if [ -z $seed ]; then echo "-s argument not used, using default seed for random (1254)"; fi
if [ $Errors -gt 0 ]; then usage compute_NS; return 1; fi
if $timeout ; then echo "timeout mode activated"; fi
if $simple_type ; then echo "simple_type mode activated"; fi


mkdir -p -m 774 ${results}
# awk -v OFS='\t' ' ($3-$2)>1200 {print $1, ($2+($3-$2)-600), ($2+($3-$2)+600);next} {print $1, $2, $3}' $peak > ${results}/${name}_fil.bed


if $timeout ; then # if --timeout option specified in compute_NS
	local time="$(calc 15*60*$number_of_NS)s" # timeout at 15 minutes per set
	echo $time 
	timeout $time python $negative_set_script -pos $peaks -of ${name} -od ${results} -fas  $genome -bed $annotation_file -n $number_of_NS -r $seed -GC $deltaGC -l $limit_type -bs $window_size
	lenpos=$( wc -l ${results}/${name}_pos.bed | awk '{print $1}' )
	lenneg=$( wc -l ${results}/${name}_1_neg.bed | awk '{print $1}' )
	local exit_status=$?
	if [[ $exit_status -eq 124 ]] || [[ $lenpos -ne $lenneg ]] ; then # if command times out, automatically run with faster option -st
		echo -e "\nRunning NSG in default mode took too long; running with --simpletype faster option instead\n"
		python $negative_set_script -pos $peaks -of ${name} -od ${results} -fas  $genome -bed $annotation_file -n $number_of_NS -r $seed -GC $deltaGC -l $limit_type -bs $window_size -st
	else
		echo "NSG ran with standard method"
	fi
elif $simple_type ; then # if -st option specified in compute_NS
	echo "simple_type option specified"
	python $negative_set_script -pos $peaks -of ${name} -od ${results} -fas  $genome -bed $annotation_file -n $number_of_NS -r $seed -GC $deltaGC -l $limit_type -bs $window_size -st
else
	echo "no timeout nor simple_type option specified, will run in default mode"
	python $negative_set_script -pos $peaks -of ${name} -od ${results} -fas  $genome -bed $annotation_file -n $number_of_NS -r $seed -GC $deltaGC -l $limit_type -bs $window_size
fi

}



compute_DNAshape(){
	local window_size=300; local seed=1254; local deltaGC=0.03; local limit_type=1000;
local timeout=false; local extend=0; local col_coverage=0; local topPeaks=0
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-f)
				local foreground=${2};shift 2;;
			-m)
				local matrix=${2};shift 2;;
			-e)
				local extend=${2};shift 2;;
			-o)
				local output=${2};shift 2;;
			-ls)
				local learning_size=${2};shift 2;;
			-n)
				local name=${2};shift 2;;
			-c)
				local col_coverage=${2};shift 2;;
			-g)
				local genome=$2;shift 2;;
			-anf)
				local annotation_file=$2;shift 2;;
			-ws)
				local window_size=$2;shift 2;;
			-gc)
				local deltaGC=$2;shift 2;;
			-lt)
				local limit_type=$2;shift 2;;
			-top)
				local topPeaks=$2
				echo "Maximum number of peaks to consider: $2";shift 2;;
			-s)
				local seed=$2;shift 2;;
		esac
	done
	
	local Errors=0
	if [ -z $foreground ]; then echo "ERROR: no foreground/peak bed file specified, use -f"; Errors+=1; fi
	if [ -z $matrix ]; then echo "ERROR: no matrix given, use -m (only pfm or TFFM accepted)"; Errors+=1; fi
	if [ -z $output ]; then echo "ERROR: no output prefix (PATH+prefix) given, use -o"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage compute_NS; return 1; fi
	
	local lowlimit=$(calc $learning_size/2)
	mkdir -p -m 774 $output/$name/testing $output/$name/training
	if [ $(wc -l $foreground | awk '{print $1}') -le $(calc $learning_size+$lowlimit) ]; then
		echo "not enough sequences for proper training/testing"
		awk -v OFS='\t' '{print $1, $2, $3}' $foreground > $output/$name/testing/foreground_testing.bed		
		awk -v OFS='\t' '{print $1, $2, $3}' $foreground > $output/$name/training/foreground_training.bed
	else
		if [ $col_coverage != 0 ]; then
			printf "\n\n BED file will be sorted according to column $col_coverage ! \n\n"
			if [ $topPeaks != 0 ]; then
				
				cat $foreground | cut -f 1,2,3,$col_coverage \
				| sort -k4,4rn > $output/$name/training/tmpOrdered.bed
			
				head -n $learning_size $output/$name/training/tmpOrdered.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/training/foreground_training.bed
				
				sed "1,${learning_size}d" $output/$name/training/tmpOrdered.bed > $output/$name/training/tmpOrdered2.bed
				head -n $topPeaks $output/$name/training/tmpOrdered2.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/testing/foreground_testing.bed
				rm $output/$name/training/tmpOrdered2.bed
			else
				cat $foreground | cut -f 1,2,3,$col_coverage \
				| sort -k4,4rn > $output/$name/training/tmpOrdered.bed
			
				head -n $learning_size $output/$name/training/tmpOrdered.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/training/foreground_training.bed
				cat $foreground | cut -f 1,2,3,$col_coverage \
				| sort -k4,4rn | sed "1,${learning_size}d" | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/testing/foreground_testing.bed
			fi
			
			
			rm $output/$name/training/tmpOrdered.bed
		else
			
			shuf $foreground > $output/$name/training/tmpOrdered.bed
			if [ $topPeaks != 0 ]; then
				head -n $learning_size $output/$name/training/tmpOrdered.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/training/foreground_training.bed
				sed "1,${learning_size}d" $output/$name/training/tmpOrdered.bed | head -n $topPeaks | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/testing/foreground_testing.bed
			else
				head -n $learning_size $output/$name/training/tmpOrdered.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/training/foreground_training.bed
				sed "1,${learning_size}d" $output/$name/training/tmpOrdered.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/testing/foreground_testing.bed
			fi
			
			
			
			
			
			rm $output/$name/training/tmpOrdered.bed
		fi
	fi
	compute_NS -p $output/$name/training/foreground_training.bed -n $name -g $genome -od $output/$name/training -anf $annotation_file -nb 1 -ws $window_size -gc $deltaGC -lt $limit_type -s $seed --timeout
	
	awk -v OFS='\t' '{print $1,$2,$3,$1":"$2"-"$3}' $output/$name/training/${name}_pos.bed > $output/$name/training/foreground_training.bed
	awk -v OFS='\t' '{print $1,$2,$3,$1":"$2"-"$3}' $output/$name/training/${name}_1_neg.bed > $output/$name/training/background_training.bed
	
	compute_NS -p $output/$name/testing/foreground_testing.bed -n $name -g $genome -od $output/$name/testing -anf $annotation_file -nb 1 -ws $window_size -gc $deltaGC -lt $limit_type -s $seed --timeout
	
	awk -v OFS='\t' '{print $1,$2,$3,$1":"$2"-"$3}' $output/$name/testing/${name}_pos.bed > $output/$name/testing/foreground_testing.bed
	awk -v OFS='\t' '{print $1,$2,$3,$1":"$2"-"$3}' $output/$name/testing/${name}_1_neg.bed > $output/$name/testing/background_testing.bed
	
	bedtools getfasta -fo $output/$name/testing/background_testing.fas -fi $genome -bed $output/$name/testing/background_testing.bed
	bedtools getfasta -fo $output/$name/training/background_training.fas -fi $genome -bed $output/$name/training/background_training.bed
	bedtools getfasta -fo $output/$name/testing/foreground_testing.fas -fi $genome -bed $output/$name/testing/foreground_testing.bed
	bedtools getfasta -fo $output/$name/training/foreground_training.fas -fi $genome -bed $output/$name/training/foreground_training.bed
	
	
	if [[ $matrix == *".xml"* ]]; then # TFFM scores computation
		echo "Training a first order TFFM + DNA shape classifier.";
		$Python_TFFM $ComputeDNAshaped trainTFFM -T $matrix \
			-i $output/$name/training/foreground_training.fas -I $output/$name/training/foreground_training.bed \
			-b $output/$name/training/background_training.fas -B $output/$name/training/background_training.bed \
			-o $output/$name/${name}_fo_classifier -t first_order \
			-1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n;
		
		$Python_TFFM $HeatmapDNAshape -c $output/$name/${name}_fo_classifier.pkl -2
		inkscape -z -e $output/$name/${name}_ShapeheatmapTFFM.png -w 1250 -h 1000 $output/$name/${name}_fo_classifier.pkl.svg
		cp $matrix  $output/$name/usedTFFM_fo.xml
	fi
	if [[ $matrix == *".pfm"* ]]; then  # PWM scores
		echo "WIP"
		awk 'NR==1{print ">MA0000.1 "$4}' $matrix > $output/$name/JASPAR.pfm
		awk 'NR!=1{print $0}NR==2{print "[","[","[","["}END{print "]","]","]","]"}' $matrix | awk '
		{ 
			for (i=1; i<=NF; i++)  {
				if($i ~ /^[0-9]+$/){
					if ($i ~ /^[0-9]{1}$/){a[NR,i] = "  "$i}
					if ($i ~ /^[0-9]{2}$/){a[NR,i] = " "$i}
					if ($i ~ /^[0-9]{3}$/){a[NR,i] = $i}
				} else{ if ($i == "[") {a[NR,i] = " "$i} else {a[NR,i] = $i}}
			}
		}
		NF>p { p = NF }
		END {
			for(j=1; j<=p; j++) {
				str=a[1,j]
				for(i=2; i<=NR; i++){
					str=str" "a[i,j];
				}
				print str
			}
		}' >> $output/$name/JASPAR.pfm
		echo "Training a PSSM + DNA shape classifier.";
		$Python_TFFM $ComputeDNAshaped trainPSSM -f $output/$name/JASPAR.pfm \
			-i $output/$name/training/foreground_training.fas -I $output/$name/training/foreground_training.bed \
			-b $output/$name/training/background_training.fas -B $output/$name/training/background_training.bed \
			-o $output/$name/${name}_PSSM_classifier \
			-1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n;
		$Python_TFFM $HeatmapDNAshape -c $output/$name/${name}_PSSM_classifier.pkl -2
		inkscape -z -e $output/$name/${name}_ShapeheatmapPSSM.png -w 1250 -h 1000 $output/$name/${name}_PSSM_classifier.pkl.svg
	fi
}


generate_weighted_fasta(){

while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-f)
			local file=$2
			echo "analyzed bed file is: ${2}";shift 2;;
        -c)
			local colcov=$2
			echo "column with corresponding coverage is: ${2}";shift 2;;
		-s)
			local size=$2
			echo "number of seqs (that will be used by KMAC) for the training set is $2"; shift 2;;
        -o)
			local out_dir=$2
			echo "output directory is: ${2}";shift 2;;
        -p)
			local prefix=$2
			echo "prefix is: ${2}";shift 2;;
		-g)
			local genome=$2
			echo "reference genome to extract seq from is:$2"; shift 2;;
		-h)
			usage generate_weighted_fasta; return;;
		--help)
			usage generate_weighted_fasta; return;;
		*)
			echo "Error in arguments"
			echo $1; usage; return;;
	esac
done

local Errors=0
if [ -z $file ]; then echo "ERROR: -f argument needed"; Errors+=1; fi
if [ -z $size ]; then echo "ERROR: -s argument needed"; Errors+=1; fi
if [ -z $out_dir ]; then echo "ERROR: -o argument needed"; Errors+=1; fi
if [ -z $prefix ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $colcov ]; then echo "-c argument nor used, sequence will not be weigthed in KMAC analysis"; local colcov=0; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas" ; fi
if [ $Errors -gt 0 ]; then usage generate_weighted_fasta; return 1; fi

mkdir -p $out_dir/$prefix
echo -e "\ngenerating weighted fasta"
head $file
if [ $colcov == 0 ]; then
		cat $file | head -n $size > $out_dir/$prefix/${prefix}.weighted.bed
		bedtools getfasta -fi $genome -bed $out_dir/$prefix/${prefix}.weighted.bed -fo $out_dir/$prefix/tmp.fas -fo $out_dir/$prefix/${prefix}.weighted.fas
	else
		echo "WARNING THIS WEIGHTING OF SEQUENCE HAS TO BE TRIPLE CHECKED"
	    cat $file | cut -f 1,2,3,$colcov \
	| sort -k4,4rn | head -n $size > $out_dir/$prefix/${prefix}.weighted.bed
	bedtools getfasta -fi $genome -bed $out_dir/$prefix/${prefix}.weighted.bed -fo $out_dir/$prefix/tmp.fas
	cat $out_dir/$prefix/tmp.fas | awk -F ":" '{print ">"$3,":"$4,$1}' \
		| sed 's/^> : //g' | sed 's/ >/ /g' \
		| sed 's/ :/:/g' > $out_dir/$prefix/${prefix}.weighted.fas
fi


}

compute_kmer(){

while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-g)
			local genome=$2
			echo "genome to use set to: ${2}";shift 2;;
		-anf)
			local annotation=$2
			echo "genome annotation file set to: ${2}";shift 2;;
		-b1)
			local bed1=$2
			echo "bed file for positive sequences: ${2}";shift 2;;
        -b2)
			local bed2=$2
			echo "bed file negative sequences (or another experimentale condition in a comparative analysis: ${2}";shift 2;;
		-c1)
			local colWeight1=$2
			echo "column number in the input bed file (-b1) that gives sequences weight (RPKM) (use 0 if no weight): ${2}";shift 2;;
		-c2)
			local colWeight2=$2
			echo "when the second bed IS NOT a negative set:column number, in the input bed file (e.g. output by initial_comp), that gives sequences weight (RPKM): ${2}";shift 2;;
		-n1)
			local name1=$2
			echo "name of first condition (usualy the positive set): $2"; shift 2;;
		-n2)
			local name2=$2
			echo "name of second condition (usualy the negative set): $2"; shift 2;;
		-pns)
			local path2negSets=$2
			echo "-pns is used if b2 not provided. Path to the dir where negative sets corresponding to condition 'name 1' bed files are stored: $2"; shift 2;;
		-ls)
			local learning_size=$2
			echo "size for learning set: ${2}";shift 2;;
		-kwi)
			local kmerWindow=$2
			echo "see KMAC manual: $2"; shift 2;;
		-kmi)
			local kmerMin=$2
			echo "minumum kmer length is: $2"; shift 2;;
		-kma)
			local kmerMax=$2
			echo "minumum kmer length is: $2"; shift 2;;
		-ktop)
			local ktop=$2
			echo "k_top is: $2"; shift 2;;
		-wei)
			local seqWeight=$2
			echo "the sequence weight type: $2 (0: no weighting, all weights=1.0; 1: as
 in the input flie; 2: square root; 3: natural logarithm)"; shift 2;;
        -o)
			local out_dir=$2
			echo "output directory is: ${2}";shift 2;;
		-r)
			local RUN=$2
			echo "run name: ${2}";shift 2;;
		-h)
			usage compute_kmer; return;;
		--help)
			usage compute_kmer; return;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_kmer; return;;
	esac
done

local Errors=0
if [ -z "${genome+x}" ]; then echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas" ; fi
if [ -z $bed1 ]; then echo "ERROR: -b1 argument needed"; Errors+=1; fi

if [ -z "${bed2+x}" ] && [ -z "${path2negSets+x}" ];  
then echo "ERROR: no bed2 (-b 2) and no negative set (-pns) are provided. One of this option is needed"; Errors+=1; fi

if [ -n "${bed2+x}" ]; then if [ -z "${colWeight2+x}" ] || [ -z $name2 ] ; then echo "ERROR: -c2 and -n2 are both needed to treat the second bed"; Errors+=1; fi; fi


if [ -z $colWeight1 ]; then echo "ERROR: -c1 argument needed"; Errors+=1; fi
if [ -z $name1 ]; then echo "ERROR: -n1 argument needed"; Errors+=1; fi
if [ -z $learning_size ]; then echo "ERROR: -ls argument needed"; Errors+=1; fi
if [ -z $kmerWindow ]; then echo "ERROR: -kwi argument needed"; Errors+=1; fi
if [ -z $kmerMin ]; then echo "ERROR: -kmi argument needed"; Errors+=1; fi
if [ -z $kmerMax ]; then echo "ERROR: -kma argument needed"; Errors+=1; fi
if [ -z $out_dir ]; then echo "ERROR: -o argument needed"; Errors+=1; fi
if [ -z $RUN ]; then echo "ERROR: -r argument needed"; Errors+=1; fi

if [ -z "${seqWeight+x}" ]; then echo "-wei argument not used, assuming no weighting"; local seqWeight=0 ; fi
if [ -z "${ktop+x}" ]; then echo "-kmerTop argument not used, assuming default 10"; local ktop=10 ; fi
if [ -z "${annotation+x}" ]; then echo "-anf argument not used, assuming it is tair10.bed"; local annotation=/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.bed ; fi


if [ $Errors -gt 0 ]; then echo "error in compute_kmer arguments (usage has yet to come)"; return 1; fi

# generate input files:
if [ -z "${bed2+x}" ] #to separate case conditionX_vs_NS and conditionX_vs_conditionY
	then 
	# case with NS
	local name2=${name1}_NS
	local tmpPath=$out_dir/Infiles/$name2
	mkdir -p $tmpPath
	compute_NS -p $bed1 -n $name1 -g $genome -od $tmpPath -anf $annotation -nb 10 -gc 0.05 --timeout -s 256
	cat $tmpPath/*bed > $tmpPath/neg_train_set.bed
	bedtools getfasta -fi $genome -fo $tmpPath/neg_train_set.fas -bed $tmpPath/neg_train_set.bed
	fastaNEG=$tmpPath/neg_train_set.fas
	
	#generate sequence for positive set
	generate_weighted_fasta -f $bed1 -c $colWeight1 -s $learning_size -o $out_dir/Infiles/ -p $name1 -g $genome
	
	else
	# case without NS (needs a second bed!)
	local learning_size2=$(wc -l $bed2 | cut -d " " -f 1)
	generate_weighted_fasta -f $bed2 -c $colWeight2 -s $learning_size2 -o $out_dir/Infiles -p $name2 -g $genome
	fastaNEG=$out_dir/Infiles/$name2/$name2.weighted.fas
	
	local learning_size1=$(wc -l $bed1 | cut -d " " -f 1)
	generate_weighted_fasta -f $bed1 -c $colWeight1 -s $learning_size1 -o $out_dir/Infiles -p $name1 -g $genome
	fi

fastaPOS=$out_dir/Infiles/$name1/$name1.weighted.fas

# run KMAC
echo "This function uses the KMAC program (GEM suite) to identify kmer sets that are enriched in a group of sequences compare to another group. KMAC pub: https://pubmed.ncbi.nlm.nih.gov/29654070/"

#local RUN=${kmerWindow}_${kmerMin}_${kmerMax}_${seqWeight}
#mkdir -p ${out_dir}

kmac(){
java -Xmx8G -jar $GEM KMAC \
       --pos_seq $fastaPOS \
       --k_win $kmerWindow \
       --k_min $kmerMin \
       --k_max $kmerMax \
	   --swt $seqWeight \
       --k_top $ktop \
       --out_name $RUN \
       --print_aligned_seqs ON \
       --gc -1 \
       --neg_seq $fastaNEG
       
}
kmac

local tmpKMAC=$out_dir/KMAC/${name1}_vs_${name2} #not that when running the pipeline for atlas the out_dir here is /nobackup/rb261841/, so it need to be change in the X_kmer.sh because the path change when I copy the results from the nobackup to the team workspace...

mkdir -p $tmpKMAC
local tmpPWD=$(pwd)

sed -i "s~$tmpPWD~$tmpKMAC~g" ${RUN}_outputs/${RUN}.ksm_list.txt
head ${RUN}_outputs/${RUN}.ksm_list.txt
if [ -d "$tmpKMAC/${RUN}_outputs" ]; then
  rm -Rf $tmpKMAC/${RUN}_outputs
fi
mv ${RUN}_outputs $tmpKMAC/.
printf "\nKMAC has finished\n\n"
}


compute_kmer_OLD(){

while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-g)
			local genome=$2
			echo "genome to use set to: ${2}";shift 2;;
		-b1)
			local bed1=$2
			echo "bed file for positive sequences: ${2}";shift 2;;
        -b2)
			local bed2=$2
			echo "bed file negative sequences (or another experimentale condition in a comparative analysis: ${2}";shift 2;;
		-c1)
			local colWeight1=$2
			echo "column number in the input bed file (-b1) that gives sequences weight (RPKM) (use 0 if no weight): ${2}";shift 2;;
		-c2)
			local colWeight2=$2
			echo "when the second bed IS NOT a negative set:column number, in the input bed file (e.g. output by initial_comp), that gives sequences weight (RPKM): ${2}";shift 2;;
		-n1)
			local name1=$2
			echo "name of first condition (usualy the positive set): $2"; shift 2;;
		-n2)
			local name2=$2
			echo "name of second condition (usualy the negative set): $2"; shift 2;;
		-pns)
			local path2negSets=$2
			echo "-pns is used if b2 not provided. Path to the dir where negative sets corresponding to condition 'name 1' bed files are stored: $2"; shift 2;;
		-ls)
			local learning_size=$2
			echo "size for learning set: ${2}";shift 2;;
		-kwi)
			local kmerWindow=$2
			echo "see KMAC manual: $2"; shift 2;;
		-kmi)
			local kmerMin=$2
			echo "minumum kmer length is: $2"; shift 2;;
		-kma)
			local kmerMax=$2
			echo "minumum kmer length is: $2"; shift 2;;
		-ktop)
			local ktop=$2
			echo "k_top is: $2"; shift 2;;
		-wei)
			local seqWeight=$2
			echo "the sequence weight type: $2 (0: no weighting, all weights=1.0; 1: as in the input flie; 2: square root; 3: natural logarithm)"; shift 2;;
        -o)
			local out_dir=$2
			echo "output directory is: ${2}";shift 2;;
		-r)
			local RUN=$2
			echo "run name: ${2}";shift 2;;
		-h)
			usage compute_kmer; return;;
		--help)
			usage compute_kmer; return;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_kmer; return;;
	esac
done

local Errors=0
if [ -z "${genome+x}" ]; then echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas" ; fi
if [ -z $bed1 ]; then echo "ERROR: -b1 argument needed"; Errors+=1; fi

if [ -z "${bed2+x}" ] && [ -z "${path2negSets+x}" ];  
then echo "ERROR: no bed2 (-b 2) and no negative set (-pns) are provided. One of this option is needed"; Errors+=1; fi

if [ -n "${bed2+x}" ]; then if [ -z "${colWeight2+x}" ] || [ -z $name2 ] ; then echo "ERROR: -c2 and -n2 are both needed to treat the second bed"; Errors+=1; fi; fi


if [ -z $colWeight1 ]; then echo "ERROR: -c1 argument needed"; Errors+=1; fi
if [ -z $name1 ]; then echo "ERROR: -n1 argument needed"; Errors+=1; fi
if [ -z $learning_size ]; then echo "ERROR: -ls argument needed"; Errors+=1; fi
if [ -z $kmerWindow ]; then echo "ERROR: -kwi argument needed"; Errors+=1; fi
if [ -z $kmerMin ]; then echo "ERROR: -kmi argument needed"; Errors+=1; fi
if [ -z $kmerMax ]; then echo "ERROR: -kma argument needed"; Errors+=1; fi
if [ -z $out_dir ]; then echo "ERROR: -o argument needed"; Errors+=1; fi
if [ -z $RUN ]; then echo "ERROR: -r argument needed"; Errors+=1; fi

if [ -z "${seqWeight+x}" ]; then echo "-wei argument not used, assuming no weighting"; local seqWeight=0 ; fi
if [ -z "${ktop+x}" ]; then echo "-kmerTop argument not used, assuming default 10"; local ktop=5 ; fi

if [ $Errors -gt 0 ]; then echo "error in compute_kmer arguments (usage has yet to come)"; return 1; fi

# generate input files:
if [ -z "${bed2+x}" ] #to separate case conditionX_vs_NS and conditionX_vs_conditionY
	then 
	# case with NS
	echo "by default the kmer enrichment will be done against part of the negative set (N-1) and the last neg set will be used for the ROC curve"
	local name2=${name1}_NS
	numNS=$(ls -ltrh ${path2negSets}/${name1}_*_neg.bed | wc -l)
	local tmpPath=$out_dir/Infiles/$name2
	mkdir -p $tmpPath
	#below I take NS from 2 to N, the number 1 is kept for the ROCs (because PWM and TFFM are run against NS number 1)
	cat ${path2negSets}/${name1}_[2-$numNS]_neg.bed > $tmpPath/neg_train_set.bed
	bedtools getfasta -fi $genome -fo $tmpPath/neg_train_set.fas -bed $tmpPath/neg_train_set.bed
	fastaNEG=$tmpPath/neg_train_set.fas
	# subset the 1st negative so that it has the same number of sequences as the pos test set:
	cat ${path2negSets}/${name1}_1_neg.bed | shuf | sed "1,${learning_size}d" > ${path2negSets}/${name1}_neg_testing_set.bed
	
	generate_weighted_fasta -f $bed1 -c $colWeight1 -s $learning_size -o $out_dir/Infiles/ -p $name1 -g $genome
	else
	# case without NS (needs a second bed!)
	local learning_size2=$(wc -l $bed2 | cut -d " " -f 1)
	generate_weighted_fasta -f $bed2 -c $colWeight2 -s $learning_size2 -o $out_dir/Infiles -p $name2 -g $genome
	fastaNEG=$out_dir/Infiles/$name2/$name2.weighted.fas
	
	local learning_size1=$(wc -l $bed1 | cut -d " " -f 1)
	generate_weighted_fasta -f $bed1 -c $colWeight1 -s $learning_size1 -o $out_dir/Infiles -p $name1 -g $genome
	fi

fastaPOS=$out_dir/Infiles/$name1/$name1.weighted.fas

# run KMAC
echo "This function uses the KMAC program (GEM suite) to identify kmer sets that are enriched in a group of sequences compare to another group. KMAC pub: https://pubmed.ncbi.nlm.nih.gov/29654070/"

#local RUN=${kmerWindow}_${kmerMin}_${kmerMax}_${seqWeight}
#mkdir -p ${out_dir}

kmac(){
java -Xmx8G -jar $GEM KMAC \
       --pos_seq $fastaPOS \
       --k_win $kmerWindow \
       --k_min $kmerMin \
       --k_max $kmerMax \
	   --swt $seqWeight \
       --k_top $ktop \
       --out_name $RUN \
       --print_aligned_seqs ON \
       --gc -1 \
       --neg_seq $fastaNEG
       
}
kmac

local tmpKMAC=$out_dir/KMAC/${name1}_vs_${name2} #not that when running the pipeline for atlas the out_dir here is /nobackup/rb261841/, so it need to be change in the X_kmer.sh because the path change when I copy the results from the nobackup to the team workspace...

mkdir -p $tmpKMAC
local tmpPWD=$(pwd)

sed -i "s~$tmpPWD~$tmpKMAC~g" ${RUN}_outputs/${RUN}.ksm_list.txt

mv ${RUN}_outputs $tmpKMAC/.
printf "\nKMAC has finished\n\n"
}


compute_ROCS(){
# compute_ROCS -p -ns -m -n -g -od -pc
local pocc=false; local colors=('#40A5C7' '#F9626E' '#F0875A' '#307C95' '#BB4A52' '#B46544'); local number_of_ksm=1; local ksm_score_meth="BEST"; local dimer="Mono"; local offset_left=0; local offset_right=0
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=("${!2}")
			echo "peaks files set to: ${peaks[@]}";shift 2;;
		-ns)
			local negative_sets=("${!2}")
			echo "negative files set to: ${negative_sets[@]}";shift 2;;
		-m)
			local matrices=("${!2}")
			echo "matrix files set to: ${matrices[@]}";shift 2;;
		-n)
			local names=("${!2}")
			echo "names associated set to: ${names[@]}";shift 2;;
		-g)
			local genome=$2
			echo "fasta of the genome set to: ${2}";shift 2;;
		-od)
			local results=$2
			echo "output directory set to: ${2}";shift 2;;
		-pc)
			local pocc=true
			echo "Pocc mode activated, pfm will be used to compute PWM score and Pocc";shift 1;;
		-d)
			local dimer=$2
			echo "PFM will be tested as a dimer too: conformations: ${2}"; shift 2;;
		-ol)
			local offset_left=$2
			echo "offset on the left set to: ${2}";shift 2;;
		-or)
			local offset_right=$2
			echo "offset on the right set to: ${2}";shift 2;;
		-nKSM)
			local number_of_ksm=$2
			echo "number of kmer set motif to use in the ksm search: ${2}";shift 2;;
		-sKSM)
			local ksm_score_meth=$2
			echo "kmer score method (SUM, BEST or MEAN): ${2}";shift 2;;
		-color)
			local colors=("${!2}")
			echo "colors set to: ${colors[@]}";shift 2;;
		-h)
			usage compute_ROCS; return;;
		--help)
			usage compute_ROCS; return;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_ROCS; return;;
	esac
done
local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed or needs to be a LIST"; Errors+=1; fi
if [ -z $negative_sets ]; then echo "ERROR: -ns argument needed or needs to be a LIST"; Errors+=1; fi
if [ -z $matrices ]; then echo "ERROR: -m argument needed or needs to be a LIST"; Errors+=1; fi
if [ -z $names ]; then echo "ERROR: -n argument needed or needs to be a LIST"; Errors+=1; fi
if [ -z $results ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas" ; fi
if [ -z $colors ]; then echo "-color argument not used, using default or needs to be a LIST"; fi
#two options for the kmer search
if [ -z $number_of_ksm ]; then echo "-nKSM argument not used, assuming you are using the first kmer set motif (set to 0 if you want all kmer set motifs (not recommended))" ; fi
if [ -z $ksm_score_meth ]; then echo "-sKSM argument not used, using default: BEST" ; fi

if [ ${#names[@]} -ne ${#peaks[@]} ] || [ ${#names[@]} -ne ${#matrices[@]} ] || [ ${#names[@]} -ne ${#negative_sets[@]} ] || [ ${#names[@]} -ge ${#colors[@]} ]; then
	if [ ${#peaks[@]} -eq ${#negative_sets[@]} ] && [ ${#peaks[@]} -ne ${#matrices[@]} ] && [ ${#names[@]} -eq ${#peaks[@]} ] && [ ${#colors[@]} -ge ${#matrices[@]} ] ; then
		for ((i=1;i<${#matrices[@]};i++)); do
			peaks+=("${peaks[0]}")
			negative_sets+=("${negative_sets[0]}")
			names+=("${names[0]}")
		done
	else
		if [ ${#peaks[@]} -eq ${#negative_sets[@]} ] && [ ${#peaks[@]} -ne ${#matrices[@]} ] && [ ${#names[@]} -eq ${#matrices[@]} ] && [ ${#colors[@]} -ge ${#matrices[@]} ] ; then
			for ((i=1;i<${#matrices[@]};i++)); do
			peaks+=("${peaks[0]}")
			negative_sets+=("${negative_sets[0]}")
		done
		else
			echo "ERROR: -n & -m have to be lists of same length"; Errors+=1
		fi
	fi
fi

if [ $Errors -gt 0 ]; then usage compute_ROCS; return 1; fi

for ((i=0;i<${#colors[@]};i++)); do
	[[ ${colors[$i]} =~ ^#.* ]] || colors[$i]="#${colors[$i]}"
done

# pocc_pfm=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/compute_POcc.py
# tffmscores=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/get_best_score_tffm.py
# scores_prog=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/scores.py
# plot_ROCS_prog=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/plots_ROCS_multiple.py

local scores=()
local endnames=()
i=0
mkdir -p -m 774 ${results}/scores

# for peak in ${peaks[@]}
for ((i=0;i<${#peaks[@]};i++))
do
	local negative_set=${negative_sets[i]}
	local name=${names[i]// /_}
	local matrice=${matrices[i]}
	local peak=${peaks[i]}
    echo -e "\tpeaks $peak"
    echo -e "\tNegset $negative_set"
    echo -e "\tname $name"
    echo -e "\tmatrice $matrice"
	if [[ $peak != *".fa"* ]]; then
		if [[ $matrice == *".xml"* ]]; then
			awk -v OFS="\t" '{print $1,$2-2,$3}' $peak > $results/pos_reworked.bed
			awk -v OFS="\t" '{print $1,$2-2,$3}' $negative_set > $results/neg_reworked.bed
			local peak=$results/pos_reworked.bed
			local negative_set=$results/neg_reworked.bed
		fi
		if [[ $matrice == *".pkl"* ]]; then
			awk -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' $peak > $results/pos_reworked.bed
			awk -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' $negative_set > $results/neg_reworked.bed
			local peak=$results/pos_reworked.bed
			local negative_set=$results/neg_reworked.bed
		fi
		bedtools getfasta -fi $genome -bed $peak -fo $results/pos_set.fa
		local peak=$results/pos_set.fa
		bedtools getfasta -fi $genome -bed $negative_set -fo $results/neg_set.fa
		local negative_set=$results/neg_set.fa
	fi

	if [[ $matrice == *".xml"* ]]; then # TFFM scores computation
		$Python_TFFM $tffmscores -o ${results}/scores/tffm_scores_pos.tsv -pos $peak -t ${matrice}
		$Python_TFFM $tffmscores -o ${results}/scores/tffm_scores_neg.tsv -pos $negative_set -t ${matrice}
		paste <(awk '$1!="None"{print $8;next}{print "0.0"}' ${results}/scores/tffm_scores_pos.tsv) <(awk '$1!="None"{print $8;next}{print "0.0"}' ${results}/scores/tffm_scores_neg.tsv) >"${results}/scores/TFFM_${name}_scores.tsv"
		scores+=("${results}/scores/TFFM_${name}_scores.tsv")
		endnames+=("${name}_TFFM")
		echo "there"
		$R_36 $F_score $results/scores TFFM_${name}_scores.tsv TFFM $name
		echo "passed"
	fi
	if [[ $matrice == *".txt"* ]]; then # K-mer scores computation
		local tmpS=${results}/scores
		
		if [ $number_of_ksm == 0 ]; then cp $matrice ${results}/scores/tmp_ksm_list; else cat $matrice | head -n $number_of_ksm > ${results}/scores/tmp_ksm_list; fi
		
		cat $peak | tr ':' '_' > ${results}/scores/tmpPos.fsa # coz gem does not like ":" ...
		java -Xmx8G -jar $GEM KSM --fasta ${results}/scores/tmpPos.fsa --ksm ${results}/scores/tmp_ksm_list --out $tmpS/ksm.scan.pos
		cat $negative_set | tr ':' '_' > ${results}/scores/tmpNeg.fsa
		java -Xmx8G -jar $GEM KSM --fasta ${results}/scores/tmpNeg.fsa --ksm ${results}/scores/tmp_ksm_list --out $tmpS/ksm.scan.neg
		rm ${results}/scores/tmpNeg.fsa ${results}/scores/tmpPos.fsa ${results}/scores/tmp_ksm_list
		
		python $parse_KSM_scores -p $tmpS/ksm.scan.pos.motifInstances.txt -n $tmpS/ksm.scan.neg.motifInstances.txt -o $tmpS/${name}_KSM_${ksm_score_meth}_scores.tsv -m $ksm_score_meth
		
		scores+=("${results}/scores/${name}_KSM_${ksm_score_meth}_scores.tsv")
		endnames+=("${name}_KSM")
		$R_36 $F_score $results/scores ${name}_KSM_${ksm_score_meth}_scores.tsv KSM $name
		
	fi
	if [[ $matrice == *".pkl"* ]]; then # DNAshape
		if [[ $matrice == *"_fo_classifier"* ]]; then # DNAshape done with first order TFFM
			echo "Applying the trained detailed TFFM + DNA shape classifier on foreground sequences.";
			$Python_TFFM $ComputeDNAshaped applyTFFM -T $(dirname $matrice)/usedTFFM_fo.xml -i $peak -I $results/pos_reworked.bed -c $matrice -o ${results}/scores/DNAshape_pred_pos.txt -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0.000001
			echo "Applying the trained detailed TFFM + DNA shape classifier on background sequences.";
			$Python_TFFM $ComputeDNAshaped applyTFFM -T $(dirname $matrice)/usedTFFM_fo.xml -i $negative_set -I $results/neg_reworked.bed -c $matrice -o ${results}/scores/DNAshape_pred_neg.txt -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0.000001
			paste <(awk 'NR!=1{print $6}' ${results}/scores/DNAshape_pred_pos.txt | sort -nr ) <(awk 'NR!=1{print $6}' ${results}/scores/DNAshape_pred_neg.txt | sort -nr ) >"${results}/scores/tab_${name}_TFFMshape.tsv"
			scores+=("${results}/scores/tab_${name}_TFFMshape.tsv")
			endnames+=("${name}_TFFMShape")
			$R_36 $F_score $results/scores ${results}/scores/tab_${name}_TFFMshape.tsv TFFMShape $name
		fi
		if [[ $matrice == *"_PSSM_classifier"* ]]; then # DNAshape done with PSSM/PFM
			echo "Applying the trained detailed PSSM + DNA shape classifier on foreground sequences.";
			$Python_TFFM $ComputeDNAshaped applyPSSM -f $(dirname $matrice)/JASPAR.pfm -i $peak -I $results/pos_reworked.bed -c $matrice -o ${results}/scores/DNAshapePSSM_pred_pos.txt -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0.000001
			echo "Applying the trained detailed PSSM + DNA shape classifier on background sequences.";
			$Python_TFFM $ComputeDNAshaped applyPSSM -f $(dirname $matrice)/JASPAR.pfm -i $negative_set -I $results/neg_reworked.bed -c $matrice -o ${results}/scores/DNAshapePSSM_pred_neg.txt -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0.000001
			paste <(awk 'NR!=1{print $6}' ${results}/scores/DNAshapePSSM_pred_pos.txt | sort -nr ) <(awk 'NR!=1{print $6}' ${results}/scores/DNAshapePSSM_pred_neg.txt | sort -nr ) >"${results}/scores/tab_${name}_PSSMshape.tsv"
			scores+=("${results}/scores/tab_${name}_PSSMshape.tsv")
			endnames+=("${name}_PSSMShape")
			$R_36 $F_score $results/scores ${results}/scores/tab_${name}_PSSMshape.tsv PSSMShape $name
		fi
	fi
	if [[ $matrice == *".pfm"* ]]; then  # PWM scores
		python $scores_prog -m ${matrice} -f $peak -o ${results}/scores/
		python $scores_prog -m ${matrice} -f $negative_set -o ${results}/scores/
		paste <(sort -u -k1,1 -k8,8nr ${results}/scores/$(basename $peak).scores | sort -u -k1,1 | awk '{print $8}'  | sort -nr )  <(sort -u -k1,1 -k8,8nr ${results}/scores/$(basename $negative_set).scores | sort -u -k1,1 | awk '{print $8}' | sort -nr ) | awk '{print $0}' >"${results}/scores/tab_${name}.tsv"
		scores+=("${results}/scores/tab_${name}.tsv")
		endnames+=("${name}_PWM")
		$R_36 $F_score $results/scores tab_${name}.tsv PFM $name
		
		if $pocc ; then # PWM Pocc
			python $pocc_pfm -s ${results}/scores/$(basename $peak).scores -o ${results}/scores/Pocc_pos.pocc
			python $pocc_pfm -s ${results}/scores/$(basename $negative_set).scores -o ${results}/scores/Pocc_neg.pocc
			paste <( sort -nr ${results}/scores/Pocc_pos.pocc )  <( sort -nr ${results}/scores/Pocc_neg.pocc ) > "${results}/scores/${name}_Pocc_scores.tsv"
			scores+=("${results}/scores/${name}_Pocc_scores.tsv")
			endnames+=("${name}_Pocc")
		fi
		
		if [ ${dimer} != "Mono" ]; then
			local dist=${dimer:2}
			local conf=${dimer:0:2}
			local dimer_name="${name}_$dimer"
			bash ${dimer_builder} ${matrice} ${dimer_name} $dist $conf $offset_left $offset_right > ${results}/scores/$(basename $matrice)_${dimer}.pfm
			
			python $scores_prog -m ${results}/scores/$(basename $matrice)_${dimer}.pfm -f $peak -o ${results}/scores/
			python $scores_prog -m ${results}/scores/$(basename $matrice)_${dimer}.pfm -f $negative_set -o ${results}/scores/
			paste <(sort -u -k1,1 -k8,8nr ${results}/scores/$(basename $peak).scores | sort -u -k1,1 | awk '{print $8}'  | sort -nr )  <(sort -u -k1,1 -k8,8nr ${results}/scores/$(basename $negative_set).scores | sort -u -k1,1 | awk '{print $8}' | sort -nr ) | awk '{print $0}' >"${results}/scores/tab_${name}_${dimer}.tsv"
			scores+=("${results}/scores/tab_${name}_${dimer}.tsv")
			endnames+=("${name}_${dimer}")
			
			
		fi
		
	fi
	
done

# local scores=`join_by " " "${scores[@]}"`
# local endnames=`join_by " " "${endnames[@]}"`
# local list_colors=`join_by " " "${colors[@]}"`
echo ${colors[@]}
python $plot_ROCS_prog -s ${scores[@]} -n ${endnames[@]} -o ${results} -of ROC.svg -c ${colors[@]}
inkscape -z -e ${results}/ROC.png -w 1000 -h 1000 ${results}/ROC.svg
}

compute_space(){
local Errors=0; local matrix_type="ASYMMETRIC"; local thresholds_dir="null"; local thresholds=(0); local maxy=0; local maxSpace=50; local minSpace=0; local offset_left=0; local offset_right=0; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=("${!2}")
			echo "peaks files set to: ${peaks[@]}";shift 2;;
		-m)
			local matrices=("${!2}")
			echo "matrix files set to: ${matrices[@]}";shift 2;;
		-n)
			local names=("${!2}")
			echo "names associated set to: ${names[@]}";shift 2;;
		-th)
			local thresholds=("${!2}")
			echo "thresholds set to: ${thresholds[@]}";shift 2;;
		-thf)
			local thresholds_dir=$2
			echo "thresholds directory set to: ${2}"; shift 2;;
		-od)
			local results=$2
			echo "output directory set to: ${2}";shift 2;;
		-maxy)
			local maxy=$2
			echo "maximum enrichment to display set to: ${2}";shift 2;;
        -maxs)
			local maxSpace=$2
			echo "maximum spacing to compute set to: ${2}";shift 2;;
        -mins)
			local minSpace=$2
			echo "minimum spacing to compute set to: ${2}";shift 2;;
        -ol)
			local offset_left=$2
			echo "offset on the left set to: ${2}";shift 2;;
        -or)
			local offset_right=$2
			echo "offset on the right set to: ${2}";shift 2;;
        -g)
			local genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
		-sym)
			local matrix_type="SYMMETRIC"; shift 1;;
		-h)
			echo "you asked for help" ;Errors+=1; shift 1;;
		--help)
			echo "you asked for help" ;Errors+=1; shift 1;;
		*)
			echo "Error in arguments"
			echo $1; Errors+=1; break;;
	esac
done


if [ -z $peaks ]; then echo "ERROR: -p argument needed or may not be a LIST"; Errors+=1; fi
if [ -z $matrices ]; then echo "ERROR: -m argument needed  or may not be a LIST"; Errors+=1; fi
if [ -z $names ]; then echo "ERROR: -n argument needed  or may not be a LIST"; Errors+=1; fi
# if [ -z $thresholds ] || []; then echo "ERROR: -th argument needed  or may not be a LIST"; Errors+=1; fi
if [ -z $results ]; then echo "ERROR: -od argument needed"; Errors+=1; fi

if [ ${thresholds[0]} -eq 0 ] && [ ${thresholds_dir} == "null" ]; then echo "ERROR either -th or -thd arguments needed; -th arguments needs to be a list"; Errors+=1; fi
if [ ${thresholds_dir} != "null" ] && [ ! -d ${thresholds_dir} ]; then echo "ERROR directory for thresholds acquisition does not exist, please check your argument -thd"; Errors+=1; fi

if [ ${matrix_type} == "SYMMETRIC" ]; then echo "palindromic mode enabled"; fi

if [[ ${#peaks[@]} -ne ${#matrices[@]} ]] || [[ ${#peaks[@]} -ne ${#names[@]} ]]; then
	echo "ERROR: LISTs for -p, -n and -m arguments should have the same length"; Errors+=1
fi

if [ -z $maxy ]; then echo "-maxy argument not used, no limits"; fi
if [ -z $maxSpace ]; then echo "-maxs argument not used, using 50"; fi
if [ -z $minSpace ]; then echo "-mins argument not used, using 0"; fi
if [ -z $offset_left ]; then echo "-ol argument not used, no offset"; fi
if [ -z $offset_right ]; then echo "-or argument not used, no offset"; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas" ; fi

if [ $Errors -gt 0 ]; then usage compute_space; return 1; fi

for ((i=0;i<${#peaks[@]};i++))
do
	
	local name=${names[i]// /_}
	local matrice=${matrices[i]}
	local peak=${peaks[i]}
	local length_mat=0
	mkdir -p -m 774 $results/${name} $results/${name}/scores
	
	if [[ $peak == *"score"* ]];then
		local pos_file=$peak # no preparation needed
	elif [[ $peak != *".fa"* ]]; then
		if [ ! -f $results/${name}/${name}_pos_set.fa ];then
			bedtools getfasta -fi $genome -bed $peak -fo $results/${name}/${name}_pos_set.fa # need to tranform in fasta
		fi
		local peak=$results/${name}/${name}_pos_set.fa
	fi
	if [[ $peak == *".fa"* ]]; then # need to compute scores
		if [[ $matrice == *".xml"* ]]; then
			if [ ! -f ${results}/${name}/scores/${name}_tffm_scores_pos.tsv ];then
				$Python_TFFM $tffm_all_scores -o ${results}/${name}/scores/${name}_tffm_scores_pos.tsv -pos $peak -t ${matrice}
			fi
			if [ ${thresholds[0]} -eq 0 ]; then
				if [ -f ${thresholds_dir}/ROCs/${name}/scores/Fscores_${name}_TFFM.txt ]; then thresholds=($(head -n1 ${thresholds_dir}/ROCs/${name}/scores/Fscores_${name}_TFFM.txt | tr ' ' "\n")); fi
				if [ -f ${thresholds_dir}/${name}/scores/Fscores_${name}_TFFM.txt ]; then thresholds=($(head -n1 ${thresholds_dir}/${name}/scores/Fscores_${name}_TFFM.txt | tr ' ' "\n")); fi
				if [ -f ${thresholds_dir}/scores/Fscores_${name}_TFFM.txt ]; then thresholds=($(head -n1 ${thresholds_dir}/scores/Fscores_${name}_TFFM.txt | tr ' ' "\n")); fi
				if [ -f ${thresholds_dir}/Fscores_${name}_TFFM.txt ]; then thresholds=($(head -n1 ${thresholds_dir}/Fscores_${name}_TFFM.txt | tr ' ' "\n")); fi
			fi
			local pos_file=${results}/${name}/scores/${name}_tffm_scores_pos.tsv
		elif [[ $matrice == *".pfm"* ]]; then
			local length_mat=$(awk -v OFS="[\t ]" '{if(NR==1){if($5=="SIMPLE"){typM="simple"} else {typM="dependency"};count=-1;next};if(typM=="simple"){count++};if(typM=="dependency"){if($1=="DEPENDENCY"){count-=2;exit};count++}}END{print count}' $matrice)
			
			if [ ! -f ${results}/${name}/scores/$(basename $peak).scores ];then
				python ${scores_prog} -m ${matrice} -f $peak -o ${results}/${name}/scores/
			fi
			
			if [ ${thresholds[0]} -eq 0 ]; then
				if [ -f ${thresholds_dir}/ROCs/${name}/scores/Fscores_${name}_PFM.txt ]; then thresholds=($(head -n1 ${thresholds_dir}/ROCs/${name}/scores/Fscores_${name}_PFM.txt | tr ' ' "\n")); fi
				if [ -f ${thresholds_dir}/${name}/scores/Fscores_${name}_PFM.txt ]; then thresholds=($(head -n1 ${thresholds_dir}/${name}/scores/Fscores_${name}_PFM.txt | tr ' ' "\n")); fi
				if [ -f ${thresholds_dir}/scores/Fscores_${name}_PFM.txt ]; then thresholds=($(head -n1 ${thresholds_dir}/scores/Fscores_${name}_PFM.txt | tr ' ' "\n")); fi
				if [ -f ${thresholds_dir}/Fscores_${name}_PFM.txt ]; then thresholds=($(head -n1 ${thresholds_dir}/Fscores_${name}_PFM.txt | tr ' ' "\n")); fi
			fi
			local pos_file=${results}/${name}/scores/$(basename $peak).scores
		fi
	fi
# 	local th=`join_by " " "${thresholds[@]}"` DEPRECATED by JL 13/01/2022 while making snakemake
	# compute spacing conformations in peaks
	python $spacing_mk -o $results/${name}/${name} -smax $maxSpace -smin $minSpace -pos $pos_file -th ${thresholds[@]} -ol $offset_left -or $offset_right -lm $length_mat
	# compute Z-score and plots
	conda run -n py38_R36 Rscript $Zscore_spacing ${results}/${name}/${name}_spacing.tsv $matrix_type ${thresholds[@]} ${results}/${name}
# 	/home/prog/R/R-4.1.0/bin/Rscript $Zscore_spacing ${results}/${name}/${name}_spacing.tsv $matrix_type ${thresholds[@]} ${results}/${name}  DEPRECATED by JL 13/01/2022 while making snakemake
done

}

spacing_2TFs () {
local Errors=0; local matrix_type="ASYMMETRIC"; local thresholds_dir="null"; local thresholds=(0); local maxy=0; local maxSpace=50; local minSpace=0; local offset_left=0; local offset_right=0; local offset_left2=0; local offset_right2=0; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
        -p)
			local peaks=$2
			echo "name of bed files with peaks coordinates: ${2}";shift 2;;
		-n)
			local name=$2
			echo "a name used as prefix for output directory and files: ${2}";shift 2;;
		-g)
			local genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
        -ma)
			local PFMa=$2
			echo "PFM file for first TF: ${2}";shift 2;;
		-mb)
			local PFMb=$2
			echo "PFM file for second TF: ${2}";shift 2;;
		-od)
			local outdir=$2
			echo "name of out directory: ${2}";shift 2;;
		-tha)
			local thresholds_a=("${!2}")
			echo "an array of three PWM score thresholds for first TF: ${2}";shift 2;;
		-thb)
			local thresholds_b=("${!2}")
			echo "an array of three PWM score thresholds for second TF: ${2}";shift 2;;
		-maxs)
			local maxSpace=$2
			echo "maximum spacing to compute set to: ${2}";shift 2;;
        -mins)
			local minSpace=$2
			echo "minimum spacing to compute set to: ${2}";shift 2;;
		-ol)
			local offset_left=$2
			echo "1st matrix offset on the left set to: ${2}";shift 2;;
        -or)
			local offset_right=$2
			echo "1st matrix offset on the right set to: ${2}";shift 2;;
		-ol2)
			local offset_left2=$2
			echo "2nd matrix offset on the left set to: ${2}";shift 2;;
        -or2)
			local offset_right2=$2
			echo "2nd matrix offset on the right set to: ${2}";shift 2;;
		-sym)
			local matrix_type="SYMMETRIC"; shift 1;;
		-h)
			usage ; return;;
		--help)
			usage ; return;;
		*)
			echo "Error in arguments"
			echo $1; usage download_SRA ; return;;
	esac
done

local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed";Errors+=1;fi
if [ -z $name ]; then echo "ERROR: -n argument needed";Errors+=1;fi
if [ -z $genome ]; then echo "ERROR: -g argument needed";Errors+=1;fi
if [ -z $PFMa ]; then echo "ERROR: -ma argument needed";Errors+=1;fi
if [ -z $PFMb ]; then echo "ERROR: -mb argument needed";Errors+=1;fi
if [ -z $outdir ]; then echo "ERROR: -od argument needed";Errors+=1;fi
if [ -z $thresholds_a ]; then echo "ERROR: -tha argument needed or needs to be a LIST";Errors+=1;fi
if [ -z $thresholds_b ]; then echo "ERROR: -thb argument needed or needs to be a LIST";Errors+=1;fi
if [ -z $maxSpace ]; then echo "-maxs argument not used, using 50"; fi
if [ -z $minSpace ]; then echo "-mins argument not used, using 0"; fi
# if [ -z $offset_left ]; then echo "-ol argument not used, no offset"; fi
# if [ -z $offset_right ]; then echo "-or argument not used, no offset"; fi
if [ $Errors -gt 0 ]; then usage spacing_2TFs ; return 1; fi


tha1=${thresholds_a[0]}; tha2=${thresholds_a[1]}; tha3=${thresholds_a[2]}
thb1=${thresholds_b[0]}; thb2=${thresholds_b[1]}; thb3=${thresholds_b[2]}

#local thresholds_a=`join_by " " "${thresholds_a[@]}"`
#local thresholds_b=`join_by " " "${thresholds_b[@]}"`

mkdir -p $outdir/scores_TFa $outdir/scores_TFb
bedtools getfasta -fi $genome -bed  $peaks > $outdir/${name}.fasta
peaks=$outdir/${name}.fasta


if [ ! -f $outdir/scores_TFa/${name}.fasta.scores ]; then $Python_TFFM $scores_prog -m $PFMa -f $peaks -o $outdir/scores_TFa; fi
if [ ! -f $outdir/scores_TFb/${name}.fasta.scores ]; then $Python_TFFM $scores_prog -m $PFMb -f $peaks -o $outdir/scores_TFb; fi


run="space_a"${tha1}_${tha2}_${tha3}"_b"${thb1}_${thb2}_${thb3}
spaceOut=$outdir/$run
mkdir -p $spaceOut
lm_a=$(cat $outdir/scores_TFa/${name}.fasta.scores | cut -f 7 | head -1)
printf "\nlength of motif $lm_a\n"
lm_b=$(cat $outdir/scores_TFb/${name}.fasta.scores |cut -f 7 | head -1)
printf "\nlength of motif $lm_b\n"

echo "get interdistance for spacing 2 TFs"
python $get_interdistances_2TF -pos $outdir/scores_TFa/${name}.fasta.scores -neg $outdir/scores_TFa/${name}.fasta.scores -pos2 $outdir/scores_TFb/${name}.fasta.scores -neg2 $outdir/scores_TFb/${name}.fasta.scores -ol $offset_left -or $offset_right -ol2 $offset_left2 -or2 $offset_right2 -th $tha1 $tha2 $tha3 -th2 $thb1 $thb2 $thb3 -smin $minSpace -smax $maxSpace -lm $lm_a -lm2 $lm_b -wi -o $spaceOut/out

sed "1d" $spaceOut/out_spacing_pos21.tsv | awk -v OFS="\t" '{print $1,$2,$4,$3,$6,$5,$8,$7}' > $spaceOut/out_spacing_pos21_converted.tsv
cat $spaceOut/out_spacing_pos12.tsv $spaceOut/out_spacing_pos21_converted.tsv |sed "1d" | sed "s/-/:/" | sed "s/:/\t/g" | sed "s/_/\t/" | awk -v OFS="\t" '{print $0,$3-$2+1}' | sed "1ichr\tStart\tEnd\tConf\tSpace\tScore1\tScore2\tmatricePosition1\tmatricePosition2\tcorrectedPosition1\tcorrectedPosition2\tSize" > $spaceOut/out_spacing_pos.tsv

/home/prog/R/R-3-5-0/bin/Rscript $Zacing_2TFs $spaceOut/out_spacing_pos.tsv $matrix_type $tha1 $tha2 $tha3 $thb1 $thb2 $thb3 $spaceOut

}

add_coverage(){

while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-b)
			local bedgraphs=("${!2}")
			echo "bedgraphs to add coverage to the table: ${bedgraphs[@]}";shift 2;;
        -n)
			local names=("${!2}")
			echo "names of added exp to: ${names[@]}";shift 2;;
        -t)
			local table=$2
			echo "table ouput by initial_comp set to: ${2}";shift 2;;
        -od)
			local out_dir=$2
			echo "name of out directory set to: ${2}";shift 2;;
		-h)
			usage add_coverage; return;;
		--help)
			usage add_coverage; return;;
		*)
			echo "Error in arguments"
			echo $1; usage add_coverage; return;;
	esac
done

local Errors=0
if [ -z $bedgraphs ]; then echo "ERROR: -b argument needed or needs to be a LIST";Errors+=1;fi
if [ -z $names ]; then echo "ERROR: -n argument needed or needs to be a LIST";Errors+=1;fi
if [ -z $table ]; then echo "ERROR: -t argument needed";Errors+=1;fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed";Errors+=1;fi

if [ ${#names[@]} -ne ${#bedgraphs[@]} ]; then
	echo "ERROR: -n & -b have to be lists of same length"; Errors+=1
fi

if [ $Errors -gt 0 ]; then usage add_coverage; return 1; fi

awk -v OFS="\t" 'NR==1 && $2!="begin" && $2!="start"{print $1,$2,$3,$4}NR>1{print $1,$2,$3,$4}' $table > $out_dir/tmp_peaks.bed
local peak_file=$out_dir/tmp_peaks.bed
local tmp_table=$out_dir/tmp_table.bed
cp $table $tmp_table

local list_files=()
local i=0

for bdg in ${bedgraphs[@]};
do
	echo $bdg
	echo ${names[$i]}
	echo 'here'
	bedtools intersect -a $peak_file -b $bdg -wa -wb -sorted -loj | awk  -v OFS="\t" '$6 != "-1" {print $0} $6=="-1" {print $1,$2,$3,$4,$1,$2,$3,0}' > $out_dir/${names[$i]}.inter
	echo 'here2'
	python $compute_coverage -i $out_dir/${names[$i]}.inter -m
	echo 'here3'
	awk  '{print (1000*$5)/($3-$2)}' $out_dir/${names[$i]}.inter.cov | sed "1i${names[$i]}" > $out_dir/${names[$i]}.inter.cov.normed
	list_files+=("$out_dir/${names[$i]}.inter.cov.normed")
	local i=$i+1
	echo 'here4'
done
# local files=`join_by " " "${list_files[@]}"`
#paste $tmp_table $files > $table
paste $tmp_table ${list_files[@]} > $out_dir/table.tsv
rm $peak_file $tmp_table

}


add_score(){
local pocc=false
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-m)
			local matrices=("${!2}")
			echo "matrice to computes scores column(s) set to: ${matrices[@]}";shift 2;;
        -n)
			local name=("${!2}")
			echo "name for score column(s) set to: ${name[@]}";shift 2;;
        -t)
			local table=$2
			echo "table set to: ${2}";shift 2;;
        -od)
			local out_dir=$2
			echo "output directory set to: ${2}";shift 2;;
        -po)
			local pocc=true
			echo "pocc mode activated";shift 2;;
		-h)
			usage add_coverage; return;;
		--help)
			usage add_coverage; return;;
		*)
			echo "Error in arguments"
			echo $1; usage add_coverage; return;;
	esac
done

local Errors=0
if [ -z $matrices ]; then echo "ERROR: -m argument needed or needs to be a LIST";Errors+=1;fi
if [ -z $name ]; then echo "ERROR: -n argument needed or needs to be a LIST";Errors+=1;fi
if [ -z $table ]; then echo "ERROR: -t argument needed";Errors+=1;fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed";Errors+=1;fi

if [ ${#name[@]} -ne ${#matrices[@]} ]; then
	echo "ERROR: -n & -m have to be lists of same length"; Errors+=1
fi

if [ $Errors -gt 0 ]; then usage add_coverage; return 1; fi
mkdir -p $out_dir
# awk -v OFS="\t" 'NR!=1{print $1, int($2+(($3-$2)/2)-25), int($2+(($3-$2)/2)+25), $4}' $table > $out_dir/tmp_peaks.bed
local peak_file=$out_dir/tmp_peaks.bed
local tmp_table=$out_dir/tmp_table.bed
local peak_fasta=$out_dir/tmp_fasta.fa
# sort -k1,1 -k2,2n $table > $tmp_table
cat $table > $tmp_table

if [[ $matrix == *".xml"* ]]; then
	awk -v OFS="\t" 'NR!=1{print $1, $2-2, $3, $4}' $table > $out_dir/tmp_peaks.bed
else
	awk -v OFS="\t" 'NR!=1{print $1, $2, $3, $4}' $table > $out_dir/tmp_peaks.bed
fi
bedtools getfasta -fi /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas -bed $peak_file -fo $peak_fasta

local list_files=()
local i=0
for matrix in ${matrices[@]};
do
	echo $matrix
	if [[ $matrix == *".xml"* ]]; then # TFFM scores computation
		$Python_TFFM $tffmscores -o ${out_dir}/tffm_scores.tsv -pos $peak_fasta -t ${matrix}
		sed 's/-/:/' ${out_dir}/tffm_scores.tsv | awk -v OFS="\t" -v FS="[:\t]" '$1!="None"{print $1,$2,$2+$4-2,$2+$5-1,$10;next}{print 0,0,0,0,0}' | awk -v OFS="\t" '{print $3,$4,$5}' | sed "1istartBS\tstopBS\t${name[$i]}" > ${out_dir}/${name[$i]}_scores.tsv 
		list_files+=("${out_dir}/${name[$i]}_scores.tsv")
	fi
	if [[ $matrix == *".txt"* ]]; then # K-mer scores computation
		echo "WIP"
	fi
	if [[ $matrix == *".pfm"* ]]; then  # PWM scores
		python $scores_prog -m ${matrix} -f $peak_fasta -o ${out_dir}/
# 		echo "${name[$i]}" > "${out_dir}/tab_${name[$i]}.tsv"
		sort -u -k1,1 -k8,8nr ${out_dir}/$(basename $peak_fasta).scores | sort -u -k1,1 | sed 's/-/:/' | awk -v OFS="\t" -v FS="[:\t]" '{print $1,$2,$2+$4-1,$2+$5-1,$10}' | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $3,$4,$5}' | sed "1istartBS\tstopBS\t${name[$i]}" > "${out_dir}/tab_${name[$i]}.tsv"
		list_files+=("${out_dir}/tab_${name[$i]}.tsv")
		
		# NOTE Use this to add type (exon, prom etc) to your table
# 		sort -u -k1,1 -k8,8nr ${out_dir}/$(basename $peak_fasta).scores | sort -u -k1,1 | sed 's/-/:/' | awk -v OFS="\t" -v FS="[:\t]" '{print $1,$2+$4,$2+$5,$1":"$2"-"$3,$10}' > ${out_dir}/tmp_${name[$i]}.tsv
# 		bedtools intersect -a ${out_dir}/tmp_${name[$i]}.tsv -b /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.bed -wao |sort -u -k4,4 -k10,10nr | sort -u -k4,4 | sed 's/-/:/' | awk -v FS="[:\t]" -v OFS="\t" '{print $4,$5,$6,$7,$11,$2,$3}' | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $5,$6,$7}' | sed "1itype\tstartBS\tstopBS" > ${out_dir}/tmp_${name[$i]}2.tsv
# 		list_files+=("${out_dir}/tmp_${name[$i]}2.tsv")
		
		if $pocc ; then # PWM Pocc
			python $pocc_pfm -s ${out_dir}/$(basename $peak_fasta).scores -o ${out_dir}/Pocc.pocc
			sort -nr ${out_dir}/Pocc.pocc > ${out_dir}/${name[$i]}_Pocc_scores.tsv
			list_files+=("${out_dir}${name[$i]}_Pocc_scores.tsv")
		fi
	fi
	local i=${i}+1
done
local files=`join_by " " "${list_files[@]}"`
paste $tmp_table $files > $out_dir/table.bed

# rm $peak_file $tmp_table

}


cooking_meth(){
echo "entered"
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-g)
			local genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
        -p)
			local peaksCov=$2
			echo "Tab-seaparated file for bound region with chr, start, end, covDAP, covAMP: ${2}";shift 2;;
		-m)
			local matrix_pfm=$2
			echo "TF pfm matrix: ${2}";shift 2;;
		#-l)
		#	local motifLength=$2
		#	echo "Motif length set to: ${2}";shift 2;;
		-c)
			local cutoff=$2
			echo "PFM score cutoff set to: ${2}";shift 2;;
		-c2)
			local cutoff2=$2
			echo "second step PFM score cutoff set to: ${2}";shift 2;;
		-sym)
			local symmetry=$2
			echo "[yes] or [no] to specifies wether or not the motif is symmetric: ${2}";shift 2;;
        -o)
			local outdir=$2
			echo "output directory set to: ${2}";shift 2;;
		-h)
			usage cooking_meth; return;;
		--help)
			usage cooking_meth; return;;
		*)
			echo "Error in arguments"
			echo $1; usage cooking_meth; return;;
	esac
done

local Errors=0
if [ -z $peaksCov ]; then echo "ERROR: -p argument needed";Errors+=1;fi
if [ -z $genome ]; then echo "ERROR: -g argument needed";Errors+=1;fi
if [ -z $matrix_pfm ]; then echo "ERROR: -m argument needed";Errors+=1;fi
#if [ -z $motifLength ]; then echo "ERROR: -l argument needed";Errors+=1;fi
if [ -z $cutoff ]; then echo "ERROR: -c argument needed";Errors+=1;fi
if [ -z $cutoff2 ]; then echo "ERROR: -c2 argument needed";Errors+=1;fi
if [ -z $outdir ]; then echo "ERROR: -o argument needed";Errors+=1;fi
if [ -z $symmetry ]; then echo "ERROR: -sym argument needed";Errors+=1;fi
if [ $Errors -gt 0 ]; then echo "error somewhere"; return 1; fi
tmp545=$(cat $matrix_pfm | wc -l)
motifLength=$(expr $tmp545 - 2)
echo $motifLength

# search the pfm against all bound regions
mkdir -p -m 774 $outdir/pfm_search


re='^[+-]?[0-9]+([.][0-9]+)?$'
if ! [[ $cutoff =~ $re ]] ; then
	if [ -f $cutoff ]; then
		local cutoff=$(awk '{print $3}' $cutoff2 )
		local cutoff2=$cutoff
	fi
fi

if [[ $peaksCov == *"_RiL_RiP.tsv"* ]]; then
	local header=$(awk 'NR==1{print tolower($7)}' $peaksCov)
	echo $header
	if [[ $header == *"ampdap"* ]]; then
		awk -v OFS="\t" 'NR!=1{print $1,$2,$3,$8,$7}' $peaksCov > $outdir/tmp_peaks.tsv
	else
		awk -v OFS="\t" 'NR!=1{print $1,$2,$3,$7,$8}' $peaksCov > $outdir/tmp_peaks.tsv
	fi
	local peaksCov=$outdir/tmp_peaks.tsv
fi
head $peaksCov
bedtools getfasta -fi $genome -bed $peaksCov -fo $outdir/pfm_search/all_peaks.fasta
python $scores_prog -m $matrix_pfm -f $outdir/pfm_search/all_peaks.fasta -o $outdir/pfm_search
pfmResults=$outdir/pfm_search/all_peaks.fasta.scores


python $full_methylation -p $peaksCov -g $genome -m $methMap -s $pfmResults -o $outdir -c $cutoff -l $motifLength

# the Rscript below will:
# compute methylation statistics
# plot these statistics
# and save these statistics as vectors into an RData object, so that it can be load into an R env to create pretty figures !
#exit 0
/home/prog/R/R-3-5-0/bin/Rscript $plot_meth_full $outdir $motifLength $symmetry $cutoff2
#$plot_meth_full $outdir $motifLength $symmetry $cutoff2
/home/prog/R/R-3-5-0/bin/Rscript $figs_meth_violin $outdir $matrix_pfm $symmetry
#$figs_meth_violin $outdir $matrix_pfm $symmetry
}


cons_scores(){
local explicit_mode=false
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
local negfile="none" #to avoid leaving an unset variable
	case $1 in
		-f)
			local filein=$2
			echo "input peak file: ${2}";shift 2;;
		-o)
			local outdir=$2
			echo "output directory set to: ${2}";shift 2;;
		-m)
			local matrix=$2
			echo "TF matrix: ${2}"; shift 2;;
		-e)
			local extend=$2
			echo "Extension required on both sides: ${2}"; shift 2;;
		-n)
			local negfile=$2
			echo "Bed file containing negative gene positions: ${2}"; shift 2;;
		-nset)
			local negset=$2
			echo "Number of negative set required: ${2}"; shift 2;;
		-seed)
			local seed=$2
			echo "Specified seed: ${2}"; shift 2;;
		--explicit)
			local explicit_mode=true
			echo "explicit mode required: extension will be integrated in output files"; shift 1;;
		*)
			echo "Error in arguments"
			echo $1; usage comparison; exit;;
	esac
done

if [ -z $seed ]; then local seed=34562; fi
if [[ $negfile == "none" ]]; then echo "no negative regions file given, shuffling will exclude input peak file sequences"; else echo -e "negative regions file: ${negfile}\nshuffling positive regions within negative regions given"; fi
if [ -z $negset ]; then echo "number of negative sets non specified, using default (1)"; local negset=1; fi


## usage
# cons_scores -f /home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/results/orthologs/conservation/tests/ChIP_DAP_DEG_peakscoord.bed -o /home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/results/orthologs/conservation/LFY/ChIP_DAP_DEG -m /home/312.6-Flo_Re/312.6.1-Commun/data/LFY.pfm -e 1000 -n /home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/results/neg_controls_DAP00001/all_nc_ATXG_positions.bed --explicit

## debug files for LFY
# outdir=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/results/orthologs/conservation/tests
# filein=$outdir/ChIP_DAP_DEG_peakscoord.bed
# matrix=/home/312.6-Flo_Re/312.6.1-Commun/data/LFY.pfm
# extend=1000
# negfile=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/results/neg_controls_DAP00001/all_nc_ATXG_positions.bed


## Necessary paths and files:
local A_thaliana_FASTA=/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas
# local scores_prog=/home/312.6-Flo_Re/312.6.1-Commun/scripts/TFgenomicsAnalysis/bin/scores.py
local get_bestscore_prog=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/scripts/select_bestscore.py
local genome_file=/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.txt
local plot_scores_prog=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/scripts/plot_cons_scores_v2.r
local bwtool_path=/home/312.3-StrucDev/312.3.1-Commun/bin/sl6
local R_36=/home/312.6-Flo_Re/312.6.1-Commun/Programs/Anaconda3/envs/py38_R36/bin/Rscript

mkdir -p $outdir ## if output directory doesn't exist, create it



## Extract filename from initial peaks file
local basename=$(basename $filein)
local file_extension="${basename##*.}"
local filename="${basename%.*}"
echo "filename: ${filename}"
echo "file extension: ${file_extension}"



## check if chrN coordinates in peak file are given correctly
if [[ $(awk 'NR==1{print $1}' $filein) != "chr"[0-9] ]]; then
	if [[ $(awk 'NR==1{print $1}' $filein) == [0-9] ]]; then
		awk -v OFS="\t" '{print "chr"$1,$2,$3}' $filein | sort -k1,3 | uniq > $outdir/${filename}.bed
		local filein=$outdir/${filename}.bed
		echo "chromosome coord: \"N\""
	elif [[ $(awk 'NR==1{print $1}' $filein) == "Chr"[0-9] ]]; then
		sed 's/Chr/chr/g' $filein | sort -k1,3 | uniq > $outdir/${filename}.bed
		local filein=$outdir/${filename}.bed
		echo "chromosome coord: \"ChrN\""
	else
		echo "unknown chromosome formatting mode, exiting..."
		exit 0
	fi
else
	echo "chromosome coordinates correctly formatted"
	if [[ ! -f $outdir/${filename}.bed ]]; then
		sort -k1,3 $filein | uniq > $outdir/${filename}.bed
	fi
	local filein=$outdir/${filename}.bed
fi



## If phylop and phastcons files are still bdg -> convert to bw
if [[ ! -f /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhyloP.bw ]]; then
	## Retrieve only chr1-5 from phastcons and phylop files
	grep -v "chrC" /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhastCons_chrN.bedGraph | grep -v "chrM" > /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhastCons_chrN.bedGraph
	
	grep -v "chrC" /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhyloP_chrN.bedGraph | grep -v "chrM" > /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhyloP_chrN.bedGraph
	
	## Convert bedgraph to bigwig for both files
	echo 'convert bedgraph files to bigwig'
	local bgtbw_path=/home/312.6-Flo_Re/312.6.1-Commun/Romain/ucscGenomeBrowser
	$bgtbw_path/bedGraphToBigWig /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhastCons_chrN.bedGraph /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.txt /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhastCons.bw
	
	$bgtbw_path/bedGraphToBigWig /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhyloP_chrN.bedGraph /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.txt /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhyloP.bw

fi

local phastcons=/home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhastCons.bw
local phylop=/home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhyloP.bw



## Get fasta sequences of input peak file, calculate scores and get best score position
if [[ ! -f $outdir/${filename}_TFBS_coord.bed ]]; then
	bedtools getfasta -fi $A_thaliana_FASTA -bed $filein -fo $outdir/${filename}.fa
	
	echo "PWM scores calculation"
	python $scores_prog -m $matrix -f $outdir/${filename}.fa -o $outdir
	
	
	## Only keep best site per sequence (greatest score)
	# 	python $get_bestscore_prog -f $outdir/${filename}.fa.scores -o $outdir
	sort -k1,1 -k8,8nr $outdir/${filename}.fa.scores \
		| awk -v OFS="\t" '{print $2,$4,$8,$1}' \
		| uniq -f3 \
		| awk -v OFS="\t" '{print $4,$1,$3}' > $outdir/${filename}_bestscore.bed
	
	
	## add best score position info to original peak bed file
	sed 's/:/\t/g' $outdir/${filename}_bestscore.bed | sed 's/-/\t/' > $outdir/${filename}_bestscore_position.bed
	# old method after first sed: awk -v OFS="\t" '{gsub("-","\t",$2)}1' 
	
	
	## Cut coordinates around best TFBS 
	echo "retrieving TFBS coordinates..."
	if [[ $matrix == "/home/312.6-Flo_Re/312.6.1-Commun/data/LFY.pfm" ]]; then # calculating matrix length
		local matrixlength=19
	else
		local matrixlength=$(tail -n+3 $matrix | wc -l)
	fi

	## Retrieve TFBS coordinates
	awk -v OFS="\t" -v fullmatrixlength=$matrixlength '{print $1,$2+$4-1,$2+$4+fullmatrixlength-1}' $outdir/${filename}_bestscore_position.bed | sort -k1,1 -k2,3n -u > $outdir/${filename}_TFBS_coord.bed
	
	
fi



################ Extend and retrieve avg score per position ################ 
####### METHOD FOR SHUFFLE AS IN PAPER #######
echo "now using real shuffling method"

if [[ ! -f $outdir/${filename}_shuffled_PhyloP.txt ]]; then
	echo "bwtool is going to be used now for positive controls"
# 	echo "bwtool for TFBS coordinates"
	## Use these coordinates as input for bwtool
	$bwtool_path/bwtool aggregate "$extend:$extend" $outdir/${filename}_TFBS_coord.bed $phastcons $outdir/${filename}_shuffled_PhastCons.txt

	$bwtool_path/bwtool aggregate "$extend:$extend" $outdir/${filename}_TFBS_coord.bed $phylop $outdir/${filename}_shuffled_PhyloP.txt
fi


if [[ ! -f $outdir/${filename}_shuffled_PhyloP_medstdev.txt ]]; then
	echo "bwtool is going to be used to retrieve median stdev and number of values used, for all positive controls"
# 	echo "bwtool for TFBS coordinates"
	## Use these coordinates as input for bwtool
	$bwtool_path/bwtool aggregate "$extend:$extend" $outdir/${filename}_TFBS_coord.bed $phastcons $outdir/${filename}_shuffled_PhastCons_medstdev.txt -expanded

	$bwtool_path/bwtool aggregate "$extend:$extend" $outdir/${filename}_TFBS_coord.bed $phylop $outdir/${filename}_shuffled_PhyloP_medstdev.txt -expanded
fi


## if there are over 50k peaks in bed file, use only one negative set
if [[ $(wc -l <$outdir/${filename}_TFBS_coord.bed) -gt 50000 ]]; then
	echo $(wc -l <$outdir/${filename}_TFBS_coord.bed)
	negset=1
	echo $negset
fi

# exit 0

for ((i=1;i<=$negset;i++)); do
	if [[ ! -f $outdir/${filename}_shuffled_nset_all_PhyloP.txt ]]; then 
		## shuffle of these coordinates
		echo "shuffling..."
		((seed=seed+i))
		
		if [[ $negfile == "none" ]]; then
			bedtools shuffle -i $outdir/${filename}_TFBS_coord.bed -g $genome_file -chrom -seed ${seed} | sort -k1,3 | uniq > $outdir/${filename}_TFBS_coord_shuffled_nset${i}.bed
		else 
			bedtools shuffle -i $outdir/${filename}_TFBS_coord.bed -g $genome_file -incl $negfile -chrom -seed ${seed} | sort -k1,3 | uniq > $outdir/${filename}_TFBS_coord_shuffled_nset${i}.bed
		fi
		## bwtool of the coordinates
		echo "looking for phastcons and phylop coordinates for set ${i}"
		$bwtool_path/bwtool aggregate "$extend:$extend" $outdir/${filename}_TFBS_coord_shuffled_nset${i}.bed $phastcons $outdir/${filename}_shuffled_nset${i}_PhastCons.txt

		$bwtool_path/bwtool aggregate "$extend:$extend" $outdir/${filename}_TFBS_coord_shuffled_nset${i}.bed $phylop $outdir/${filename}_shuffled_nset${i}_PhyloP.txt
	fi
done


if [[ ! -f $outdir/${filename}_shuffled_nset_all_PhyloP.txt ]]; then
	if [[ $negset != "1" ]]; then
		multijoin $outdir/${filename}_shuffled_nset_all_PhastCons.txt $outdir/${filename}_shuffled_nset[0-9]_PhastCons.txt
		ls $outdir | grep -P "^${filename}_shuffled_nset[0-9]+_PhastCons.txt" | awk -v dir=$outdir '{print dir"/"$1}' | xargs -d"\n" rm
		# rm $outdir/${filename}_shuffled_nset[0-9]_PhastCons.txt 
		
		multijoin $outdir/${filename}_shuffled_nset_all_PhyloP.txt $outdir/${filename}_shuffled_nset[0-9]_PhyloP.txt
		ls $outdir | grep -P "^${filename}_shuffled_nset[0-9]+_PhyloP.txt" | awk -v dir=$outdir '{print dir"/"$1}' | xargs -d"\n" rm
		# rm $outdir/${filename}_shuffled_nset[0-9]_PhyloP.txt 
		
		## remove negsets
		ls $outdir | grep -P "^${filename}_TFBS_coord_shuffled_nset[0-9]+.bed$" | awk -v dir=$outdir '{print dir"/"$1}' | xargs -d"\n" rm
	else
		mv $outdir/${filename}_shuffled_nset1_PhastCons.txt $outdir/${filename}_shuffled_nset_all_PhastCons.txt
		mv $outdir/${filename}_shuffled_nset1_PhyloP.txt $outdir/${filename}_shuffled_nset_all_PhyloP.txt
	fi
fi

echo "now plotting..."
local R_cons_plots=/home/312.6-Flo_Re/312.6.1-Commun/Programs/Anaconda3/envs/LFYUFO_figs/bin/Rscript
## plot
if [[ ! -f $outdir/${filename}_shuffled_phastcons_phylop.pdf ]] || [[ ! -f $outdir/${filename}_shuffled_phastcons_phylop_${extend}.pdf ]]; then
	if [[ $explicit_mode == "true" ]]; then
		
		$R_cons_plots $plot_scores_prog -n ${filename}_shuffled -i $outdir -o $outdir -e $extend --explicit
# 		conda run -n LFYUFO_figs Rscript $plot_scores_prog -n ${filename}_shuffled -i $outdir -o $outdir -e $extend --explicit
	else
		
		$R_cons_plots $plot_scores_prog -n ${filename}_shuffled -i $outdir -o $outdir -e $extend
# 		conda run -n LFYUFO_figs Rscript $plot_scores_prog -n ${filename}_shuffled -i $outdir -o $outdir -e $extend
	fi
fi

}
