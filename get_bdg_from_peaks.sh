## script to get RIP-normalized bdg files for IGB peak visualization
source /home/312.6-Flo_Re/312.6.1-Commun/scripts/TFgenomicsAnalysis/compil_functions.sh

main_dir=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY-UFO/results
mapping_dir=$main_dir/Mapping
dir_peakcalling=$main_dir/PeakCalling
dir_comparisons=$main_dir/comparisons




### Produce RIP bedgraph files
get_RIP_RIL_bdgs (){
# get_RIP_bdgs -n <STRING> -p <FILE> -bam <PATH> -snames <LIST> -od <PATH> -pd <PATH> -m <STRING> -s <INT>
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-n)
			local name_cons=$2
			echo "consensus name set to: ${2}";shift 2;;
		-p)
			local peaks=$2
			echo "merged peaks file set to:${2}";shift 2;;
		-bam)
			local bam_dir=$2
			echo "bam directory set to: ${2}";shift 2;;
		-sn)
			local samples_names=("${!2}")
			echo "list of sample's name set to: ${samples_names[@]}";shift 2;;
		-od)
			local out_dir=$2
			 echo "output directory set to: ${2}";shift 2;;
		-pd)
			local peaks_dir=$2
			 echo "directory where original peakcalling data is stored: ${2}";shift 2;;
		-m)
			local mode=$2
			 echo "normalization mode requested: ${2}";shift 2;;
		-s)
			local seedrandom=$2
			echo "random seed set to: ${2}"; shift 2;;
		-h)
			usage peakcalling_MACS2; return;;
		--help)
			usage peakcalling_MACS2; return;;
		*)
			echo "Error in arguments"
			echo $1; usage peakcalling_MACS2; return;;
	esac
done



## first compute rpkmRIL or RIP for the set of original rep
mkdir -p $out_dir
# echo $peaks_dir
# echo $bam_dir
# echo $peaks
compute_rpkmrip_rpkmril -p $peaks -bd bam_dir -pd peaks_dir -sn samples_names[@] -m $mode -o $out_dir/${mode}/$name_cons

# exit 0

# exit 0
local out_dir=$out_dir/${mode}/$name_cons

# exit 0



i=1

## for every rep sample, create a new bdg from
## the initial bam file with the desired scaling
for name in ${samples_names[@]}; do
	echo $name
	
	local log=$peaks_dir/${name}/${name}.log
	
	local fragment_length=$(cat $log | grep "predicted fragment length is" | awk -v OFS="\t" '{print $13}')
	local in_dir=$bam_dir/${name}
# 	

	local bam=$(find $in_dir -name "*.filtered.sorted.nodup.bam" -type f | grep -v "control")
	
	local RIP=$(awk -v i=$i 'NR==i{print $0}' $out_dir/RIP.txt)
	echo "RIP: ${RIP}"

	local libsize=$(awk -v i=$i 'NR==i{print $0}' $out_dir/tmpTotalTags.txt) # used to be tmpFiltTags
	echo "libsize: $libsize" 
	
	
	if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
		R2_value=$(awk -v FS=" " 'NR>1{print $2}' ${bam%.filtered.sorted.nodup.bam}.minimal.stats)
		echo $R2_value
	fi
	
	if [[ $mode == "inPeaks" ]]; then
		local scale=$(calc 1000000/$RIP)
		m="RIP"
	elif [[ $mode == "inLibs" ]]; then
		local scale=$(calc 1000000/$libsize)
		m="RIL"
	fi
	
	
	if [[ ${R2_value} == "NA" ]]; then ## if single end
		echo "single end"
		if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
			echo $bam
			
			if [[ ! -f $out_dir/${in_dir##*/}.bamtobed.bed ]]; then
				bedtools bamtobed -i $bam > $out_dir/${in_dir##*/}.bamtobed.bed
			fi
			
			if [[ ! -f $out_dir/${in_dir##*/}.ext.bed ]]; then
				awk -v fraglen=${fragment_length} -v OFS="\t" '$6=="+"{print $1,$2,$3+fraglen,$4,$5,$6}$6=="-"{print $1,$2-fraglen,$3,$4,$5,$6}' $out_dir/${in_dir##*/}.bamtobed.bed | awk -v OFS="\t" '$2<=0{print $1,1,$3,$4,$5,$6;next}{print $0}' > $out_dir/${in_dir##*/}.ext.bed
			fi
			
			if [[ ! -f $out_dir/${in_dir##*/}_${m}_cov.bdg ]]; then
				sed -i 's/chr/Chr/g' $out_dir/${in_dir##*/}.ext.bed  ## Laura 27/05/2021
				bedtools genomecov -bga -scale $scale -i $out_dir/${in_dir##*/}.ext.bed -g /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.size | sed 's/Chr/chr/g' > $out_dir/${in_dir##*/}_${m}_cpm.bdg ## Laura 27/05/2021
				
				cp $out_dir/${in_dir##*/}_${m}_cpm.bdg $out_dir/${in_dir##*/}_${m}_cov.bdg
				rm $out_dir/${in_dir##*/}_${m}_cpm.bdg
			fi
			
		fi
	
	else ## if paired end
		echo "paired end"
		if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
			bedtools genomecov -bga -scale $scale -ibam $bam > $out_dir/${in_dir##*/}_${m}_cpm.bdg
			cp $out_dir/${in_dir##*/}_${m}_cpm.bdg $out_dir/${in_dir##*/}_${m}_cov.bdg
			rm $out_dir/${in_dir##*/}_${m}_cpm.bdg
		fi
# 			bedtools genomecov -bga -scale $scale -i $out_dir/${in_dir##*/}.ext.bed -g /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.size | sed 's/Chr/chr/g' > $out_dir/${in_dir##*/}_${m}_cpm.bdg ## Laura 27/05/2021
				
# 			cp $out_dir/${in_dir##*/}_${m}_cpm.bdg $out_dir/${in_dir##*/}_${m}_cov.bdg
			
		
		
	fi
	((i++))
done


## retrieve all bdgs and get average coverage
list_bdg=$(find $out_dir -name "*_cov.bdg" -type f | grep -v "control")
echo ${list_bdg[@]}

bedtools unionbedg -i ${list_bdg[@]} > $out_dir/${name_cons}_${m}_cov_sep.bdg

awk -v OFS="\t" '{moy=0;nb=0;for(i=4;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/${name_cons}_${m}_cov_sep.bdg > $out_dir/${name_cons}_${m}_cov.bdg

rm $out_dir/${name_cons}_${m}_cov_sep.bdg

}




