# homer
# Parallel HOMER known transcription factor motif analyses for multiple DMRichR comparisons
# Ben Laufer

module load homer

homerDMR(){
	i=${1}
	
	echo
	cd ${i}/Extra/HOMER
	echo "Running HOMER in ${PWD}"
	
    echo 
    echo "Testing DMRs from ${i}"
	mkdir both
	cp ../../DMRs/DMRs.bed .
	
	call="findMotifsGenome.pl \
	DMRs.bed \
	mm10 \
	both/ \
	-bg background.bed \
	-cpg \
	-size given \
	-p 15 \
	-nomotif"

	echo $call
	eval $call
}
export -f homerDMR

cd /share/lasallelab/Ben/PEBBLES/DNA/DMRs/
mkdir homerLogs

parallel --dry-run --will-cite --results homerLogs -j 4 "homerDMR {}" ::: brain_male brain_female placenta_male placenta_female
parallel --verbose --will-cite --results homerLogs -j 4 "homerDMR {}" ::: brain_male brain_female placenta_male placenta_female
