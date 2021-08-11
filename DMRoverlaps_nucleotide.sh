# GAT Enrichments
# Ben Laufer

cd /share/lasallelab/Ben/PEBBLES/DNA/DMRs/GAT/

echo "Sorting bed files"

parallel --dry-run --will-cite "sort -k1,1 -k2,2n {} > {.}_sorted.bed" ::: *.bed
parallel --verbose --will-cite "sort -k1,1 -k2,2n {} > {.}_sorted.bed" ::: *.bed

echo "Done sorting bed files"

module load gat/1.3.4 

echo "Running enrichment analysis"

call="gat-run.py \
--segments=placenta_consensus_sigRegions_mm10_sorted.bed \
--annotations=brain_consensus_sigRegions_mm10_sorted.bed \
--workspace=placenta_consensus_regions_mm10_sorted.bed \
--counter=nucleotide-overlap \
--num-samples=10000 \
--num-threads=20 \
--log=BrainAndPlacenta.log \
> BrainAndPlacenta_GAT_results.tsv"

echo $call
eval $call

echo "Done running enrichment analysis"

echo "Removing temporary files"
rm *_sorted.bed
echo "Done removing temporary files"
