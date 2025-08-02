#!/bin/bash
set -e

primers=data/Primers_120plex.bed
genome=data/ucsc.hg38.fa
adapter_seq="CCTACACGACGCTCTTCCGATCT"
mkdir -p results/summary
summary_file="results/summary/summary.tsv"

echo -e "sample\ton_target_pct\toff_target_pct\ttotal_reads\toff_target_regions\toff_with_primers\toff_in_repeats\tmean_gc_off\tmedian_dist_to_primer\tmin_dist_to_primer\tnonspecific_primers_overlap" > $summary_file

mkdir -p results/gc_content

if [ ! -f "$genome.bwt" ]; then
	echo "Indexing genome for BWA"
	bwa index $genome
fi

echo "Extracting primer sequences"
bedtools getfasta -fi $genome -bed $primers -fo results/primers.fa

echo "Checking the alignment of primers across the genome"
bwa mem $genome results/primers.fa > results/primers.sam

echo "Converting SAM to BAM and sorting"
samtools view -Sb results/primers.sam | samtools sort -o results/primers.sorted.bam
samtools index results/primers.sorted.bam

echo "Identifying non-specific primers"
samtools view results/primers.sorted.bam | cut -f1 | sort | uniq -c > results/primers_alignment_counts.txt
while read count coord; do
  [ "$count" -gt 1 ] && echo "$coord"
done < results/primers_alignment_counts.txt > results/nonspecific_primers.txt

echo "Analyzing GC content"
bedtools nuc -fi $genome -bed $primers > results/gc_content/primers_gc.txt

for sample_dir in data/VC_RnD_SPET_13_*; do
	sample_name=$(basename "$sample_dir")
	echo "Analyzing $sample_name"
	mkdir -p "results/$sample_name"
	python analyze_offtarget.py "$sample_dir"
	sort -k1,1 -k2,2n "results/$sample_name/top_regions.bed" > "results/$sample_name/top_regions.sorted.bed"
	bedtools intersect -a data/repeats.bed -b "results/$sample_name/top_regions.sorted.bed" -wo > "results/$sample_name/offtarget_repeats_overlap.bed"
	sort -k1,1 -k2,2n $primers > "results/primers.sorted.bed"
	bedtools intersect -a "results/$sample_name/top_regions.sorted.bed" \
		-b "results/primers.sorted.bed" -wa -wb > "results/$sample_name/offtarget_with_primers.bed"
	bedtools nuc -fi $genome -bed "results/$sample_name/top_regions.bed" > "results/$sample_name/top_regions.gc.txt"
	bedtools closest -a "results/$sample_name/top_regions.sorted.bed" -b "results/primers.sorted.bed" -d > "results/$sample_name/closest_primers.bed"
	cp "results/$sample_name/closest_primers.bed" "results/$sample_name/closest_primers.tsv"
	awk -F'\t' '($NF <= 500) && ($NF != -1)' "results/$sample_name/closest_primers.tsv" > "results/$sample_name/closest_primers_500bp.tsv"
	bedtools closest -a "results/$sample_name/top_regions.sorted.bed" -b "results/primers.sorted.bed" -d | awk '$NF <= 1000' > "results/$sample_name/filtered_closest_primers.bed"
	

	qc_json="$sample_dir/QC_RnD_SPET_13_*.json"
	on_target_pct=$(jq -r '.["On target"]["aftere UMI processing"]["M3"]' $qc_json)
	off_target_pct=$(awk -v x="$on_target_pct" 'BEGIN {print 100 - x}'| sed 's/,/./')
	distances=$(cut -f8 "results/$sample_name/closest_primers.bed" | grep -vE '^-1$')
    	if [ -n "$distances" ]; then
		sorted_dist=$(echo "$distances" | sort -n)
    		count=$(echo "$sorted_dist" | wc -l)
		min_dist=$(echo "$sorted_dist" | head -n1)
		half=$((count/2))
    		if [ $((count % 2)) -eq 1 ]; then
        		median_dist=$(echo "$sorted_dist" | tail -n +$((half + 1)) | head -n1)
    		else
        		val1=$(echo "$sorted_dist" | tail -n +$half | head -n1)
        		val2=$(echo "$sorted_dist" | tail -n +$((half + 1)) | head -n1)
        		median_dist=$(( (val1 + val2) / 2 ))
    		fi
	else
    		median_dist=0
    		min_dist=0
	fi
	mean_gc=$(tail -n +2 "results/$sample_name/top_regions.gc.txt" | cut -f5 | paste -sd+ - | bc -l)
	count=$(tail -n +2 "results/$sample_name/top_regions.gc.txt" | wc -l)
	mean_gc=$(echo "scale=4; $mean_gc / $count" | bc -l)
	mean_gc=$(echo "$mean_gc" | sed 's/^\./0./')
	repeat_hits=$(cut -f4-6 "results/$sample_name/offtarget_repeats_overlap.bed" | sort -u | wc -l)
	nonspecific_hits=$(grep -f results/nonspecific_primers.txt "results/$sample_name/offtarget_with_primers.bed" | wc -l)
	intersect_hits=$(wc -l < "results/$sample_name/offtarget_with_primers.bed")
	total_reads=$(jq '.["Total reads PROCESSED"]["total"]' $qc_json)
	off_bed=$(find "$sample_dir" -name "RnD_SPET_*_off_target_M3.bed" | head -1)
	num_regions=$(cut -f1-3 "$off_bed" | sort -u | wc -l)
	echo -e "$sample_name\t$on_target_pct\t$off_target_pct\t$total_reads\t$num_regions\t$intersect_hits\t$repeat_hits\t$mean_gc\t$median_dist\t$min_dist\t$nonspecific_hits" >> $summary_file

	echo "Finished $sample_name"

done
(head -n1 results/summary/summary.tsv && tail -n +2 results/summary/summary.tsv | sort -V) > results/summary/summary_sorted.tsv
echo "Finished all samples"
