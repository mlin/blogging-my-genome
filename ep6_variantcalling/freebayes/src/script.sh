#!/bin/bash

main() {
	#
	# Fetch and index genome
	#
	echo "Fetching & indexing reference genome"
	dx download "$genome_fastagz" -o genome.fa.gz
	gunzip genome.fa.gz
	samtools faidx genome.fa

	#
	# Fetch and index mappings
	#
	echo "Fetching & indexing BAM"
	input=""
	for i in "${!sorted_bams[@]}"; do
	  dx download "${sorted_bams[$i]}" -o "input-$i.bam"
	  samtools index "input-$i.bam"
	  input="$input -b input-$i.bam"
	  samtools view -H "input-$i.bam" | awk -v OFS="\t" '$1 == "@SQ" {split($2,a,":"); split($3,b,":"); print a[2], b[2]; }' >> /tmp/tseqs.raw
	done

	# for parallelization: make list of target sequences in order of
	# decreasing size
	sort -nr -k 2,2 /tmp/tseqs.raw | cut -f1 | huniq > /tmp/tseqs

	sleep 1
	echo "Target sequences:"
	cat /tmp/tseqs
	sleep 1

	#
	# Run freebayes
	#
	echo "Executing parallel freebayes"
	xargs -a /tmp/tseqs -n 1 -P 8 -r -I XXX --verbose \
		freebayes $input -v XXX.vcf -f genome.fa -r XXX $advanced_options

	#
	# Upload results
	#
	name=`dx describe "${sorted_bams[0]}" --name`
	name="${name%.bam}"

	echo "Uploading results"
	file_id=`cat /tmp/tseqs.raw | cut -f1 | huniq | xargs -I XXX cat XXX.vcf | bgzip -c | dx upload -o "$name.vcf.gz" --brief -`
	dx-jobutil-add-output "variants_vcfgz" "$file_id"
}
