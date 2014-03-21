#!/bin/bash

main() {
	# unpack vcftools
	echo "Unpacking vcftools"
	mkdir -p /usr/local/vcftools
	sh -c "cd /usr/local/vcftools && tar xf /tmp/vcftools_built.tar"
	export PERL5LIB=/usr/local/vcftools/perl

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
	done

	# for parallelization: make list of reference sequences in order of
	# decreasing size
	# TODO: even better would be order of decreasing mapped reads, which
	#       we could get from 'samtools idxstat'
	sort -k 2rn,2rn genome.fa.fai | cut -f1 > /tmp/rseqs

	sleep 1
	echo "Reference sequences:"
	cat /tmp/rseqs
	sleep 1

	#
	# Run freebayes
	#
	echo "Executing parallel freebayes"
	xargs -a /tmp/rseqs -n 1 -P 8 -r -I XXX --verbose \
		freebayes $input -v XXX.vcf -f genome.fa -r XXX $advanced_options

	ls -1sh *.vcf

	# for vcf-concat: make ordered list of .vcf file names
	cut -f1 genome.fa.fai | awk '{printf("%s.vcf\n",$0)}' > /tmp/rseqs.vcfs

	#
	# Upload results
	#
	name=`dx describe "${sorted_bams[0]}" --name`
	name="${name%.bam}"

	echo "Uploading results"
	file_id=`xargs -a /tmp/rseqs.vcfs -n 999999 /usr/local/vcftools/bin/vcf-concat | bgzip -c | dx upload -o "$name.vcf.gz" --brief -`
	dx-jobutil-add-output "variants_vcfgz" "$file_id"
}
