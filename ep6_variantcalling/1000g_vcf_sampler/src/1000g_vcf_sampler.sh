#!/bin/bash
# 1000g_vcf_sampler 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://wiki.dnanexus.com/Developer-Portal for tutorials on how
# to modify this file.

sample_chrom() {
    url=$(printf 'https://s3.amazonaws.com/1000genomes/release/20110521/ALL.chr%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz' $2)
    curl -# "$url" | zcat \
        | awk "substr(\$0,1,1) == \"#\" { print \$0; next } 100*rand() < $1 { print \$0 }" \
        > "$2.vcf"
}

main() {
    set -ex

    echo "autosomes_only: $autosomes_only"

    chromosomes=(`seq 1 22`)
    maybe_autosomes=".autosomes"
    if [[ "$autosomes_only" != "true" ]]; then
        chromosomes=("${chromosomes[@]}" X)
        maybe_autosomes=""
    fi
    echo "chromosomes: ${chromosomes[@]}"

    export -f sample_chrom
    printf -- '%s\n' "${chromosomes[@]}" | xargs -P 8 -n 1 -i bash -e -c 'sample_chrom "$@"' _ "$percentage" {}

    ls -tlh

    vcfgz=$(printf -- '%s.vcf\n' "${chromosomes[@]}" | xargs -n 999999 vcfcombine | bgzip -c \
            | dx upload -o "1000G_genotypes${maybe_autosomes}_${percentage}percent.vcf.gz" --brief -)

    dx-jobutil-add-output vcfgz "$vcfgz" --class=file
}