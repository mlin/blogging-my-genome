#!/bin/bash

main() {

    set -x

    dx download "$your_vcf" -o /your.vcf.gz
    vcfname=$(dx describe --name "$your_vcf")
    vcfname=${vcfname%.vcf.gz}

    dx download "$pop_vcf" -o /pop.vcf.gz

    if [ -z "$Rlibs" ]; then
        R --no-save --quiet < /libs.R
        tar -C / -cf /1000g_unphased_pca.Rlibs.tar /usr/local/lib/R/site-library
        gzip /1000g_unphased_pca.Rlibs.tar
        output_Rlibs=$(dx upload /1000g_unphased_pca.Rlibs.tar.gz --brief)
        dx-jobutil-add-output Rlibs "$output_Rlibs" --class=file
    else
        dx download "$Rlibs" -o /1000g_unphased_pca.Rlibs.tar.gz
        tar -C / -zxf /1000g_unphased_pca.Rlibs.tar.gz
    fi

    mkdir /figures
    sed -i "s/YOUR_NAME/$your_name/" /generate.R
    R --no-save --quiet < /generate.R

    while IFS= read -r -d $'\0' file; do
        dx-jobutil-add-output figures "$(dx upload $file --brief)" --class=file --array
    done < <(find /figures -type f -print0)
}
