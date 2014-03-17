#!/bin/bash

main() {
    echo "Value of bam: '$bam'"
    echo "Value of features: '$features'"

    dx download -o features "$features"

    bam_name=$(dx describe --name "$bam")
    dx cat "$bam" | bedtools coverage -abam stdin -b features -d > raw
    cat raw | cut -f1-3 | sort -S 4G -u | wc -l
    report=$(cut -f5 raw | gzip -c | dx upload -o "${bam_name%.bam}.feature_basewise_coverage.gz" --brief -)

    dx-jobutil-add-output report "$report" --class=file
}
