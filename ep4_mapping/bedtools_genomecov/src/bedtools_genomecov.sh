#!/bin/bash

main() {
    echo "Value of bam: '$bam'"

    bam_name=$(dx describe --name "$bam")
    report=$(dx cat "$bam" | bedtools genomecov -ibam stdin | dx upload -o "${bam_name%.bam}.genomecov" --brief -)

    dx-jobutil-add-output report "$report" --class=file
}
