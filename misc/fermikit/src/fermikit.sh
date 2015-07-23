#!/bin/bash

main() {
    set -ex -o pipefail

    echo "" | gzip -c > reads.fq.gz
    for i in ${!reads[@]}
    do
        dx cat "${reads[$i]}" >> reads.fq.gz
    done

    ls -lh

    echo "SHELL := bash
.SHELLFLAGS := -e -o pipefail -c" > prefix.mak
    fermi.kit/fermi2.pl unitig -s${genome_size} -t$(nproc) -l${read_length} -p prefix reads.fq.gz >> prefix.mak
    cat prefix.mak
    make -f prefix.mak
    echo $?

    ls -lh

    if [ -z "$output_name" ]; then
        output_name=$(dx describe --name "${reads[0]}")
        output_name=${output_name%.gz}
        output_name=${output_name%.fastq}
        output_name=${output_name%.fq}
    fi

    unitigs=$(zcat prefix.mag.gz | pigz -c | dx upload --destination "${output_name}.mag.gz" --brief -)
    dx-jobutil-add-output unitigs "$unitigs" --class=file

    fmd_index=$(dx upload --destination "${output_name}.fmd" --brief prefix.flt.fmd)
    dx-jobutil-add-output fmd_index "$fmd_index" --class=file

    logs=$(tar -cvz *.log | dx upload --destination "${output_name}.fermikit_logs.tar.gz" --brief -)
    dx-jobutil-add-output logs "$logs" --class=file
}
