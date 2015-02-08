#!/bin/bash

main() {
    set -e -x -o pipefail

    if [ -z "$output_name "]; then
        output_name=$(dx describe "${reads[0]}" --name)
        output_name=${output_name%.fastq.gz}
    fi

    # parallelized trimadap of all the reads
    pids=()
    trimmed_reads=()
    for fqgz in "${reads[@]}"; do
        nm=$(dx describe "$fqgz" --name)
        nm=${nm%.fastq.gz}
        dx cat "$fqgz" | trimadap | gzip -1 > "${nm}.trimmed.fastq.gz" &
        pids+=($!)
        trimmed_reads+=("${nm}.trimmed.fastq.gz")
    done
    for pid in "${pids[@]}"; do wait $pid || exit $?; done

    ls -lh *.trimmed.fastq.gz

    # run fermi2
    fermi2.pl unitig -t $(nproc) -p unitig -s "$genome_size" -l "$read_length" "zcat $trimmed_reads"  | tee unitig.mk
    make -f unitig.mk
    ls -lh

    # upload outputs
    zip -c logs.zip unitig*.log
    dx-jobutil-add-output logs --class=file \
        $(dx upload --brief --destination "${output_name}.fermi2_unitig.logs.zip" logs.zip)
    dx-jobutil-add-output fmd_index --class=file \
        $(dx upload --brief --destination "${output_name}.fmd" unitig.flt.fmd)
    dx-jobutil-add-output unitigs --class=file \
        $(dx upload --brief --destination "${output_name}.mag.gz" unitig.mag.gz)
    
}
