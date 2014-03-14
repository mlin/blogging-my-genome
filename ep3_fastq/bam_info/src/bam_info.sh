#!/bin/bash
# bam_info 0.0.1
# Generated by dx-app-wizard.

bam_info() {
    set -ex

    echo "-- FLAGSTAT"
    dx cat "$1" | tee input.bam | pigz -dc | /samtools flagstat -

    echo "

-- COUNTS BY READ GROUP"

    # list read groups
    read_groups=($(/samtools view -H input.bam | awk '$1 == "@RG" {split($2,a,":"); print a[2]; }'))
    for rg in "${read_groups[@]}"
    do
        echo -n "${rg}  "
        cat input.bam | pigz -dc | /samtools view -r "$rg" -c -
    done
    
    rm -f input.bam.bai
    samtools index input.bam || true
    # samtools index may fail if bam isn't coordinate-sorted
    if [ -f input.bam.bai ] ; then
        echo "

-- IDXSTAT"
        samtools idxstats input.bam
    fi

    echo "

-- FULL BAM HEADER"
    /samtools view -H input.bam
}

main() {
    echo "Value of bam: '$bam'"

    # boilerplate stuff needed to get dnanexus samtools running in the ubuntu
    # 12.04 execution environment
    rm /etc/apt/apt.conf.d/99dnanexus
    add-apt-repository -y ppa:ubuntu-toolchain-r/test
    apt-get update
    apt-get install -y gcc-4.8 g++-4.8
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 50
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.8 50
    bzip2 -d /samtools.bz2
    chmod +x /samtools

    for one_bam in "${bam[@]}"
    do
        bam_fid=$(dx-jobutil-parse-link "$one_bam")
        bam_name=$(dx describe --name "$bam_fid")
        echo "$bam_name"
        bam_info "$bam_fid" > report.txt
        bam_name=$(basename "$bam_name" .bam)
        report_fid=$(dx upload -o "${bam_name}.bam_info.txt" --brief report.txt)
        dx-jobutil-add-output --class array:file report "$report_fid"
    done
}
