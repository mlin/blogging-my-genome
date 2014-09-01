#!/bin/bash
# gemini_load 0.0.1

set -e -x -o pipefail

run_vep() {
    # install apt dependencies
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-suggests --no-install-recommends libdbi-perl bioperl

    # deploy VEP cache
    mkdir -p /vep/cache
    pushd /vep/cache
    dx cat file-BGYpvVj0p0PYy1FG7GyQ0j6Q | tar x # homo_sapiens_vep_75.tar

    # deploy VEP code
    cd /vep
    unzip -q ensembl-tools-75.zip
    dx download file-BPGyxyj0FVg3GyPPVx9QpfyP -o ensembl-api.tar.gz
    tar zxf ensembl-api.tar.gz
    export PERL5LIB=${PERL5LIB}:/vep/ensembl/modules:/vep/ensembl-compara/modules:/vep/ensembl-variation/modules:/vep/ensembl-funcgen/modules

    # download input VCF
    dx download "$1" -o input.vcf.gz

    # run VEP
    cd ensembl-tools-release-*/scripts/variant_effect_predictor
    perl variant_effect_predictor.pl --cache --offline --dir /vep/cache/.vep --species homo_sapiens \
                                     --sift b \
                                     --polyphen b \
                                     --symbol \
                                     --numbers \
                                     --biotype \
                                     --total_length \
                                     --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \
                                     -i /vep/input.vcf.gz --vcf --fork $(nproc)
    # VEP command-line ref: http://gemini.readthedocs.org/en/latest/content/functional_annotation.html#stepwise-installation-and-usage-of-vep

    # stage annotated VCF
    cp variant_effect_output.txt /home/dnanexus/my.vcf
    popd
}

setup_gemini() {
    wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
    mkdir -p /usr/local /usr/local/share/gemini
    python gemini_install.py /usr/local /usr/local/share/gemini
    export PATH=$PATH:/usr/local/gemini/bin

    gemini update --dataonly --extra cadd_score
    gemini update --dataonly --extra gerp_bp
}

main() {
    echo "Value of vcfgz: '$vcfgz'"
    echo "Value of vep: '$vep'"
    echo "Value of advanced_options: '$advanced_options'"

    setup_gemini &

    if [ "$vep" = "true" ]; then
        run_vep "$vcfgz"
        advanced_options="-t VEP $advanced_options"
    else
        dx cat "$vcfgz" | zcat > /home/dnanexus/my.vcf
    fi

    wait
    export PATH=$PATH:/usr/local/gemini/bin
    cd /home/dnanexus
    gemini load --cores $(nproc) -v my.vcf $advanced_options my.gemini.db

    name=`dx describe "$vcfgz" --name`
    name="${name%.gz}"
    name="${name%.vcf}"

    dbgz=$(gzip -c my.gemini.db | dx upload --brief -o "${name}.gemini.db.gz" -)
    dx-jobutil-add-output dbgz "$dbgz" --class=file
}
