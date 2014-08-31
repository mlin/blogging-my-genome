#!/bin/bash
# gemini_load 0.0.1

main() {
    echo "Value of vcfgz: '$vcfgz'"
    echo "Value of advanced_options: '$advanced_options'"

    # TODO: install & run VEP according to
    # http://gemini.readthedocs.org/en/latest/content/functional_annotation.html#stepwise-installation-and-usage-of-vep

    wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
    mkdir -p /usr/local /usr/local/share/gemini
    python gemini_install.py /usr/local /usr/local/share/gemini
    export PATH=$PATH:/usr/local/gemini/bin

    gemini update --dataonly --extra cadd_score
    gemini update --dataonly --extra gerp_bp

    dx cat "$vcfgz" | zcat > my.vcf
    gemini load --cores $(nproc) -v my.vcf $advanced_options my.gemini.db

    name=`dx describe "$vcfgz" --name`
    name="${name%.gz}"
    name="${name%.vcf}"

    dbgz=$(gzip -c my.gemini.db | dx upload --brief -o "${name}.gemini.db.gz" -)
    dx-jobutil-add-output dbgz "$dbgz" --class=file
}
