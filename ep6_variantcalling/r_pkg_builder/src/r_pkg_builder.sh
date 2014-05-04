#!/bin/bash

main() {

    set -x

    if [ -n "${apt[@]}" ]; then
        DEBIAN_FRONTEND=NONINTERACTIVE apt-get install --no-install-suggests --no-install-recommends -y ${apt[@]}
    fi
    
    touch /cran.txt
    for cran_pkg in "${cran[@]}"; do
        echo $cran_pkg >> /cran.txt
    done
    touch /bioconductor.txt
    for bioc_pkg in "${bioconductor[@]}"; do
        echo $bioc_pkg >> /bioconductor.txt
    done

    R --no-save --quiet < /libs.R
    tar -C / -cf /Rlibs.tar /usr/local/lib/R/site-library
    gzip /Rlibs.tar
    output_Rlibs=$(dx upload /Rlibs.tar.gz --brief --destination "$output_name")
    dx-jobutil-add-output Rlibs "$output_Rlibs" --class=file
}
