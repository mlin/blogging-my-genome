#!/bin/bash
# gemini_console 0.0.1

set -e -o pipefail

setup_gemini() {
    wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
    mkdir -p /usr/local /usr/local/share/gemini
    python gemini_install.py /usr/local /usr/local/share/gemini
}

main() {
    setup_gemini &
    dx cat "$dbgz" | zcat > my.db
    wait
    export PATH=$PATH:/usr/local/gemini/bin

    mkdir gemini_console_outputs

    echo "Your GEMINI database has been downloaded to /home/dnanexus/my.db
Sleeping for $minutes minutes. Now, 'dx ssh $DX_JOB_ID' and try:
gemini query -q \"select count(*) from variants where type='snp' and is_lof=1 and (aaf_1kg_all<=0.01 or aaf_1kg_all is null)\" my.db
Furthermore, any files you place in the /home/dnanexus/gemini_console_outputs/
directory will be output from the job upon its completion ('killall sleep' to
complete the job immediately.)"

    sleep $(expr $minutes \* 60) || true

    for fn in gemini_console_outputs/*; do
        dx-jobutil-add-output --class file --array gemini_console_outputs "$(dx upload --brief $fn)"
    done
}
