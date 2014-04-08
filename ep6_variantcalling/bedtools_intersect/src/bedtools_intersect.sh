#!/bin/bash

main() {
	aname=$(dx describe "$a" --name)

	case "$aname" in
	*.gz)
			pp="gunzip -c"
			auncomp="${aname%.gz}"
			abase="${auncomp%.*}"
			if [ -z "$abase" ]; then
				abase="$auncomp"
				aext=".gz"
			else
				aext=".${auncomp##*.}.gz"
			fi
			aparam="-a"
			;;
    *.bam) 
            pp="gunzip -c"
            abase="${aname%.bam}"
            aext=".bam"
            aparam="-a"
            ;;
    *)
            pp="tee /dev/null"
            abase="${aname%.*}"
            if [ -z "$abase" ]; then
				abase="$aname"
				aext=".gz"
			else
				aext=".${abase##*.}.gz"
			fi
			aparam="-abam"
            ;;
    esac

    dx download -o b "$b"

    res=$(dx cat "$a" | $pp | bedtools intersect $aparam stdin -b b $advanced_options \
    		| bgzip -c | dx upload -o "${abase}${filename_insert}${aext}" --brief -)

    dx-jobutil-add-output result "$res" --class=file
}
