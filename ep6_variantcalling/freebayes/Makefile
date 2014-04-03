all: resources/usr/bin/freebayes resources/usr/bin/vcfcombine

resources/usr/bin/freebayes:
	cd upstream && git submodule update --init --recursive && make
	cp upstream/bin/freebayes resources/usr/bin

resources/usr/bin/vcfcombine: resources/usr/bin/freebayes
	make -C upstream/vcflib vcfcombine
	cp upstream/vcflib/bin/vcfcombine resources/usr/bin

clean:
	cd upstream && make clean || true
	cd upstream/vcflib && make clean || true

.PHONY: all freebayes vcfcombine clean
