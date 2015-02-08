all:
	make -C upstream/seqtk -j `nproc`
	make -C upstream/bfc -j `nproc`
	make -C upstream/ropebwt2 -j `nproc`
	make -C upstream/fermi2 -j `nproc`
	mkdir -p resources/usr/local/bin
	cp upstream/seqtk/trimadap upstream/bfc/bfc upstream/ropebwt2/ropebwt2 upstream/fermi2/fermi2 upstream/fermi2/fermi2.pl resources/usr/local/bin

.PHONY: all
