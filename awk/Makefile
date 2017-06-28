
.PHONY: all test

all:
	rm -f nwalign
	ln -s nwalign.awk nwalign

test:
	export PATH=$$PWD:$$PATH ; cram nwalign_awk.t
