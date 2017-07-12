
all:
	export GOPATH=${PWD} && go install nwalign

test:
	$(eval DATA := "`readlink -f ../data`")
	export BLOSUM_DATA=${DATA} && cd src/libnwalign && go test

clean:
	rm -rf bin pkg
