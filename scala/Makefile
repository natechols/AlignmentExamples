
.PHONY: build clean test

build:
	mkdir -p src/main/resources
	cp ../data/BLOSUM* src/main/resources/
	sbt compile pack

clean:
	sbt clean
	rm -rf project
	rm -rf target
	rm -rf src/main/resources

test:
	mkdir -p src/main/resources
	cp ../data/BLOSUM* src/main/resources/
	sbt test
