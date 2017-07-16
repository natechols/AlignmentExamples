
.PHONY: build clean test

build:
	mkdir -p src/main/resources
	cp ../data/BLOSUM* src/main/resources/
	javac src/main/java/com/github/natechols/nwalign/NeedlemanWunsch.java

clean:
	rm -rf src/main/resources
	find . -name "*.class" | xargs rm -f
