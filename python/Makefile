
.PHONY: install clean test

install:
	python setup.py install

clean:
	find . -name "*.pyc" -exec rm {} \;
	rm -rf nwalign/resources
	rm -rf dist
	rm -rf build
	rm -rf *.egg-info

test:
	nosetests test/test_nwalign.py

js:
	python -c "from nwalign import blosum; blosum.to_js('blosum.js')"
