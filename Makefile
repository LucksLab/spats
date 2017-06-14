
none : 
	@echo Target required
	@exit 1


PYENV = PYTHONPATH=.
TOOLS_DIR = tools
TEST_PKG = spats_shape_seq.tests

# unit tests

.PHONY: unit
unit :
	nosetests

.PHONY: test
test : unit

VERSION = $(shell grep version= setup.py | sed -e 's/.*="\(.*\)",.*/\1/g')
pip_dist: unit
	@read -p "Submit new release, version ${VERSION}? " ANS; if [ "$$ANS" == "y" ]; then echo "Submitting..."; else echo "Aborted."; exit 1; fi
	@rm -rf build dist
	python setup.py sdist
	twine upload dist/spats_shape_seq-${VERSION}.tar.gz

.PHONY: docs
docs:
	@cd doc  &&  make html

cjb: pkg/cjb.zip
	@unzip pkg/cjb.zip

bin/UIClient.app: pkg/UIClient.zip
	@mkdir -p bin
	@rm -rf bin/UIClient.app
	@cd bin && unzip ../pkg/UIClient.zip
	@touch bin/UIClient.app

.PHONY: viz
viz: cjb bin/UIClient.app
	@./bin/UIClient.app/Contents/MacOS/UIClient &
	@PYTHONPATH=.:. python ${TOOLS_DIR}/runviz.py

.PHONY: jbpy-pkg
jbpy-pkg:
	@mkdir -p tmp
	@rm -rf tmp/cjb
	@cp -rf ~/fw/trees/jbpy/cjb tmp/
	@rm `find tmp/cjb -name "*~"`
	@rm `find tmp/cjb -name "*.pyc"`
	@cd tmp  &&  zip -r cjb.zip cjb
	@mv tmp/cjb.zip pkg/cjb.zip

# for some subset of tests:
# make u.[module]
# make u.[module]:[class]
# make u.[module]:[class].[function]
u.%:
	nosetests $(patsubst u.%, ${TEST_PKG}.%, $@)


# runs a method in tests/misc.py
t.%:
	@${PYENV} python $(patsubst t.%, ${TOOLS_DIR}/misc.py %, $@)


# profiling

# profiles a method in tests/misc.py
p.%:
	@mkdir -p tmp
	@${PYENV} python -m cProfile -o tmp/runprof.out $(patsubst p.%, ${TOOLS_DIR}/misc.py %, $@)

prof:
	@${PYENV} python ${TOOLS_DIR}/prof.py|head -50

pstats:
	@${PYENV} python -m pstats tmp/runprof.out

# use @profile on the function that you want line profiled
lp.%:
	@${PYENV} kernprof -l $(patsubst lp.%, ${TOOLS_DIR}/misc.py %, $@)

lprof:
	@${PYENV} python -m line_profiler misc.py.lprof
