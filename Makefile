
none : 
	@echo Target required
	@exit 1


PYENV = PYTHONPATH=.
TOOLS_DIR = tools
PKG_NAME = spats_shape_seq
TEST_PKG = ${PKG_NAME}.tests

# unit tests

.PHONY: unit
unit :
	nosetests ${PKG_NAME}
	@if [ ! -f native/bin/unittest ] ; then echo "\n**Warning: native tests not executed.\n"; fi

.PHONY: test
test : unit

VERSION = $(shell python -c 'import spats_shape_seq; print spats_shape_seq._VERSION')
IS_PRODUCTION = $(shell python -c 'import spats_shape_seq; print "1" if spats_shape_seq._PRODUCTION else "0"')
IS_PUBLIC_RELEASE = $(shell python -c 'import spats_shape_seq; print "1" if spats_shape_seq._PUBLIC_RELEASE else "0"')
.PHONY: show_ver
show_ver:
	@echo -n Version: ${VERSION}
ifeq ($(IS_PUBLIC_RELEASE), 1)
	@echo -n \ Public
else
	@echo -n \ Private
endif
ifeq ($(IS_PRODUCTION), 1)
	@echo \ Production
else
	@echo \ Beta
endif

pip_dist: unit
	@read -p "Submit new release, version ${VERSION}? " ANS; if [ "$$ANS" == "y" ]; then echo "Submitting..."; else echo "Aborted."; exit 1; fi
	@rm -rf build dist
	python setup.py sdist
	twine upload dist/${PKG_NAME}-${VERSION}.tar.gz

local_install:
	@(pip show -q ${PKG_NAME} && sudo pip uninstall -y -q ${PKG_NAME}) || echo
	@rm -rf build dist
	@python setup.py sdist
	@sudo pip install dist/${PKG_NAME}-${VERSION}.tar.gz

local_uninstall:
	@sudo pip uninstall -y -q ${PKG_NAME}

.PHONY: docs
docs:
	@cd doc  &&  make html  ||  (echo "Docs failed. Make sure sphinx is installed on your system:\n\n$ sudo pip install sphinx\n"  &&  exit 1)

.PHONY: showdocs
showdocs: docs
	@open doc/build/html/index.html

.PHONY: native
native:
	@cd native  &&  make all

cjb: pkg/cjb.zip
	@unzip -q pkg/cjb.zip

bin/UIClient.app: pkg/UIClient.zip
	@mkdir -p bin
	@rm -rf bin/UIClient.app
	@cd bin && unzip -q ../pkg/UIClient.zip
	@touch bin/UIClient.app

.PHONY: vizprep
vizprep: clean cjb bin/UIClient.app

.PHONY: viz
viz:
	@PYTHONPATH=. python ${PKG_NAME}/tool.py viz

.PHONY: clean
clean:
	@rm -rf cjb
	@rm -rf bin/UIClient.app

.PHONY: vizsrv
vizsrv: cjb
	@PYTHONPATH=. python ${TOOLS_DIR}/runviz.py

.PHONY: cjb-pkg
cjb-pkg:
	@mkdir -p tmp
	@rm -rf tmp/cjb
	@cp -rf ~/fw/trees/jbpy/cjb tmp/
	@rm `find tmp/cjb -name "*~"`
	@rm `find tmp/cjb -name "*.pyc"`
	@cd tmp  &&  zip -r cjb.zip cjb
	@mv tmp/cjb.zip pkg/cjb.zip

.PHONY: UIC-pkg
UIC-pkg:
	@mkdir -p tmp
	@rm -rf tmp/UIClient.app ~/fw/trees/jbpy/UIClient/tmp/UIC.xcarchive
	@cd ~/fw/trees/jbpy/UIClient && xcodebuild -scheme UIClient archive -archivePath tmp/UIC.xcarchive
	@mv ~/fw/trees/jbpy/UIClient/tmp/UIC.xcarchive/Products/Applications/UIClient.app tmp/
	@cd tmp && zip -r UIClient.zip UIClient.app
	@mv tmp/UIClient.zip pkg/UIClient.zip


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
