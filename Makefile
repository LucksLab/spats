
none : 
	@echo Target required
	@exit 1

TEST_PATH=$(shell pwd)/$(BINDIR):$(PATH)
.PHONY: test
test : $(TARGETS)
	@export PATH="$(TEST_PATH)"  &&  cd test/Read_Mapping  &&  bash test_read_mapping.sh


PYENV = PYTHONPATH=.
TOOLS_DIR = tools
TEST_PKG = spats.tests.test_spats

# unit tests

.PHONY: unit
unit :
	${PYENV} python -m unittest spats.tests.test_spats

# for some subset of tests:
# make u.[class]
# make u.[class].[function]
u.%:
	${PYENV} python -m unittest $(patsubst u.%, ${TEST_PKG}.%, $@)


# runs a method in tests/misc.py
t.%:
	@${PYENV} python $(patsubst t.%, ${TOOLS_DIR}/misc.py %, $@)


# profiling

# profiles a method in tests/misc.py
p.%:
	@mkdir -p tmp
	@PYTHONPATH=.:src python -m cProfile -o tmp/runprof.out $(patsubst p.%, ${TOOLS_DIR}/misc.py %, $@)

prof:
	@${PYENV} python ${TOOLS_DIR}/prof.py|head -50

pstats:
	@${PYENV} python -m pstats tmp/runprof.out

# use @profile on the function that you want line profiled
lp.%:
	@${PYENV} kernprof -l $(patsubst lp.%, ${TOOLS_DIR}/misc.py %, $@)

lprof:
	@${PYENV} python -m line_profiler misc.py.lprof

