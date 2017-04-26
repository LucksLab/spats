
SCRIPTDIR = src
BINDIR = bin
CFGDIR = config

SUPPORT_SCRIPT_NAMES = analyze_spats_targets targets_analyzer spats_common
SUPPORT_SCRIPT_SRCS = $(patsubst %, $(SCRIPTDIR)/%.py, $(SUPPORT_SCRIPT_NAMES))
SUPPORT_SCRIPTS = $(patsubst %, $(BINDIR)/%.py, $(SUPPORT_SCRIPT_NAMES))
TARGET_SCRIPT_NAMES = adapter_trimmer reactivities_split spats
TARGET_SCRIPT_SRCS = $(patsubst %, $(SCRIPTDIR)/%.py, $(TARGET_SCRIPT_NAMES))
TARGET_SCRIPTS = $(patsubst %, $(BINDIR)/%, $(TARGET_SCRIPT_NAMES))

TARGETS = $(TARGET_EXES) $(TARGET_SCRIPTS) $(SUPPORT_SCRIPTS)


all : $(TARGETS)

$(BINDIR)/% : $(SCRIPTDIR)/%.py
	@mkdir -p $(BINDIR)
	cp $< $@
	@chmod 755 $@

$(BINDIR)/%.py : $(SCRIPTDIR)/%.py
	@mkdir -p $(BINDIR)
	cp $< $@

.PHONY: clean
clean :


.PHONY: clear
clear : clean
	rm -rf $(BINDIR) *.lprof tmp

TEST_PATH=$(shell pwd)/$(BINDIR):$(PATH)
.PHONY: test
test : $(TARGETS)
	@export PATH="$(TEST_PATH)"  &&  cd test/Read_Mapping  &&  bash test_read_mapping.sh

PYENV = PYTHONPATH=.:${SCRIPTDIR}:cpp/build/lib.macosx-10.11-intel-2.7

# unit tests

.PHONY: unit
unit :
	@${PYENV} python -m unittest tests.test_spats

# for some subset of tests:
# make u.[class]
# make u.[class].[function]
u.%:
	@${PYENV} python -m unittest $(patsubst u.%, tests.test_spats.%, $@)

# runs a method in tests/misc.py
t.%:
	@${PYENV} python $(patsubst t.%, tests/misc.py %, $@)


# profiling

# profiles a method in tests/misc.py
p.%:
	@mkdir -p tmp
	@PYTHONPATH=.:src python -m cProfile -o tmp/runprof.out $(patsubst p.%, tests/misc.py %, $@)

prof:
	@${PYENV} python tests/prof.py|head -50

pstats:
	@${PYENV} python -m pstats tmp/runprof.out

# use @profile on the function that you want line profiled
lp.%:
	@${PYENV} kernprof -l $(patsubst lp.%, tests/misc.py %, $@)

lprof:
	@${PYENV} python -m line_profiler misc.py.lprof

