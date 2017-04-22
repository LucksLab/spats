
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
	rm -rf $(BINDIR)

TEST_PATH=$(shell pwd)/$(BINDIR):$(PATH)
.PHONY: test
test : $(TARGETS)
	@export PATH="$(TEST_PATH)"  &&  cd test/Read_Mapping  &&  bash test_read_mapping.sh

# unit tests

.PHONY: unit
unit :
	PYTHONPATH=.:src python -m unittest tests.test_spats

# for some subset of tests:
# make u.[class]
# make u.[class].[function]
u.%:
	@PYTHONPATH=.:src python -m unittest $(patsubst u.%, tests.test_spats.%, $@)

