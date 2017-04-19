
CC = g++
AR = ar
RANLIB = ranlib

SRCDIR = src
SCRIPTDIR = src
DEPSDIR = deps
OBJDIR = obj
BINDIR = bin
DISTDIR = dist
CFGDIR = config

VERSION = $(shell cat $(CFGDIR)/version)
BASE_CFLAGS = -I$(DEPSDIR) -DPACKAGE_VERSION=\"$(VERSION)\" -Wall -Wno-strict-aliasing -m64
RELEASE_FLAGS = -O3 -DNDEBUG
DEBUG_FLAGS = -g -O0 -DDEBUG
CFLAGS = $(BASE_CFLAGS) $(RELEASE_FLAGS)
LFLAGS = $(CFLAGS) -static

SRCS = $(shell echo $(SRCDIR)/*.cpp)
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))
PAREN:=(
LIB_SRCS = $(shell grep -L "main$(PAREN)" $(SRCDIR)/*.cpp)
LIB_OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(LIB_SRCS))
LIB = $(OBJDIR)/libspats.a
TARGET_SRCS = $(shell grep -l "main$(PAREN)" $(SRCDIR)/*.cpp)
TARGET_EXES = $(patsubst $(SRCDIR)/%.cpp, $(BINDIR)/%, $(TARGET_SRCS))

SUPPORT_SCRIPT_NAMES = analyze_spats_targets targets_analyzer spats_common
SUPPORT_SCRIPT_SRCS = $(patsubst %, $(SCRIPTDIR)/%.py, $(SUPPORT_SCRIPT_NAMES))
SUPPORT_SCRIPTS = $(patsubst %, $(BINDIR)/%.py, $(SUPPORT_SCRIPT_NAMES))
TARGET_SCRIPT_NAMES = adapter_trimmer reactivities_split spats
TARGET_SCRIPT_SRCS = $(patsubst %, $(SCRIPTDIR)/%.py, $(TARGET_SCRIPT_NAMES))
TARGET_SCRIPTS = $(patsubst %, $(BINDIR)/%, $(TARGET_SCRIPT_NAMES))

TARGETS = $(TARGET_EXES) $(TARGET_SCRIPTS) $(SUPPORT_SCRIPTS)


all : $(TARGETS)

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@mkdir -p $(OBJDIR)
	@echo CC: $<
	$(CC) $(CFLAGS) -c -o $@ $<

$(LIB) : $(LIB_OBJS)
	@mkdir -p $(BINDIR)
	@echo AR: $@
	@rm -f $@
	$(AR) cr $@ $(LIB_OBJS)
	$(RANLIB) $@

$(BINDIR)/adapter_trimmer : $(SCRIPTDIR)/adapter_trimmer.py
	@mkdir -p $(BINDIR)
	cp $< $@
	@chmod 755 $@

$(BINDIR)/reactivities_split : $(SCRIPTDIR)/reactivities_split.py
	@mkdir -p $(BINDIR)
	cp $< $@
	@chmod 755 $@

$(BINDIR)/spats : $(SCRIPTDIR)/spats.py
	@mkdir -p $(BINDIR)
	cp $< $@
	@chmod 755 $@

$(BINDIR)/%.py : $(SCRIPTDIR)/%.py
	cp $< $@

$(BINDIR)/% : $(OBJDIR)/%.o $(LIB)
	@echo LINK: $(shell basename $@)
	$(CC) $(LFLAGS) -o $@ $?

.PHONY: clean
clean :
	rm -rf $(OBJDIR)

.PHONY: clear
clear : clean
	rm -rf $(DISTDIR) $(BINDIR)

TEST_PATH=$(shell pwd)/$(BINDIR):$(PATH)
.PHONY: test
test : $(TARGETS)
	@export PATH="$(TEST_PATH)"  &&  cd test/Read_Mapping  &&  bash test_read_mapping.sh

PROJECT = spats
OS = $(shell uname)
ARCH = $(shell uname -m)
FULLNAME = $(PROJECT)-$(VERSION)-$(OS)-$(ARCH)
PKG = $(FULLNAME).tgz
PKGDIR = $(PROJECT)-$(VERSION)

$(DISTDIR)/$(PKG) : $(TARGETS)
	@echo DIST: $(FULLNAME).tgz
	mkdir -p $(DISTDIR)/$(PKGDIR)
	cp $(TARGETS) $(DISTDIR)/$(PKGDIR)
	cd $(DISTDIR)  &&  tar czf $(PKG) $(PKGDIR)

dist: clear $(DISTDIR)/$(PKG)
