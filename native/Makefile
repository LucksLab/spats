
none : 
	@echo Target required
	@exit 1


HEADERS := $(shell ls *.hpp)
SOURCES := $(shell ls *.cpp | grep -v test.cpp | grep -v cotrans.cpp)
OBJS := $(patsubst %.cpp, obj/%.o, $(SOURCES))
#CFLAGS := -g
CFLAGS := -Ofast

obj/%.o : %.cpp $(HEADERS)
	@mkdir -p obj
	g++ -c $(CFLAGS) -o $@ $<

bin/%: %.cpp $(OBJS)
	@mkdir -p bin
	g++ $(CFLAGS) -lsqlite3 -lpthread -o $@ $^

debug: bin/test
	lldb test

all: bin/test bin/unittest bin/cotrans

.PHONY: clean
clean:
	rm -rf obj *~

.PHONY: clear
clear: clean
	rm -rf bin
