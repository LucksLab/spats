
HEADERS := $(shell ls *.hpp)
SOURCES := $(shell ls *.cpp | grep -v test.cpp)
OBJS := $(patsubst %.cpp, obj/%.o, $(SOURCES))
#CFLAGS := -g
CFLAGS := -O3

obj/%.o : %.cpp $(HEADERS)
	@mkdir -p obj
	g++ -c $(CFLAGS) -o $@ $<

test: test.cpp $(OBJS)
	g++ $(CFLAGS) -o $@ $^

debug: test
	lldb test
