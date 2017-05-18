
HEADERS := $(shell ls *.hpp)
SOURCES := $(shell ls *.cpp | grep -v test.cpp)
OBJS := $(patsubst %.cpp, obj/%.o, $(SOURCES))

obj/%.o : %.cpp $(HEADERS)
	@mkdir -p obj
	g++ -c -o $@ $<

test: test.cpp $(OBJS)
	g++ -o $@ $^
