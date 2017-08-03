
HEADERS := $(shell ls *.hpp)
SOURCES := $(shell ls *.cpp | grep -v test.cpp)
OBJS := $(patsubst %.cpp, obj/%.o, $(SOURCES))
#CFLAGS := -g
CFLAGS := -Ofast

obj/%.o : %.cpp $(HEADERS)
	@mkdir -p obj
	g++ -c $(CFLAGS) -o $@ $<

test: test.cpp $(OBJS)
	g++ $(CFLAGS) -lsqlite3 -lpthread -o $@ $^

debug: test
	lldb test
