PROG:=main
SRCS:=main.cpp

CXX:=g++ -Wall -std=c++14 -O3 -Ofast

OBJS:=$(SRCS:.cpp=.o)
DEPS:=$(SRCS:.cpp=.d)

all: $(PROG)

$(PROG): $(OBJS)
	$(CXX) -o $@ $^

%.o: %.cc
	$(CXX) -MMD -c $<

.PHONY: edit go

go: $(PROG)
	$(PROG)

edit:
	@geany -s -i $(SRCS) *.h &

clean:
	@rm -f $(PROG) *.o *.gch *.d core tags

-include $(DEPS)
