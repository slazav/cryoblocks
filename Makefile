TARGETS=cryoblocks tests

LDLIBS= -lm -lgsl -lhe3
CPPFLAGS=-std=c++11 -DHE3 -g
CC=g++

all: $(TARGETS)
	./tests

clean:
	rm -f $(TARGETS) *.o

cryoblocks.o tests.o: c_args.h c_links.h c_blocks.h c_dim.h
c_dim.o c_args.o cryoblocks.o tests.o: inc/err.h
tests.o: inc/assert_err.h

cryoblocks: cryoblocks.o c_dim.o c_args.o

tests: tests.o c_dim.o c_args.o


