TARGETS=cryoblocks

LDLIBS= -lm -lgsl -lhe3
CPPFLAGS=-std=c++11 -DHE3 -g
CC=g++

all: $(TARGETS)
clean:
	rm -f $(TARGETS) *.o

cryoblocks.o: c_args.h c_links.h c_blocks.h c_dim.h c_err.h
c_dim.o: c_err.h

cryoblocks: cryoblocks.o c_dim.o c_args.o

