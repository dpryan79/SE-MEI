INCLUDES = -Ihtslib
CC = gcc
OPTS = -Wall -O2 -g

.PHONY:	all htslib clean clean-all

.SUFFIXES:	.c .o 

.c.o:
	$(CC) -c $(OPTS) $(INCLUDES) $< -o $@

all: compactRepeats extractSoftclipped

htslib:
	$(MAKE) -C ../htslib

compactRepeats:	htslib compactRepeats.o
	$(CC) $(OPTS) $(INCLUDES) -o compactRepeats compactRepeats.o htslib/libhts.a -lz -lpthread

extractSoftclipped: htslib extractSoftclipped.o
	$(CC) $(OPTS) $(INCLUDES) -o extractSoftclipped extractSoftclipped.o htslib/libhts.a -lz -lpthread

clean:
	rm -f *.o compactRepeats extractSoftclipped

clean-all: clean
	make --directory=../htslib clean
