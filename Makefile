INCLUDES = -Ihtslib
CC = gcc
OPTS = -Wall -g

.PHONY:	all htslib clean clean-all

.SUFFIXES:	.c .o 

.c.o:
	$(CC) -c $(OPTS) $(INCLUDES) $< -o $@

all: compactRepeats extractSoftclipped compareGroups

htslib:
	$(MAKE) -C htslib

compareGroups:
	$(CC) $(OPTS) -o compareGroups compareGroups.c

compactRepeats:	htslib compactRepeats.o
	$(CC) $(OPTS) $(INCLUDES) -o compactRepeats compactRepeats.o htslib/libhts.a -lz -lpthread

extractSoftclipped: htslib extractSoftclipped.o
	$(CC) $(OPTS) $(INCLUDES) -o extractSoftclipped extractSoftclipped.o htslib/libhts.a -lz -lpthread

clean:
	rm -f *.o compactRepeats extractSoftclipped compareGroups

clean-all: clean
	make --directory=../htslib clean
