#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <getopt.h>

typedef struct {
    char *chrom;
    int start, stop;
} BED;

FILE ** parseList(char *s, int *l) {
    FILE **fps;
    char *p = s;
    int i = 0;

    //Count the number of files
    *l = 1;
    while(*p) {
        if(*p==',') (*l)++;
        p++;
    }
    fps = malloc(sizeof(FILE*) * (*l));
    assert(fps);

    //Fill in the file pointer list
    p = strtok(s, ",");
    while(p) {
        fps[i] = NULL;
        fps[i] = fopen(p, "r");
        if(!fps[i++]) {
            fprintf(stderr, "Error: couldn't open %s!\n", p);
            assert(1==0);
        }
        p = strtok(NULL, ",");
    }

    return fps;
}

BED * initBED() {
    BED *b = malloc(sizeof(BED));
    assert(b);
    b->chrom = NULL;
    return b;
}

//On the first call, b->chrom must be NULL!
BED * parseLine(char *s, BED *b) {
    char *p, *p2;

    if(b->chrom) free(b->chrom);

    p = strtok(s, "\t");
    b->chrom = strdup(p);
    p = strtok(NULL, "\t");
    b->start = atoi(p);
    p = strtok(NULL, "\t");
    p2 = p;
    while(*p2) {
        if(iscntrl(*p2)) *p2 = '\0'; //Remove a possible line ending
        p2++;
    }
    b->stop = atoi(p);
    return b;
}

BED * destroyBED(BED *b) {
    if(b) {
        if(b->chrom) free(b->chrom);
        free(b);
    }
    return NULL;
}

//Returns NULL and free()s b on EOF
BED * nextLine(FILE *fp, BED *b) {
    char line[1024];
    if(!(fgets(line, 1024, fp))) return destroyBED(b);
    return parseLine(line, b);
}

//Returns 1 for in group and 0 otherwise
int inGroup(FILE **fp, BED **bs, int l, BED *ref) {
    int i, rv = 0;

    for(i=0; i<l; i++) {
        //Ensure we're on the same chromosome
        while(bs[i] && strcmp(bs[i]->chrom, ref->chrom) < 0) {
            bs[i] = nextLine(fp[i], bs[i]);
        }
        if(bs[i] && strcmp(bs[i]->chrom, ref->chrom) != 0) continue;

        //Looking at overlapping positions?
        while(bs[i] && bs[i]->start < ref->start) {
            bs[i] = nextLine(fp[i], bs[i]);
            if(!bs[i]) break;
            if(strcmp(bs[i]->chrom, ref->chrom) != 0) break; //Otherwise, we jump to EOF!?!?
        }
        if(bs[i] && strcmp(bs[i]->chrom, ref->chrom) == 0 && 
            ((bs[i]->start >= ref->start && bs[i]->start < ref->stop) ||
            (bs[i]->stop > ref->start && bs[i]->stop <= ref->stop) ||
            (bs[i]->start < ref->start && bs[i]->stop > ref->stop)))
            rv = 1;
    }
    return rv;
}

void parseLabels(char *s, char **l1, char **l2) {
    char *p = s;
    *l1 = NULL;
    *l2 = NULL;
    while(*p) {
        if(*p == ',') {
            *l2 = strdup(p+1);
            *p = '\0';
            *l1 = strdup(s);
            *p = ',';
            break;
        }
        p++;
    }
    assert(*l1);
    assert(*l2);
}

void usage(char *prog) {
    fprintf(stderr, "Usage: %s [OPTIONS] group1,list group2,list all_sites > output\n", prog);
    fprintf(stderr, "\n\
This program parses two comma separated lists of BED files containing putative\n\
insert sites and looks for sites unique to a single group (whether it be in a\n\
single sample or multiple samples). The all_sites BED file contains all sites\n\
described be all samples (this is used mainly for efficiency).\n\
\n\
N.B., there should be no spaces after commas in the lists (e.g.,\n\
\"sample1.bed, sample2.bed\" would be incorrect, but \"sample1.bed,sample2.bed\"\n\
would work correctly). Note also that all BED files must be lexographically\n\
sorted.\n\
\n\
Options\n\
\n\
-L STR,STR    A comma-separated list of group labels (e.g.,\n\
              \"Control,Treatment\"). If this is not specified, the output will\n\
              use \"Group1\" and \"Group2\".\n");
}

int main(int argc, char *argv[]) {
    FILE **fp1 = NULL, **fp2 = NULL;
    FILE *ref;
    int nfile1 = 0, nfile2 = 0, i;
    int inG1, inG2;
    char c, *labels[2] = {NULL, NULL};
    BED **BED1, **BED2, *refBED;

    while((c = getopt(argc, argv, "hL:")) != -1) {
        switch(c) {
        case 'h' :
            usage(argv[0]);
            return 0;
        case 'L' :
            parseLabels(optarg, labels, labels+1);
            break;
        default :
            fprintf(stderr, "Invalid option '%i'\n", c);
            usage(argv[0]);
            return 1;
        }
    }

    if(argc-optind != 3) {
        usage(argv[0]);
        return 1;
    }
    fp1 = parseList(argv[optind++], &nfile1);
    fp2 = parseList(argv[optind++], &nfile2);
    ref = fopen(argv[optind], "r");

    //In case we didn't specify -L
    if(!labels[0]) labels[0] = strdup("Group1");
    if(!labels[1]) labels[1] = strdup("Group2");

    BED1 = malloc(sizeof(BED*) * nfile1);
    BED2 = malloc(sizeof(BED*) * nfile2);
    assert(BED1);
    assert(BED2);
    refBED = initBED();
    for(i=0; i<nfile1; i++) {
        BED1[i] = initBED();
        BED1[i] = nextLine(fp1[i], BED1[i]);
    }
    for(i=0; i<nfile2; i++) {
        BED2[i] = initBED();
        BED2[i] = nextLine(fp2[i], BED2[i]);
    }

    while((refBED = nextLine(ref, refBED))) {
        inG1 = inGroup(fp1, BED1, nfile1, refBED);
        inG2 = inGroup(fp2, BED2, nfile2, refBED);
        if(inG1 && !inG2)
            printf("%s\t%i\t%i\t%s\n", refBED->chrom, refBED->start, refBED->stop, labels[0]);
        if(!inG1 && inG2)
            printf("%s\t%i\t%i\t%s\n", refBED->chrom, refBED->start, refBED->stop, labels[1]);
    }

    for(i=0; i<nfile1; i++) {
        BED1[i] = destroyBED(BED1[i]);
        fclose(fp1[i]);
    }
    free(BED1);
    free(fp1);
    for(i=0; i<nfile2; i++) {
        BED2[i] = destroyBED(BED2[i]);
        fclose(fp2[i]);
    }
    free(BED2);
    free(fp2);
    free(labels[0]);
    free(labels[1]);
    fclose(ref);
    return 0;
}
