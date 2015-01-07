#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

int sort_func(const void *a, const void *b) {
    int32_t a1 = *((int32_t *)a);
    int32_t b1 = *((int32_t *)b);
    if(a1<b1) return -1;
    if(a1==b1) return 0;
    return 1;
}

int countElements(int32_t *tids, int l) {
    int total = 0, i;
    int32_t last = tids[0];
    for(i=1; i<l; i++) {
        if(last != tids[i]) {
            total++;
            last = tids[i];
        }
    }
    return last;
}

double countElements2(int32_t *tids, int pos, int l) {
    double out = 1.0;
    int i;
    int32_t ltid = tids[pos];
    for(i=pos+1; i<l; i++) {
        if(tids[i] != ltid) return out;
        out += 1.0;
    }
    return out;
}

void processStack(bam1_t **reads, int l, double threshold, bam_hdr_t *header, int rlen) {
    int32_t tids[l], *tid_id;
    int i, n_elements, pos = 0, n_above = 0, pos2;
    double *counts;
    char *p;

    for(i=0; i<l; i++) tids[i] = reads[i]->core.tid;
    qsort(tids, l, sizeof(int32_t), sort_func);
    n_elements = countElements(tids, l);
    tid_id = malloc(sizeof(int32_t)*n_elements);
    counts = malloc(sizeof(double)*n_elements);
    for(i=0; i<n_elements; i++) {
        counts[i] = countElements2(tids,pos,l)/((double) l);
        tid_id[i] = tids[pos];
        pos += counts[i];
        if(counts[i] > threshold) n_above++;
    }
    if(n_above == 1) {
        for(i=0; i<n_elements; i++) if(counts[i]>threshold) break;
        p = strtok(bam_get_qname(reads[0]), ":");
        printf("%s\t", p);
        p += strlen(p)+1;
        pos2 = atoi(p);
        printf("%i\t%i\t%s\n", ((pos2-rlen<0) ?0:pos2-rlen), pos2+rlen, header->target_name[tid_id[i]]);
    }
    free(tid_id);
    free(counts);
}

void usage(char *prog) {
    printf("%s [OPTIONS] input.bam\n", prog);
    printf("\n\
The goal of this program is to process alignments to repeat regions (LINE/L1 and\n\
SINE/Alu repeats, primarily) and determine if they might reasonably be assigned\n\
to a single repeat type. This is done by (1) seeing if the alignments from a\n\
single read are to only a single type of repeat and, if not, (2) seeing if only\n\
a small fraction are aligned to another repeat class (see the -f option).\n\
\n\
Options:\n\
\n\
-l INT   The read length, which is also how much we flank each insertion site (default 50bp).\n\
-f float The threshold for rejecting an assignment to a single repeat type. The\n\
         default is 0.01, so 1%%.\n");
}

int main(int argc, char *argv[]) {
    bam_hdr_t *header = NULL;
    bam1_t **reads = NULL, *read = NULL;
    htsFile *ifile = NULL;//, *obam= NULL;
    double threshold = 0.01;
    char *lname = NULL;
    int m_stack, l_stack = 0, i, rlen = 50;
    int rv = 0;

    if(argc==1) {
        usage(argv[0]);
        return 0;
    }
    for(i=1; i<argc; i++) {
        if(strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return 0;
        } else if(strcmp(argv[i], "-f") == 0) {
            threshold = atof(argv[++i]);
        } else if(strcmp(argv[i], "-l") == 0) {
            rlen = atoi(argv[++i]);
        } else if(ifile == NULL) {
            ifile = hts_open(argv[i], "rb");
            header = sam_hdr_read(ifile);
        } else {
            usage(argv[0]);
            return -1;
        }
    }

    if(ifile == NULL) {
        usage(argv[0]);
        return -1;
    }

    //Allocate the stack that'll hold the mappings
    m_stack = 100;
    reads = malloc(sizeof(bam1_t*) * m_stack);
    if(reads == NULL) {
        fprintf(stderr, "Couldn't allocate enough room for the 'reads' array. Out of memory?\n");
        rv =-2;
        goto err;
    }
    for(i=0; i<m_stack; i++) reads[i] = bam_init1();

    read = bam_init1();
    while(sam_read1(ifile, header, read) >= 0) {
        if(!lname) lname = strdup(bam_get_qname(read));
        if(read->core.flag & BAM_FUNMAP) continue;
        if(strcmp(lname, bam_get_qname(read)) == 0) {
            if(l_stack >= m_stack) {
                m_stack *= 2;
                if((reads = realloc(reads, sizeof(bam1_t *) * m_stack)) == NULL) {
                    fprintf(stderr, "Couldn't reallocate enough room to hold %i reads. Perhaps there's not enough memory!", m_stack);
                    rv = -1;
                    goto err;
                }
                for(i=l_stack; i<m_stack; i++) reads[i] = bam_init1();
                reads[l_stack] = bam_init1();
            }
            bam_copy1(reads[l_stack], read);
            l_stack++;
        } else {
            processStack(reads, l_stack, threshold, header, rlen);
            free(lname);
            lname = strdup(bam_get_qname(read));
            bam_copy1(reads[0], read);
            l_stack = 1;
        }
    }
    processStack(reads, l_stack, threshold, header, rlen);

err:
    for(i=0; i<m_stack; i++) bam_destroy1(reads[i]);
    free(reads);
    sam_close(ifile);
    return rv;
}
