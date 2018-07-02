/* bwtool_msum - simultaneously output same regions of multiple bigWigs */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <jkweb/common.h>
#include <jkweb/obscure.h>
#include <jkweb/linefile.h>
#include <jkweb/hash.h>
#include <jkweb/options.h>
#include <jkweb/sqlNum.h>
#include <jkweb/basicBed.h>
#include <jkweb/bigWig.h>
#include <beato/bigs.h>
#include "bwtool.h"
#include <beato/cluster.h>
#include "bwtool_shared.h"

void usage_msum()
/* Explain usage of msum program and exit. */
{
    errAbort(
      "bwtool msum - output sum of same regions over multiple bigWigs\n"
      "usage:\n"
      "   bwtool msum input1.bw input2.bw input3.bw ...\n"
      );
    }

void output_pbws_sums(struct perBaseWig *pbw_list, FILE *out)
/* outputs sums across of perBaseWigs all at the same section */
{
    if (pbw_list && pbw_list->len > 0)
    {
        int i = 0;
        struct perBaseWig *pbw;
        struct slDouble *c;
        //per base loop, assumes bedGraph output
        unsigned long prev_sum = 0;
        unsigned long prev_start = 0;
        for (i = 0; i < pbw_list->len; i++)
        {
            double sum = 0.0;
            for (pbw = pbw_list; pbw != NULL; pbw = pbw->next)
                sum += pbw->data[i];
            //condense bases with same sum
            if(sum != prev_sum) 
            {
                if(i > 0)
                {
                    print_str(pbw_list->chrom, out);
                    putc_unlocked('\t', out);
                    print_double(pbw_list->chromStart+prev_start, out);
                    putc_unlocked('\t', out);
                    print_double(pbw_list->chromStart+i, out);
                    putc_unlocked('\t', out);
                    if(prev_sum < 10)
                        putc_unlocked(prev_sum+'0', out);
                    else
                        print_double(prev_sum, out);
                    putc_unlocked('\n',out);
                }
                prev_start = i;
                prev_sum = sum;
            }
        }
        if(i > 0)
        {
            print_str(pbw_list->chrom, out);
            putc_unlocked('\t', out);
            print_double(pbw_list->chromStart+prev_start, out);
            putc_unlocked('\t', out);
            print_double(pbw_list->chromStart+i, out);
            putc_unlocked('\t', out);
            if(prev_sum < 10)
                putc_unlocked(prev_sum+'0', out);
            else
                print_double(prev_sum, out);
            putc_unlocked('\n',out);
        }
    }
}

void bwtool_msum(struct hash *options, char *regions, double fill,
        struct slName **p_files, char *tmp_dir, char *output_file)
/* bwtool_msum - main for msum program */
{
    struct metaBig *mb;
    struct metaBig *mb_list = NULL;
    struct bed *bed;
    struct slName *file;
    int num_sections = 0;
    int i = 0;
    boolean verbose = (hashFindVal(options, "verbose") != NULL) ? TRUE : FALSE;
    struct slName *labels = NULL;
    struct slName *files = *p_files;
    FILE *out = (output_file) ? mustOpen(output_file, "w") : stdout;
    /*char wigfile[512];
    safef(wigfile, sizeof(wigfile), "%s.tmp.wig", output_file);
    FILE *out = mustOpen(wigfile, "w");*/

    /* open the files one by one */
    if (slCount(files) == 1)
        check_for_list_files(&files, &labels, 0);
    for (file = files; file != NULL; file = file->next)
    {
        mb = metaBigOpenWithTmpDir(file->name, tmp_dir, regions);
        slAddHead(&mb_list, mb);
    }
    slReverse(&mb_list);
    num_sections = slCount(mb_list->sections);
    for (bed = mb_list->sections; bed != NULL; bed = bed->next)
    {
        struct perBaseWig *pbw_list = NULL;
        /* load each region */
        if (verbose)
            fprintf(stderr, "section %d / %d: %s:%d-%d\n", i++, num_sections, bed->chrom, bed->chromStart, bed->chromEnd);
        for (mb = mb_list; mb != NULL; mb = mb->next)
        {
            struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, bed->chrom, bed->chromStart, bed->chromEnd, FALSE, fill);
            /* if the load returns null then NA the whole thing. */
            /* this isn't very efficient but it's the easy way out. */
            if (!pbw)
                pbw = alloc_perBaseWig(bed->chrom, bed->chromStart, bed->chromEnd);
            slAddHead(&pbw_list, pbw);
        }
        slReverse(&pbw_list);
        output_pbws_sums(pbw_list, out);
        perBaseWigFreeList(&pbw_list);
    }
    /* close the files */
    carefulClose(&out);
    //now convert output to BigWig
    //writeBw(wigfile, output_file, mb->chromSizeHash);
    while ((mb = slPopHead(&mb_list)) != NULL)
        metaBigClose(&mb);
    slNameFreeList(p_files);
}
