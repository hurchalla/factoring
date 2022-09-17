
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifndef _MSC_VER
#  include <sys/time.h>
#endif
#include <stdint.h>
#include <math.h>

#include <type_traits>


#include "microecm_c.h"         // Interface to C version of microecm.
#include "get_single_factor.h"  // Convenience interface to C++ version of microecm.




#ifdef _MSC_VER
#ifndef _TIMES_H
#define _TIMES_H

#ifdef _WIN32
#include <sys/timeb.h>
#include <sys/types.h>
#include <winsock2.h>

int gettimeofday(struct timeval* t,void* timezone);

// from linux's sys/times.h

//#include <features.h>

#define __need_clock_t
#include <time.h>


/* Structure describing CPU time used by a process and its children.  */
struct tms
  {
    clock_t tms_utime;          /* User CPU time.  */
    clock_t tms_stime;          /* System CPU time.  */

    clock_t tms_cutime;         /* User CPU time of dead children.  */
    clock_t tms_cstime;         /* System CPU time of dead children.  */
  };

/* Store the CPU time used by this process and all its
   dead children (and their dead children) in BUFFER.
   Return the elapsed real time, or (clock_t) -1 for errors.
   All times are in CLK_TCKths of a second.  */
clock_t times (struct tms *__buffer);

typedef long long suseconds_t ;

#endif
#endif

int gettimeofday(struct timeval* t,void* timezone)
{       struct _timeb timebuffer;
        _ftime( &timebuffer );
        t->tv_sec=timebuffer.time;
        t->tv_usec=1000*timebuffer.millitm;
		return 0;
}

clock_t times (struct tms *__buffer) {

	__buffer->tms_utime = clock();
	__buffer->tms_stime = 0;
	__buffer->tms_cstime = 0;
	__buffer->tms_cutime = 0;
	return __buffer->tms_utime;
}
#endif   //_MSC_VER



double my_difftime(struct timeval * start, struct timeval * end)
{
  double secs;
  double usecs;

  if (start->tv_sec == end->tv_sec) {
    secs = 0;
    usecs = end->tv_usec - start->tv_usec;
  }
  else {
    usecs = 1000000 - start->tv_usec;
    secs = end->tv_sec - (start->tv_sec + 1);
    usecs += end->tv_usec;
    if (usecs >= 1000000) {
      usecs -= 1000000;
      secs += 1;
    }
  }

  return secs + usecs / 1000000.;
}




template <typename T>
T get_ecm_factor_C(T x)
{
    uint64_t loc_lcg = 0;
    int is_arbitrary = 0;
    return getfactor_uecm(x, is_arbitrary, &loc_lcg);
}
template <typename T>
T get_ecm_factor_CPP(T x)
{
    bool expect_arbitrary_factors = false;
    return hurchalla::get_single_factor_ecm(x, expect_arbitrary_factors);
}



int main(int argc, char** argv)
{
    FILE *in;
    double t_time;
    int i = 0, num, correct, avg_curves;
    struct timeval gstart;
    struct timeval gstop;
    int nf;
    int num_files;
    char filenames[100][80];
    static uint64_t knowns1[100000], knowns2[100000];

    // --- tinyecm.c can generate these files if you dont have them ---
    //
    // testing large semiprime inputs.
    {
        strcpy(filenames[i++], "semiprimes_42bit.dat"); // 0
        strcpy(filenames[i++], "semiprimes_44bit.dat"); // 1
        strcpy(filenames[i++], "semiprimes_46bit.dat"); // 2 
        strcpy(filenames[i++], "semiprimes_48bit.dat"); // 3 
        strcpy(filenames[i++], "semiprimes_50bit.dat"); // 4 
        strcpy(filenames[i++], "semiprimes_52bit.dat"); // 5 
        strcpy(filenames[i++], "semiprimes_54bit.dat"); // 6 
        strcpy(filenames[i++], "semiprimes_56bit.dat"); // 7 
        strcpy(filenames[i++], "semiprimes_58bit.dat"); // 8 
        strcpy(filenames[i++], "semiprimes_60bit.dat"); // 9 
        strcpy(filenames[i++], "semiprimes_62bit.dat"); // 10 
        strcpy(filenames[i++], "semiprimes_64bit.dat"); // 11 
        num_files = i;
    }


#ifndef UECM_SEMIPRIME_FILENUM
#  error "you must define UECM_SEMIPRIME_FILENUM"
#endif
    nf = UECM_SEMIPRIME_FILENUM;
//    assert(nf < num_files);

    char buf[1024];
    int curves;
    uint64_t known1, known2;
    uint64_t in64;

    // read in the first few test inputs to determine
    // the size of inputs in this file (we assume they
    // are all approximately the same).
    in = fopen(filenames[nf], "r");
    if (in == NULL)
    {
        printf("Error: can't open %s\n", filenames[nf]);
        return 1;
    }
    for (i = 0; i < 1; i++)
    {
        if (0 == fgets(buf, 1024, in)) {
            printf("error fgets1\n");
            return 1;
        }
        sscanf(buf, "%lu, %lu, %lu", &in64, &known1, &known2);
    }
    fclose(in);


    correct = 0;
    num = 100000;

    uint64_t inputs[100000];
    in = fopen(filenames[nf], "r");
    if (in == NULL)
    {
        printf("Error: can't open %s\n", filenames[nf]);
        return 1;
    }
    printf("testing file: %s\n", filenames[nf]);

    for (i = 0; i < num; i++)
    {
        uint64_t outf;

        if (0 == fgets(buf, 1024, in))  {
            printf("error fgets3\n");
            return 3;
        }
        sscanf(buf, "%lu, %lu, %lu", &in64, &known1, &known2);
        inputs[i] = (uint64_t)known1 * (uint64_t)known2;
        knowns1[i] = known1;
        knowns2[i] = known2;
    }
    fclose(in);
    for (int tries = 0; tries < 16; ++tries)
    {
#if 1
        // -------------------
        // benchmark C version
        // -------------------
        correct = 0;
        gettimeofday(&gstart, NULL);

        for (i = 0; i < num; i++)
        {
            uint64_t outf;
            {
                outf = get_ecm_factor_C(inputs[i]);
                if (outf == 0)
                    printf("error bad ecm factor1\n");
                if ((outf == knowns1[i]) ||
                    (outf == knowns2[i]))
                {
                    correct++;
                }
            }
        }

        gettimeofday(&gstop, NULL);

        t_time = my_difftime(&gstart, &gstop);

        printf("(C microecm) got %d of %d correct in %2.3f sec\n",
                correct, num, t_time);
//      printf("percent correct = %.2f\n", 
//                100.0*(double)correct / (double)num);
//      printf("average time per input = %.2f ms\n", 1000 * t_time / (double)num);
#endif

#if 1
        // ---------------------
        // benchmark C++ version
        // ---------------------
        correct = 0;
        gettimeofday(&gstart, NULL);

        for (i = 0; i < num; i++)
        {
            uint64_t outf;
            {
                outf = get_ecm_factor_CPP(inputs[i]);
                if (outf == 0)
                    printf("error bad ecm factor1\n");
                if ((outf == knowns1[i]) ||
                    (outf == knowns2[i]))
                {
                    correct++;
                }
            }
        }

        gettimeofday(&gstop, NULL);

        t_time = my_difftime(&gstart, &gstop);

        printf("(C++ port) got %d of %d correct in %2.3f sec\n",
                correct, num, t_time);
#endif
    }

    return 0;
}
