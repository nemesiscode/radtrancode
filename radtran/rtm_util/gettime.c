/*
 *
 *_TITL:	gettime
 *
 *_DESC:	Obtain the system time for use in determining the total time
 *		used by the program since execution.
 *
 *_HIST:	PDP	8oct2002	Original version. The incentive of
 *					creating gettime.c was to replace 
 *					the "etime" command (used for
 *					determining the system time and
 *					hence the total program-execution 
 *					time) which isn't recognised by
 *					the Intell FORTRAN compiler for
 *					use on Linux 7.1+
 *					systems. However, the routines
 *					within this code work for both
 *					Linux and OSF systems.

 *
 */

/* Fortran interface */
#include "f2c.h"

/* Date and time header */
#include <time.h>

/* Standard headers */
#include <string.h>


/* gettime_ function */
void gettime_(double *tme1)
{
  /* time_t is a long */
  time_t tp;

  /* return the time; & = address */
  time(&tp);
  /* explicitly-cast tp as a double */
  *tme1 = (double) tp;
} 

/* PROGRAM FOR TESTING PURPOSES ...
 *
 * #include <stdio.h>
 * #include <stdlib.h>
 *
 * /* gettime_ function *
 * main()
 * {
 *   int diff;
 *   time_t tp, tme1, tme2;
 *
 *   /* return the first time; & = address *
 *   time(&tp);
 *   tme1 = tp;
 *   printf("first time: %d\n", tme1);
 *
 *   /* sleep for ten seconds*
 *   sleep(10);
 *
 *   /* return the second time *
 *   time(&tp);
 *   tme2 = tp;
 *   printf("second time: %d\n", tme2);
 *
 *   /* determine the difference between tme2 - tme1 *
 *   diff= difftime(tme2, tme1);  
 *   printf("tme1, tme2: %d\n%d", tme1, tme2);
 *   printf("time: %d\n", diff);
 * }
 */
 
/********************************************************************/
