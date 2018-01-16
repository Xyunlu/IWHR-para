/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: md_timer_mpi.c,v $
 *
 * $Author: tuminaro $
 *
 * $Date: 1999/09/30 17:05:01 $
 *
 * $Revision: 5.2 $
 *
 * $Name:  $
 *====================================================================*/
#ifndef lint
static char *cvs_timerint_id =
  "$Id: md_timer_mpi.c,v 5.2 1999/09/30 17:05:01 tuminaro Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/
#ifdef MPI
#include <mpi.h>
#endif

/* MPI timer */

extern double second(void);
double second(void)

{

#ifdef MPI
  return (MPI_Wtime());
#else 
  return (0.0);
#endif

} /* second */
