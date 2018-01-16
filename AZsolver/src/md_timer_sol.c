/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: md_timer_sol.c,v $
 *
 * $Author: kddevin $
 *
 * $Date: 1996/12/05 20:26:11 $
 *
 * $Revision: 5.0 $
 *
 * $Name:  $
 *====================================================================*/
#ifndef lint
static char *cvs_timersol_id =
  "$Id: md_timer_sol.c,v 5.0 1996/12/05 20:26:11 kddevin Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


#include <sys/time.h>
#include <sys/resource.h>

/* SUN SOLARIS timer */

double second()

{
  extern long int clock(void);

  return( (double) clock() / 1.0e+06);
}
