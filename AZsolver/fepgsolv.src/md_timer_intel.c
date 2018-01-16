/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: md_timer_intel.c,v $
 *
 * $Author: tuminaro $
 *
 * $Date: 1999/09/30 17:04:27 $
 *
 * $Revision: 5.1 $
 *
 * $Name:  $
 *====================================================================*/
#ifndef lint
static char *cvs_timerint_id =
  "$Id: md_timer_intel.c,v 5.1 1999/09/30 17:04:27 tuminaro Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


/* INTEL timer */
extern double dclock(void);
extern double second(void);

double second(void)

{


  return (dclock());

} /* second */
