/*
 * Copyright 2013, Regents of Liu Shaopeng
 *
 * partsub.c
 * 
 * This file contains code for partitioning routines that from ParMetis
 *
 * Started 7/13/2013
 * Liu Shaopeng
 *
 * $Id: partsub.c,v 1.0 2013/07/13 21:47:20 karypis Exp $
 *
 */
#include <parmetisbin.h>
#include <string.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "c" {
#endif
//   For Fortran interface  //
//#define partsub PARTSUB
#define partsub0 partsub0_
extern int partsub0(idx_t *, idx_t *, idx_t *, idx_t *);
#ifdef __cplusplus
  }
#endif

#define NCON    3

/***********************************************************************************/
/*! This function tests the various graph partitioning and ordering routines */
/***********************************************************************************/
int partsub0(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *part)
{

  idx_t ncon, nparts, npes, mype, opt2, realcut;
  idx_t ier, icheck=-1;
  idx_t nvtxs;
  //graph_t graph;
  //idx_t *part;
  idx_t numflag=0, wgtflag=0, options[10], edgecut, ndims;
  real_t ipc2redist, *tpwgts = NULL, ubvec[MAXNCON];
  MPI_Status stat;

  MPI_Comm comm;

  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  if(mype == icheck) 
    printf("\n MYPE is %d, Begin call ParMetis subroutine....\n", mype);

  nvtxs = vtxdist[mype+1]-vtxdist[mype];

  if(mype == icheck) {
    printf("\n MYPE is %d, vtxdist == \n", mype);
    for(int ipe=0; ipe<=npes; ipe ++)
      printf(" %d ", vtxdist[ipe]);
    printf("\n xadj === \n");
    for(int i=0; i<=nvtxs;i++)
      printf( " %d ", xadj[i]);
    printf("\n adjncy === \n");
    for(int i=0; i<xadj[nvtxs];i++)
      printf( " %d ", adjncy[i]);
  }
  
  //if(mype == 0) 
    //printf("\n MYPE is %d, before Barrier....\n", mype);
  gkMPI_Barrier(comm);
  //if(mype == 0) 
    //printf("\n MYPE is %d, after Barrier....\n", mype);

  //part   = imalloc(nvtxs, "TestParMetis_V3: part");
  if(mype == icheck) 
    printf("\n MYPE is %d, nvtxs = %d \n", mype, nvtxs);
  tpwgts = rmalloc(npes, "ParMetis_V3: tpwgts");
  rset(1, 1.05, ubvec);

  /* Compute a good partition and move the graph. Do so quietly! */
  options[0] = 1;
  options[1] = 1;
  options[2] = 1;
  nparts = npes;
  wgtflag = 0;
  numflag = 0;
  ncon = 1; 
  rset(nparts*ncon, 1.0/(real_t)nparts, tpwgts);
  //if(mype == 0) 
    //printf("\n Begin call ParMETIS_V3_PartKway subroutine ... \n");
  ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, NULL, NULL, 
      &wgtflag, &numflag, &ncon, &npes, tpwgts, ubvec, options, &edgecut, 
      part, &comm);

  if (mype == icheck)
    printf("\nParMETIS_V3_RefineKway finished with ncon: %"PRIDX", nparts: %"PRIDX"\n", ncon, nparts);

/*
  options[3] = PARMETIS_PSR_UNCOUPLED;
  ParMETIS_V3_RefineKway(vtxdist, xadj, adjncy, NULL, 
      NULL, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, 
      &edgecut, part, &comm);
*/

/*  Reduce all part */
  /*
  if( mype == 0) {
    printf("knodeall = %d \n", vtxdist[npes]);
    idx_t *part_all, *part_dist;
    part_all = imalloc(vtxdist[npes], "Reduce all part");

    for(int ivtx = vtxdist[0]; ivtx < vtxdist[1]; ++ivtx)
      part_all[ivtx] = part[ivtx];

    for(int ipe=1; ipe<npes; ++ipe){
      int nvtxs_ipe = vtxdist[ipe+1]-vtxdist[ipe];
      part_dist = imalloc(nvtxs_ipe, "Reduce all part_dist");
      gkMPI_Recv((void *)part_dist, nvtxs_ipe, IDX_T, ipe, 0, comm, &stat);
      for(int ivtx = 0; ivtx < nvtxs_ipe; ivtx ++)
        part_all[ivtx+vtxdist[ipe]] = part_dist[ivtx];

      gk_free((void **) &part_dist, LTERM);
    }

    FILE *pfile;
    if(!(pfile=fopen("part.idx","wr"))){
      printf("Error! Open file failed !! \n");
      exit(-1);
    }
    for(int ivtx=0; ivtx<vtxdist[npes]; ivtx ++)
      fprintf(pfile," %d \n", part_all[ivtx]);
    fclose(pfile);

    gk_free((void **) &part_all, LTERM);
  }
  else{
    gkMPI_Send((void *)part, nvtxs, IDX_T, 0, 0, comm);
  }
  */
  //gk_free((void **)&tpwgts, &part, LTERM);

  return 0;
}
