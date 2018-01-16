C---------------------------------------------------------------
C
       subroutine azsolv(N_update, update, update_index,
     +                  external, extern_index,
     +                  bindx, val)
C
C---------------------------------------------------------------
C
       implicit none
       include "az_aztecf.h"
       include "mpif.h"
C
C           POISSON EQUATION WILL BE SOLVED ON an n x n GRID.
C	    NOTE: n should be odd for rhs to be properly set.
C
       double precision b(0:1024),x(0:1024)
C           rhs and approximate solution
       integer    i
C
C             See Aztec User's Guide for the variables that follow:
C
       integer proc_config(0:AZ_PROC_SIZE), options(0:AZ_OPTIONS_SIZE)
       double precision params(0:AZ_PARAMS_SIZE)
       integer data_org(0:1024)
       double precision status(0:AZ_STATUS_SIZE)

       integer update(0:*), external(0:*)
       integer update_index(0:*), extern_index(0:*)
       integer bindx(0:*)
       double  precision val(0:*)
       integer N_update,ierror

       integer myid, n0, n1, j, ichk
C
       ichk = 1
C
c       call MPI_INIT(ierror)
C
C           # of unknowns updated on this node
c       n = 3
C
C      get number of processors and the name of this processor
C
       call AZ_set_proc_config(proc_config, MPI_COMM_WORLD)
C
C      Define paritioning:matrix rows in ascending order assigned
C      to this node
C
c       nrow = n*n

       myid = proc_config(21)
c       print *,"myproc:", myid
c       print *,"24:", proc_config(24)
c       call AZ_read_update(N_update,update,proc_config,nrow,1,0)
C============================================================C
       if( myid .eq. ichk) then
         print *,"myproc:", proc_config(21),'N_update=', N_update
         print *,'update:', (update(i), i=0, N_update-1)
       endif
C============================================================C
C
C      create the matrix: each processor creates only rows
C      appearing in update[] (using global col. numbers).
C
c       bindx(0) = N_update+1
c       do 250 i = 0, N_update-1
c          call create_matrix_row_5pt(update(i),i,val,bindx)
c250    continue
C============================================================C
       if( myid .eq. ichk) then
          print *,'bindx =', (bindx(i),i=0,bindx(N_update)-1)
          print *,'val =', (val(i),i=0,bindx(N_update)-1)
c          print *,'diag_bindx = ', (bindx(i),i=0,N_update)
c          print *,'diag:', (val(i),i=0,N_update-1)
c          print *,'maxa =', bindx(N_update)-bindx(0)
c          print *,'col:'
c          do i=0, N_update-1
c            print *,'-----------', i, '---------', update(i)
c            n0 = bindx(i)
c            n1 = bindx(i+1)-1
c            print *,(bindx(j),val(j),j=n0, n1)
c          enddo
       endif
C============================================================C
C
C      convert matrix to a local distributed matrix */
C
       call AZ_transform(proc_config,external,bindx,val,update,
     $                   update_index,extern_index,data_org,
     $                   N_update,0,0,0,0,AZ_MSR_MATRIX)
C============================================================C
       if( myid .eq. ichk ) then
          print *,'external:', (external(i),i=0,9) !global
          print *,'extern_index:', (extern_index(i),i=0,9) !local
          print *,'update:',(update(i),i=0,N_update-1)          ! global
          print *,'update_index:',(update_index(i),i=0,N_update-1)  ! local
          print *,'non bindx:', (bindx(i),i=0, N_update)
          print *,'diag:', (val(i),i=0, N_update-1)
          print *,'col:'
          do i=0, N_update-1
             n0 = bindx(i)
             n1 = bindx(i+1)-1
             print *,(bindx(j),j=n0, n1)
          enddo
          print *,'data_org:', (data_org(i),i=0,20)
          print *,'AZ_VBR_MATRIX=', AZ_VBR_MATRIX
          print *,'AZ_MSR_MATRIX=', AZ_MSR_MATRIX
          print *,'AZ_N_internal=', AZ_N_internal
          print *,'AZ_N_border =', AZ_N_border
          print *,'AZ_N_external =', AZ_N_external
          print *,'AZ_N_int_blk =', AZ_N_int_blk
          print *,'AZ_N_bord_blk =', AZ_N_bord_blk
          print *,'AZ_N_ext_blk =', AZ_N_ext_blk
          print *,'AZ_N_neigh =', AZ_N_neigh
          print *,'AZ_total_send = ', AZ_total_send
          print *,'AZ_name = ', AZ_name
          print *,'AZ_neighbors = ', AZ_neighbors
          print *,'AZ_rec_length = ', AZ_rec_length
          print *,'AZ_send_length =', AZ_send_length
          print *,'AZ_send_list =', AZ_send_list
          print *,'data_org[AZ_send_list]=', data_org(AZ_send_list)
       endif
C============================================================C

C
C      initialize AZTEC options
C
       call AZ_defaults(options, params)
C
C      Set rhs (delta function at grid center) and initialize guess
C
       do 350 i = 0, N_update-1
          x(update_index(i)) = 0.0
          b(update_index(i)) = 0.0
          if (update(i) .eq. 0) b(update_index(i)) = 1.0
350    continue
C
C      solve the system of equations using b  as the right hand side
C
       call AZ_solve(x,b, options, params, 0,bindx,0,0,
     $               0,val, data_org, status, proc_config)

C
c       call MPI_FINALIZE(ierror)
C
       stop
       end

C*********************************************************************
C*********************************************************************
C
       subroutine create_matrix_row_5pt(row,location,val,bindx)
C
       integer row,location,bindx(0:*)
       double precision val(0:*)
       integer   n
       common    /global/ n
C
C Add one row to an MSR matrix corresponding to a 5pt discrete
C approximation to the 2D Poisson operator on an n x n square.
C
C  Parameters:
C     row          == global row number of the new row to be added.
C     location     == local row where diagonal of the new row will be stored.
C     val,bindx    == (see user's guide). On output, val[] and bindx[]
C                     are appended such that the new row has been added.
C
       integer k
C
C      check neighbors in each direction and add nonzero if neighbor exits
C
       k = bindx(location)
       bindx(k)  = row + 1
       if (mod(row,n) .ne. n-1) then
          val(k) = -1.
          k = k + 1
       endif
       bindx(k)  = row - 1
       if (mod(row,n) .ne.   0) then
          val(k) = -1.
          k = k + 1
       endif
       bindx(k)  = row + n
       if (mod(row/n,n) .ne. n-1) then
          val(k) = -1.
          k = k + 1
       endif
       bindx(k)  = row - n
       if (mod(row/n,n) .ne.   0) then
          val(k) = -1.
          k = k + 1
       endif

       bindx(location+1) = k
       val(location)   = 4.
       return
       end
