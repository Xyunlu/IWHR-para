C---------------------------------------------------------------
C
       subroutine azsolv(N_update, N_extern, update,
     +                  bindx, val, b, x)
C
C---------------------------------------------------------------
C
       implicit none
       include "az_aztecf.h"
       include "mpif.h"
C
       integer N_update, N_extern
       integer bindx(0:*)
       double  precision val(0:*)
       integer update(0:*)
C           rhs and approximate solution
       double precision b(0:*), x(0:*)
       integer    i
C
C             See Aztec User's Guide for the variables that follow:
C
       integer proc_config(0:AZ_PROC_SIZE), options(0:AZ_OPTIONS_SIZE)
       double precision params(0:AZ_PARAMS_SIZE)
       integer data_org(0:1024)
       double precision status(0:AZ_STATUS_SIZE)

       integer, dimension(:), allocatable :: external
       integer, dimension(:), allocatable :: update_index, extern_index
       integer ierror

       integer myid, n0, n1, j, ichk
C
       ichk = 1
C
C      get number of processors and the name of this processor
C
       call AZ_set_proc_config(proc_config, MPI_COMM_WORLD)
C
       myid = proc_config(21)
C============================================================C
       if( myid .eq. ichk) then
c          print *,'bindx =', (bindx(i),i=0,bindx(N_update)-1)
c          print *,'val =', (val(i),i=0,bindx(N_update)-1)
c          print *,'diag_bindx = ', (bindx(i),i=0,N_update)
          print *,'diag:', (val(i),i=0,N_update-1)
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
         allocate(external(0:N_extern-1), extern_index(0:N_extern-1))
         allocate(update_index(0:N_update-1))
c         allocate(b(0:N_update-1), x(0:N_update+N_extern-1))

C
C      convert matrix to a local distributed matrix */
C
       call AZ_transform(proc_config,external,bindx,val,update,
     $                   update_index,extern_index,data_org,
     $                   N_update,0,0,0,0,AZ_MSR_MATRIX)
C============================================================C
       ichk = 1
       if( myid .eq. ichk ) then
c          print *,'external:', (external(i),i=0,N_extern-1) !global
c          print *,'extern_index:', (extern_index(i),i=0,N_extern-1) !local
c          print *,'update:',(update(i),i=0,N_update-1)          ! global
c          print *,'update_index:',(update_index(i),i=0,N_update-1)  ! local
c          print *,'non bindx:', (bindx(i),i=0, N_update)
c          do i=0, N_update-1
c            write(*,*) i,'==',update(i),'==', val(i)
c          enddo
c          print *,'col:'
c          do i=0, N_update-1
c             n0 = bindx(i)
c             n1 = bindx(i+1)-1
c             print *,(bindx(j),j=n0, n1)
c          enddo
c          print *,'data_org:', (data_org(i),i=0,20)
c          print *,'AZ_VBR_MATRIX=', AZ_VBR_MATRIX
c          print *,'AZ_MSR_MATRIX=', AZ_MSR_MATRIX
c          print *,'AZ_N_internal=', AZ_N_internal
c          print *,'AZ_N_border =', AZ_N_border
c          print *,'AZ_N_external =', AZ_N_external
c          print *,'AZ_N_int_blk =', AZ_N_int_blk
c          print *,'AZ_N_bord_blk =', AZ_N_bord_blk
c          print *,'AZ_N_ext_blk =', AZ_N_ext_blk
c          print *,'AZ_N_neigh =', AZ_N_neigh
c          print *,'AZ_total_send = ', AZ_total_send
c          print *,'AZ_name = ', AZ_name
c          print *,'AZ_neighbors = ', AZ_neighbors
c          print *,'AZ_rec_length = ', AZ_rec_length
c          print *,'AZ_send_length =', AZ_send_length
c          print *,'AZ_send_list =', AZ_send_list
c          print *,'data_org[AZ_send_list]=', data_org(AZ_send_list)
       endif
C============================================================C

C
C      initialize AZTEC options
C
       call AZ_defaults(options, params)
       options(AZ_max_iter) = 1000
       options(AZ_solver) = AZ_bicgstab
       options(AZ_precond) = AZ_sym_GS
       params(AZ_tol) = 1.0D-12
C
C      Set rhs (delta function at grid center) and initialize guess
C
       do 350 i = 0, N_update-1
          x(update_index(i)) = 0.0
350    continue
C
C      solve the system of equations using b  as the right hand side
C
       call AZ_solve(x, b, options, params, 0, bindx, 0, 0,
     $               0, val, data_org, status, proc_config)

c       if(myid .eq. 1) write(*,*) 'x:',(x(i),i=0,N_update-1)
C
c       call MPI_FINALIZE(ierror)
C
       stop
       end

C*********************************************************************
