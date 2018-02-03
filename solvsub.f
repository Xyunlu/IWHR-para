C---------------------------------------------------------------
C
       subroutine azsolv(N_update, N_extern, update,
     +                  bindx, val, rhs, x)
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
       double precision rhs(0:*), x(0:*)
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
       double precision, dimension(:), allocatable :: b, sol
       integer ierror

       integer myid, n0, n1, j, ichk
C
C      get number of processors and the name of this processor
C
       call AZ_set_proc_config(proc_config, MPI_COMM_WORLD)
C
       myid = proc_config(21)
C============================================================C
C============================================================C
         allocate(external(0:N_extern-1), extern_index(0:N_extern-1))
         allocate(update_index(0:N_update-1))

C
C      convert matrix to a local distributed matrix */
C
       call AZ_transform(proc_config,external,bindx,val,update,
     $                   update_index,extern_index,data_org,
     $                   N_update,0,0,0,0,AZ_MSR_MATRIX)
C============================================================C
C============================================================C

C
C      initialize AZTEC options
C
       call AZ_defaults(options, params)
       options(AZ_max_iter) = 100
       options(AZ_solver) = AZ_bicgstab
       options(AZ_precond) = AZ_sym_GS
       params(AZ_tol) = 1.0D-12
C
C      Set rhs (delta function at grid center) and initialize guess
C
       allocate( b(0: N_update-1), sol(0: N_update+N_extern-1) )
       do 350 i = 0, N_update-1
          sol(update_index(i)) = 0.0
          b(update_index(i)) = rhs(i)
350    continue
C
C      solve the system of equations using b  as the right hand side
C
       call AZ_solve(sol, b, options, params, 0, bindx, 0, 0,
     $               0, val, data_org, status, proc_config)

       do i=0, N_update-1
         x(i) = sol(update_index(i))
       enddo

       deallocate(external, extern_index, update_index, b, sol)
C
       return
       end

C*********************************************************************
