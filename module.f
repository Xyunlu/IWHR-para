C===========================================================C
      Module ComData
        integer MAXN, MAXN1
        parameter( MAXN=990000,MAXN1=510000 )

        include 'mpif.h'
        integer :: ierr, mype, npes, jpe
        integer, allocatable :: idata(:)
        real*8, allocatable :: rdata(:)

        ! 节点自由度列表，JRL-局部，JRG-整体方程号
        integer, allocatable :: JRL(:,:), JRG(:,:)
        real*8, allocatable :: vol(:), etl(:), utl(:), gamtl(:)

        ! 节点分区信息
        ! nod_gid: 节点整体编号
        ! nod_id: 节点分类信息，>0-内部节点；0-边界节点(external)
        integer, dimension(:), allocatable :: nod_gid,nod_id,nod_gid_sort

        ! 子区域之间通信信息标识
        integer :: num_ipe_belong       ! 本区域的所有边界节点属于num_ipe_belong个子区域
        integer, allocatable :: idx_neigh(:)
        integer, allocatable :: xadj_neigh_s(:), adjncy_neigh_s(:)
        integer, allocatable :: xadj_neigh_r(:), adjncy_neigh_r(:)

        character*24 filename,ext

      End Module ComData

      Module FEMData
        ! 子区域单元和节点的整体编号, iE_lg-单元, iN_lg-节点
        integer, dimension(:), allocatable :: iE_lg, iN_lg

        ! integer,dimension(:,:),allocatable :: IPE_Sub
        integer :: knode, knode_i, neq, neq_tol
        real, dimension(:),allocatable :: coor

        integer,dimension(:),allocatable :: mnode,nnode,node

        integer,dimension(:),allocatable :: neighs, neighr
        integer,dimension(:),allocatable :: ipool
        double precision, dimension(:), allocatable :: rpool

        integer :: kdgof

        ! 定义求解器所需数据
        integer :: nnz
        integer, dimension(:), allocatable :: na, ia
        double precision, dimension(:), allocatable :: am, rhs, sol

c        integer, dimension(:), allocatable :: maplg, imaplg
        integer :: N_update, N_external
        integer, dimension(:), allocatable :: update, external
        integer, dimension(:), allocatable :: update_index, extern_index

        integer, dimension(:), allocatable :: bindx
        real*8, dimension(:), allocatable :: val
        integer :: nupdate(0:1024)

      End Module FEMData
C===========================================================C
       module solvmodule
       contains
         subroutine DCRS2DMSR(myid,N,nextern,na,ia,a,
     +                      update,bindx,val,b,x)
         implicit none
         integer :: N,update(*),na(*),ia(*)
         double precision :: a(*)
         integer, dimension(:), allocatable :: bindx, external
         double precision, dimension(:), allocatable :: val,b,x
         integer myid,i,j,k,l,n0,n1,row,col,nnz,nextern,ise,ise1
         integer ierr, nodiag, N_extern

         nnz = 0
         do i=1, N
           n0 = na(i)
           n1 = na(i+1)
           nnz = nnz + n1 - n0
           update(i) = update(i) - 1
         enddo

c         if(myid.eq.0) then
c           print *,'In crs2dmsr: myid,N_update,N_external=', myid,N,nextern
c           print *,'update:', (update(i),i=1,N)
c           print *,'na(1), na(N+1),ia(nnz) =', na(1), na(N+1),ia(nnz)
c           print *,'a(1), a(N+1),a(nnz) =', a(1), a(N+1),a(nnz)
c         endif

         if( allocated(bindx) ) deallocate(bindx, val, external) 
         allocate(bindx(nnz+1),val(nnz+1),external(nextern))

         N_extern = 0
         do i=1, nnz
           col = ia(i)
           ise = 1
           do j=1, N
             if( col .eq. update(j)+1 ) then
               ise = 0
               exit
             endif
           enddo
           if ( ise == 0 ) cycle
           ise = 1
           do j=1, N_extern
             if( col == external(j)) then
               ise = 0
               cycle
             endif
           enddo
           if( ise == 1) then
             N_extern = N_extern + 1
             external(N_extern) = col
           endif
         enddo
c         print *,'In crs2DMSR:myid, N_extern = ', myid, N_extern
c         print *,'external:', (external(j),j=1, N_extern)

         k = N+1
         bindx(1) = k
         N_extern = 0
         do i=1, N
           row = update(i) + 1
           n0 = na(i)+1
           n1 = na(i+1)
           nodiag = 0
           do j=n0, n1
             col = ia(j)
             if( col .eq. row ) then
               val(i) = a(j)
               nodiag = 1
             else
               k = k+1
               if( k > nnz+1) then
                 print *,'Error!! k > nnz+1 !!!', k, nnz+1
                 call My_endjob(ierr)
               endif
               val(k) = a(j)
               bindx(k) = col-1
               ise = 1
               do l=1, N
                 if( col-1 .eq. update(l)) then
                   ise = 0
                   exit
                 endif
               enddo
               if(ise == 1)  then
                 ise1 = 1
                 do l=1, N_extern
                   if( col-1 == external(l)) then
                     ise1 = 0
                     exit
                   endif
                 enddo
                 if( ise1 == 1 ) then
                   N_extern = N_extern + 1
                   external(N_extern) = col -1
c                   if(myid .eq. 0) then
c                     print *,'N_extern, external=', N_extern, col-1
c                   endif
                 endif
               endif
             endif
           enddo
           if( nodiag == 0) then
             print *,'Error! Not found diag of a matrix!'
             print *,'myid,i,row=',myid,i,row
             call My_endjob(ierr)
           endif
           bindx(i+1) = k
         enddo

         if( N_extern > nextern) then
           write(*,*) 'Error in DCRS2DMRS!'
           write(*,*) 'myid, N_extern, nextern=',myid,N_extern,nextern
           write(*,*) 'external:', (external(i),i=1, nextern)
           call My_endjob(ierr)
         endif
c         print *,'myid, N_update, N_extern=', myid, N, N_extern
         if( myid.eq.0) then
c           print *,'nextern = ', nextern
c           print *,'bindx:', (bindx(i),i=1,nnz+1)
c           print *,'val:', (val(i),i=1,nnz+1)
         endif

         nextern = N_extern
         if ( allocated(b) ) deallocate(b, x)
         allocate(b(N), x(N+nextern))

         end subroutine
       end module solvmodule
