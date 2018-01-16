       module mymodule
       contains
         subroutine CRS2DMSR(myid,N,nextern,update,na,ia,a,
     +                      bindx,val,b,x)
         implicit none
         integer :: N,update(*),na(*),ia(*)
         double precision :: a(*)
         integer, dimension(:), allocatable :: bindx
         double precision, dimension(:), allocatable :: val,b,x
         integer myid,i,j,k,l,n0,n1,row,col,nnz,nextern,ise

         nnz = 0
         do i=1, N
           row = update(i)
           n0 = na(row)
           n1 = na(row+1)
           nnz = nnz + n1 - n0
           update(i) = update(i) - 1
         enddo

         allocate(bindx(nnz+1),val(nnz+1))

         bindx(1) = N+1
         k = N+1
         nextern = 0
         do i=1, N
           row = update(i)+1
           n0 = na(row)
           n1 = na(row+1)-1
           do j=n0, n1
             col = ia(j)
             if( col .eq. row) then
               val(i) = a(j)
             else
               k = k+1
               val(k) = a(j)
               bindx(k) = col-1
               ise = 1
               do l=1, N
                 if( col-1 .eq. update(l)) then
                   ise = 0
                   exit
                 endif
               enddo
               if(ise .ne. 0) nextern = nextern + 1
             endif
           enddo
           bindx(i+1) = k
         enddo

         if( myid.eq.0) then
           print *,'nextern = ', nextern
c           print *,'bindx:', (bindx(i),i=1,nnz+1)
c           print *,'val:', (val(i),i=1,nnz+1)
         endif
c         allocate(external(nextern), extern_index(nextern))
         allocate(b(N), x(N+nextern))

         end subroutine
       end module mymodule
