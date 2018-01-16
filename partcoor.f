c      implicit real*8(a-h,o-z)
C============================================================C
C     Input                                                  C
C         knode_t  : 节点总数                                C
C     Output                                                 C
C         nod_gidx : 对应节点整体编号                        C
C============================================================C
      subroutine partcoor(knodemax,knodei,nod_gidx,mnodeg,nnodeg,nodeg,
     +                    knode_t,coorg)
      USE ComData
      USE FEMData

      implicit none
C.... Input
      character*24 xyzfile, elmfile, iofile, outfile
      character*24 idfile, matfile, dispfile
      character*24 mshfile, resfile
c      character*6  ext
      logical      filflg
c.....
      integer      nod_gidx(1), knodei, knodemax, knode_t
      integer      mnodeg(2),nnodeg(2),nodeg(1)
      real         coorg(1)
C....
      integer, dimension(:), allocatable :: part
C....
c      integer, dimension(:), allocatable :: ipool
      integer, dimension(:), allocatable :: mnode_t, nnode_t
c      integer, dimension(:), allocatable :: nod_gid,nod_id,nod_gid_sort
c      integer, dimension(:), allocatable :: node, mnode, nnode
c      integer, dimension(:), allocatable :: xadj_neigh_s,adjncy_neigh_s 
c      integer, dimension(:), allocatable :: xadj_neigh_r,adjncy_neigh_r 
      integer, dimension(:), allocatable :: idx_elm_bound, num_elm_bound
      integer, dimension(:), allocatable :: nod_elm_bound
C.... FEM DATA
      character*1  material
c      double precision, dimension(:), allocatable :: coor

      integer kelem_t, kelem_max, knode_max, kelemg
      integer kelem_read,mmate_t,nmate_t
      integer i,j,k,ityp,numtyp,ne,num,num_t,inod,nod1,nod2
      integer ie, ipe, idx_pe, nodg, jnod, kcoor
      integer nnod, num_bound, num1, num2
      integer knode_1, knode_0,knode_bound
      integer is, knode_ipe,knode_neigh_s,knode_neigh_r,ir
      integer knode_send, knode_recv
      integer nodi, nodj, nne, kelem, kemate, imat
      integer is_belong, is_bound, num_elm_read
      integer nskp, numarg, ii, nbe, nleft, nright

      double precision emate(100000)
      integer          mmate(100), nmate(100)
c      integer          idx_neigh(100)
      integer          node_idx_local(100)
C.... Temp array
      integer,dimension(:),allocatable :: node_t
      integer          itmp, num_max, numg
c      double precision rtmp,rdata(100)
C.... Time
      double precision T0,T1,T_READ,T0P,T1P,T_PARTD
C....
      integer          AllocateStatus, DeAllocateStatus
C....

c      include 'mpif.h'
      integer iDBG, iCheck, mycomm, root

c      call MPI_Init(ierr)
      call MPI_COMM_DUP(MPI_COMM_WORLD, mycomm, ierr)

      call MPI_COMM_RANK(mycomm, mype, ierr)
      call MPI_COMM_SIZE(mycomm, npes, ierr)

      root=0

      T_READ = 0.D0
      T_PARTD = 0.D0
c      if(mype.eq.0 ) then
c        write(*,'(a)') 'Partcoor begin .........................'
c      endif
      material = 'y'
      numtyp = 3
      kcoor = 3
c      print *,'0000 coorg(1,*) =', (coorg(k),k=1,kcoor)

      allocate(mnode_t(numtyp))
      allocate(nnode_t(numtyp))
      allocate(mnode(numtyp))
      allocate(nnode(numtyp))
      allocate(idx_elm_bound(100000))
      allocate(num_elm_bound(numtyp+1))

C
      t0=MPI_WTime()
      allocate(part(knode_t))
      open(21,file='part.idx',form='formatted',status='old')
      do i=1, knode_t
        read(21,*) part(i)
      enddo
      close(21)

c      open(21,file='partb.idx',form='binary',status='old')
c      read(21) (part(i),i=1,knode_t)
c      close(21)
      t1=MPI_WTime()
      if(iDBG.eq.1 .and. mype.eq.0)
     +  write(*,'(i3,a,E12.5)') mype,'Read Part.idx: ', t1-t0
c
      write(ext,'(i6)') mype

      knode_max = knodemax
      kelem_t = 0
      num_t = 0
      do ityp=1, numtyp
        mnode_t(ityp) = mnodeg(ityp)
        nnode_t(ityp) = nnodeg(ityp)
        kelem_t = kelem_t + mnodeg(ityp)*nnodeg(ityp)
        num_t = num_t + mnode_t(ityp)
      enddo
      kelem_max = kelem_t/npes*15/10
      num_max = num_t/npes*15/10
      
      iDBG = 0
      if(iDBG.eq.1 .and. mype.eq.0) then
        print *,'knode_t, kelem_t = ', knode_t, kelem_t
        print *,'knode_max, kelem_max = ', knode_max, kelem_max
      endif

c      write(*,*) Mype, ', knodei=', knodei, ', knode_t=', knode_t

      allocate(nod_gid(knode_max))
      allocate(nod_id(knode_max))
      allocate(node(kelem_max), iE_lg(num_max))
      do i=1,num_max
        iE_lg(i) = 0
      enddo

C.... Array the nod id, find out the internal nodes
      do i=1,knodei
        nod_gid(i) = nod_gidx(i)
        nod_id(i) = 1
      enddo
      knode = knodei
      knode_1 = knode
C..
C..   Read element num_elm_read by num_elm_read
C..
c      if(mype.eq.0) then
c        write(*,'(a,i5)') 'Begin process element ..........', mype
c      endif

      T_READ = 0.D0
c      num_elm_read = 1000000
c      open(21,file=elmfile,form='binary',status='old')
      kelem = 0
      kelemg = 0
      numg = 0
      num_bound = 0
      num_elm_bound(1) = 0
c      print *,'numtyp =', numtyp
      do ityp=1, numtyp
        T0 = MPI_WTime()
c        read(21) num_t,nnod
        num_t = mnodeg(ityp)
        nnod = nnodeg(ityp)
        allocate(node_t(nnod))  ! temp array
c        read(21) num_t,nnod,((node_t((i-1)*nnod+j),j=1,nnod),i=1,num_t)
c        print *,'ityp =', ityp,'num_t, nnod=', num_t,nnod,',kelem_read=',kelem_read
        T1 = MPI_WTime()
        T_READ = T_READ + T1-T0
        num_elm_bound(ityp+1) = 0
c        num_t = mnode_t(ityp)
c        nnod = nnode_t(ityp)
        mnode(ityp) = 0
        nnode(ityp) = nnod
        num = 0
        nne = nnod
        if(material .eq. 'y' .or. material .eq. 'Y') nne = nne-1
        T0P = MPI_WTime()
        do ne=1,num_t
          is_belong = 0
          is_bound =0
          T0 = MPI_WTime()
          ! read(21) (node_t(j),j=1,nnod)
          do j=1, nnod
            node_t(j) = nodeg(kelemg+(ne-1)*nnod+j)
          enddo
          T1 = MPI_WTime()
          T_READ = T_READ + T1-T0
          do j=1,nne
            node_idx_local(j) = 0
            nodi = node_t(j)
C......     There are more than one node of this element belong to this block,
C           this element is a inner element of this block
            inod = 0
            if(part(nodi) .eq. mype) then
              do inod=1,knodei
                if(nod_gidx(inod) .eq. nodi) then
                  is_belong = 1
                  node_idx_local(j) = inod
                  exit
                endif
              enddo
              if( inod .gt. knodei) then
                write(*,'(a)') 'Error!Inner node not found in nod_gidx!'
                write(*,*) 'inod =', inod, ',knodei =', knodei
                write(*,*) 'nodi =', nodi,',part(nodi)=',part(nodi)
                write(*,*) 'mype = ', mype, knodei, nod_gidx(knodei)
                call endjob(ierr)
              endif
C......     There are some node of this element is not belong to this
C           block, this element is a bound element of this block
            else
              is_bound = 1
            endif
c            if(mype.eq.0) then
c              print *,'ne, nodi, part(nodi) =', ne, nodi, part(nodi)
c              print *,'is_belong, is_bound = ', is_belong, is_bound
c            endif
          enddo
C.....
          if( is_belong .eq. 0 ) goto 101
          num = num + 1
          numg = numg + 1
          if ( numg .gt. num_max) then
            write(*,*) 'Error! numg gt num_max!', numg, num_max
            stop 1
          endif
          if( is_bound .eq. 1) then
            num_bound = num_bound + 1
            num_elm_bound(ityp+1) = num_elm_bound(ityp+1) + 1
            idx_elm_bound(num_bound) = num
          endif
          do j=1,nne
            inod = node_idx_local(j)
            nodi = node_t(j)
C...       It is a inner node
            if(inod .gt. 0) then
              node(kelem+(num-1)*nnod+j) = inod
              if( is_bound .eq. 1) nod_id(inod) = 2
c              if( mype .eq. 1 .and. inod.eq.1) then
c                print *,'111, nodi, num, ne =', nodi, num, ne
c                print *,'is_bound =', is_bound
c                print *,'nod_id(inod) = ', nod_id(inod)
c              endif
              goto 111
              write(*,'(i5,a,i8)') mype, ' process: Error! -1- ', nodi
              write(*,*) 'ne = ', ne, 'part =', part(nodi)
              write(*,*) 'nod_gid =', (nod_gid(ii),ii=1,knode_1)
C...        It is a external node
            else
              do k = knode_1+1, knode
                if( nodi .eq. nod_gid(k) ) then
                  node(kelem+(num-1)*nnod+j) = k
                  goto 111
                endif
              enddo
              knode = knode + 1
              if(knode .gt. knodemax) then
                write(*,'(a,2i10)') 'knode more than knodemax!!',knode,knodemax
                stop 'Error in partcoor....'
              endif
              node(kelem+(num-1)*nnod+j) = knode
              nod_gid(knode) = nodi
              nod_id(knode) = 0
            endif
111         continue
          enddo
          node(kelem+num*nnod) = node_t(nnod) !材料编号
          iE_lg(numg) = ne ! 对应的全局单元编号
101       continue
        enddo   ! Enddo of num_t
        T1P = MPI_WTime()
        T_PARTD = T_PARTD + T1P-T0P
        mnode(ityp) = num
        kelem = kelem + num*nnod
        kelemg = kelemg + num_t*nnod
        deallocate(node_t)
c        read(21) mmate(ityp),nmate(ityp),
c     +           (rtmp,i=1,mmate(ityp)*nmate(ityp))
      enddo     ! Enndo of ityp, process element finished
c      close(21)
c      deallocate(nod_lid)
      if(iDBG.eq.1 .and. mype.eq.0) then
        write(*,'(a,E15.8,a)') 'Read elm file elasp time:', T_READ, ' s'
        write(*,'(a,E15.8,a)') 'Part elm elasp time:', T_PARTD, ' s'
      endif
c      write(*,*) mype,'==> numg =', numg, num_max
c      write(*,*) mype,'==> knode_1, knode=', knode_1, knode
c
c      write(*,*) mype,'===>',num, nnod
c      outfile='bdelm_' // trim(adjustl(ext))
c      open(21,file=outfile,form='binary',status='unknown')
c      write(21) num,nnod
c      write(21) (node(i),i=1,num*nnod)
c      close(21)

c      if(mype.eq.0) then

      num_elm_bound(1) = 0
      do ityp=1, numtyp
        num_elm_bound(ityp+1)=num_elm_bound(ityp)+num_elm_bound(ityp+1)
      enddo
c      print *,mype,'==>num_elm_bound:', (num_elm_bound(I),i=1,numtyp+1)

      iDBG = 0
      if( mype .eq. icheck .and. iDBG .eq. 1 ) then
        write(*,*) 'mype is ', mype, 'knode = ', knode
        write(*,*) 'nod_gid:', (nod_gid(i),i=1,knode)
        write(*,*) 'nod_id', (nod_id(i),i=1,knode)
        write(*,*) 'num_elm_bound =', (num_elm_bound(i),i=1,numtyp+1)
        write(*,*) 'idx_elm_bound =', (idx_elm_bound(i),i=1,num_bound)
      endif

      outfile='nodid_' // trim(adjustl(ext))
      open(21,file=outfile,form='formatted',status='unknown')
      write(21,*) (nod_id(i),i=1,knode)
      close(21)

      outfile='nodgid_' // trim(adjustl(ext))
      open(21,file=outfile,form='formatted',status='unknown')
      write(21,*) (nod_gid(i),i=1,knode)
      close(21)

c      if(mype.eq.0) then
c        write(*,'(a)') 'Begin process neighb info....'
c      endif
C.... Find neighbor proccessor
      allocate(nod_elm_bound(8*num_bound))

      kelem = 0
      knode_bound = 0
      num_ipe_belong = 0
      do ityp=1, numtyp
        num1 = num_elm_bound(ityp)+1
        num2 = num_elm_bound(ityp+1)
        num = mnode(ityp)
        nnod = nnode(ityp)
c        write(*,*) mype,',ityp,num,nnod,num1,num2==>',ityp,num,nnod,num1,num2
        nne = nnod
        if(material .eq. 'y' .or. material .eq. 'Y') nne = nne-1
        do ie=num1, num2
          ne = idx_elm_bound(ie)
          do inod=1, nne
            nodi = node(kelem+(ne-1)*nnod+inod)
            nodg = nod_gid(nodi)
            idx_pe = part(nodg)
            if( idx_pe .eq. mype ) goto 300   ! inner point
            do ipe = 1, num_ipe_belong
              if(idx_pe .eq. idx_neigh(ipe)) goto 200
            enddo
            num_ipe_belong = num_ipe_belong + 1
            idx_neigh(num_ipe_belong) = idx_pe
200         continue
            do jnod =1, knode_bound
              if(nodg .eq. nod_elm_bound(jnod)) goto 300
            enddo
            knode_bound = knode_bound + 1
            nod_elm_bound(knode_bound) = nodg
300         continue
          enddo
        enddo
        kelem = kelem + num*nnod
      enddo

c.... QSORT the neighbour proccessor number
      allocate(ipool(num_ipe_belong))
      do i=1,num_ipe_belong
        ipool(i)=i
      enddo
      nleft=0
      nright=num_ipe_belong-1
      call qsort(idx_neigh, ipool, nleft, nright)
      deallocate(ipool)
c      write(*,*) mype,' num_ipe_belong,idx_neigh:', num_ipe_belong,(idx_neigh(i),i=1,num_ipe_belong)

      iDBG = 0
      if( iDBG .eq. 1 ) then
        print *,'mype, num_ipe_belong =',mype, num_ipe_belong
        print *,'idx_neigh =', (idx_neigh(i),i=1,num_ipe_belong)
c        print *,'knode_bound =', knode_bound
c        print *,'nod_elm_bound=',(nod_elm_bound(i),i=1,knode_bound)
c        print *,'part of nod_bound:', (part(nod_elm_bound(i)),i=1,knode_bound)
      endif

      allocate(xadj_neigh_s(num_ipe_belong+1))
      allocate(adjncy_neigh_s(8*num_bound))
      allocate(xadj_neigh_r(num_ipe_belong+1))
      allocate(adjncy_neigh_r(8*num_bound))

      xadj_neigh_s(1) = 0
      xadj_neigh_r(1) = 0
      do i=1, num_ipe_belong
        knode_neigh_s = 0
        knode_neigh_r = 0
        ipe = idx_neigh(i)
        kelem = 0
        knode_ipe = 0
        do ityp=1, numtyp
          num1 = num_elm_bound(ityp)+1
          num2 = num_elm_bound(ityp+1)
          num = mnode(ityp)
          nnod = nnode(ityp)
          nne = nnod
          if(material .eq. 'y' .or. material .eq. 'Y') nne = nne-1
          do ie=num1, num2
            ne = idx_elm_bound(ie)
            do inod=1, nne
              nodi = node(kelem+(ne-1)*nnod+inod)
              nodg = nod_gid(nodi)
              if( part(nodg) .eq. ipe) goto 400
            enddo
            goto 410
400         continue
C....       The element have some nodes must be sent to ipe and recv
C...        some nodes data from ipe also
            do inod=1,nne
              nodi = node(kelem+(ne-1)*nnod+inod)
              nodg = nod_gid(nodi)
C....         this node must be sent to ipe
              if( part(nodg) .eq. mype) then
C....           check this node if exist in the list for send
                nskp = xadj_neigh_s(i)
                do is=1, knode_neigh_s
                  nodj = adjncy_neigh_s(nskp+is)
c                  if(nodj .eq. nodg) goto 420
                  if(nodj .eq. nodi) goto 420
                enddo
                knode_neigh_s = knode_neigh_s + 1
c                adjncy_neigh_s(nskp+knode_neigh_s) = nodg
                adjncy_neigh_s(nskp+knode_neigh_s) = nodi
C....         this nodes must be recv from ipe
              else if( part(nodg) .eq. ipe) then
                nskp = xadj_neigh_r(i)
                do ir=1, knode_neigh_r
                  nodj = adjncy_neigh_r(nskp+ir)
c                  if(nodj .eq. nodg) goto 420
                  if(nodj .eq. nodi) goto 420
                enddo
                knode_neigh_r = knode_neigh_r + 1
c                adjncy_neigh_r(nskp+knode_neigh_r) = nodg
                adjncy_neigh_r(nskp+knode_neigh_r) = nodi
              endif
420           continue
            enddo
410         continue
          enddo
          kelem = kelem+num*nnod
        enddo
        iDBG = 0
        icheck = 0
        if(iDBG .eq. 1 .and. mype .eq. icheck) then
          print *,'mype, i, ipe, knode_neigh_s, knode_neigh_r='
          print *,mype, i, ipe, knode_neigh_s, knode_neigh_r
        endif
        xadj_neigh_s(i+1) = xadj_neigh_s(i)+knode_neigh_s
        xadj_neigh_r(i+1) = xadj_neigh_r(i)+knode_neigh_r
      enddo

      iDBG=0
c      if(iDBG.eq.1 .and. mype.eq.icheck) then
      if(iDBG.eq.1) then
        print *,'ipe =', ipe
        print *,'idx_neigh =', (idx_neigh(i),i=1,num_ipe_belong)
        print *,'xadj_neigh_s =',(xadj_neigh_s(i),i=1,num_ipe_belong+1)
        knode_send = xadj_neigh_s(num_ipe_belong+1)
c        print *,'adjncy_neigh_s=',(adjncy_neigh_s(i),i=1,knode_send)
        print *,'xadj_neigh_r =',(xadj_neigh_r(i),i=1,num_ipe_belong+1)
        knode_recv = xadj_neigh_r(num_ipe_belong+1)
c        print *,'adjncy_neigh_r=',(adjncy_neigh_r(i),i=1,knode_recv)
      endif

C....

      allocate(nod_gid_sort(knode),ipool(knode))
      do i=1,knode
        nod_gid_sort(i)=i
        ipool(i)=nod_gid(i)
      enddo
      nleft=0
      nright=knode-1
      call qsort(ipool,nod_gid_sort,nleft,nright)
      deallocate(ipool)

c      if(mype.eq.0) then
c        write(*,'(i5,a)') mype,',Begin process coor file ...'
c      endif
      T0=MPI_WTime()
      allocate(coor(knode*kcoor), STAT=AllocateStatus)
      IF(AllocateStatus .NE. 0) STOP "* Not enough memory for coor *"
c      open(21,file=xyzfile,form='binary',status='old')
c      read(21) knode_t,kcoor
c      allocate(coor(knode_t*kcoor), STAT=AllocateStatus)
c      IF(AllocateStatus .NE. 0) STOP "* Not enough memory for coor *"
c      rewind(21)
c      print *,'111 coorg(1,*) =', (coorg(k),k=1,kcoor)
      knode_0=1
      inod=nod_gid_sort(knode_0)
      nodg=nod_gid(inod)
      do i=1,knode_t
c        read(21) (rdata(j),j=1,kcoor)
c        do j=1,kcoor
c          rdata(j) = coorg((i-1)*kcoor+j)
c        enddo
        if( i .eq. nodg) then
          do k=1,kcoor
            coor((inod-1)*kcoor+k) = coorg((i-1)*kcoor+k)
          enddo
          knode_0=knode_0+1
          if(knode_0 .gt. knode) exit
          inod=nod_gid_sort(knode_0)
          nodg=nod_gid(inod)
        endif
      enddo
c      close(21)
      deallocate(nod_gid_sort)
      T1=MPI_WTime()
c      write(*,'(i4,a,E12.5)') mype,',Read and Part coor0:', T1-T0

      xyzfile = 'coor0_' // trim(adjustl(ext))
      open(21,file=xyzfile,form='binary',status='unknown')
      write(21) knode,kcoor,
     +          ((coor((i-1)*kcoor+j),j=1,kcoor),i=1,knode)
      close(21)

C....
C....   OUTPUT the part info
C....
c      outfile = trim(prjname) // adjustl(ext)
c      print *,'outfile=', trim(outfile)
c      open(21,file=outfile,form='formatted',status='unknown')
C.... Coordinate
c      write(21,'(i10, i5)') knode, kcoor
c      do i=1, knode
c        nodg = nod_gid(i)
c        write(21,'(i8,i10,3E15.6,i3)') i,nodg,(coor((nodg-1)*kcoor+j),j=1,kcoor),nod_id(i)
c      enddo
c      write(21,*)
C.... Element and Material
c      write(21,'(2i15)') kelem, kemate
c      write(21,'(3i5)') numtyp, nbe, numtyp-nbe
c      kelem = 0
c      kemate = 0
c      do ityp=1, numtyp
c        write(21,'(i10, i5)') mnode(ityp), nnode(ityp)
c        num = mnode(ityp)
c        nnod = nnode(ityp)
c        print *,'ityp, num, nnod = ', ityp, num, nnod
c        do ie=1, num
c          write(21,'(10i10)') ie,(node(kelem+(ie-1)*nnod+j),j=1,nnod)
c        enddo
c        write(21,'(2i5)') mmate(ityp),nmate(ityp)
c        do imat=1, mmate(ityp)
c          write(21,'(i5,100E14.5)') imat,(emate(kemate+(imat-1)*nmate(ityp)+j),j=1,nmate(ityp))
c        enddo
c        kelem = kelem + num*nnod
c        kemate = kemate + mmate(ityp)*nmate(ityp)
c      enddo
c      write(21,*) 
C.... neighor
c      if(mype.eq.0) then
c        write(*,'(a)') 'Begin output neighor info ....'
c      endif
c      write(21,'(i15)') num_ipe_belong
c      do i=1,num_ipe_belong
c        nod1 = xadj_neigh_s(i)+1
c        nod2 = xadj_neigh_s(i+1)
c        write(21,'(2i10)') idx_neigh(i), nod2-nod1+1
c        write(21,'(10i10)') (adjncy_neigh_s(j),j=nod1,nod2)
c      enddo
c      write(21,'(i15)') num_ipe_belong
c      do i=1,num_ipe_belong
c        nod1 = xadj_neigh_r(i)+1
c        nod2 = xadj_neigh_r(i+1)
c        write(21,'(2i10)') idx_neigh(i), nod2-nod1+1
c        write(21,'(10i10)') (adjncy_neigh_r(j),j=nod1,nod2)
c      enddo
c      close(21)
c
      outfile = 'conn_' // adjustl(ext)
      open(21,file=outfile,form='formatted',status='unknown')
      write(21,'(i15)') num_ipe_belong
      do i=1,num_ipe_belong
        nod1 = xadj_neigh_s(i)+1
        nod2 = xadj_neigh_s(i+1)
        write(21,'(2i10)') idx_neigh(i), nod2-nod1+1
        write(21,'(10i10)') (adjncy_neigh_s(j),j=nod1,nod2)
        nod1 = xadj_neigh_r(i)+1
        nod2 = xadj_neigh_r(i+1)
        write(21,'(2i10)') idx_neigh(i), nod2-nod1+1
        write(21,'(10i10)') (adjncy_neigh_r(j),j=nod1,nod2)
      enddo
      close(21)

C..... Output for Gid postprocess
      iDBG = 1
      if(iDBG .eq. 1) then
        mshfile = 'part' // trim(adjustl(ext)) // '.flavia.msh'
        resfile = 'part' // trim(adjustl(ext)) // '.flavia.res'
c        print *,'mshfile is ', mshfile, 'resfile is ', resfile

        open(21,file=mshfile,form='formatted',status='unknown')
        write(21,'(a,a)') 'Mesh "partc8" Dimension 3 Elemtype ',
     +                    'Hexahedra Nnode 8'
        write(21,'(a)') 'Coordinates'
        do i=1, knode
          write(21,12) i,(coor((i-1)*kcoor+j),j=1,kcoor)
        enddo
        write(21,'(a)') 'End coordinates'
        kelem = 0
        write(21,'(a)') 'Elements'
        num = mnode(1)
        nnod = nnode(1)
        do i=1, num
          write(21,13) i,
     +                         node(kelem+(i-1)*nnod+1),
     +                         node(kelem+(i-1)*nnod+2),
     +                         node(kelem+(i-1)*nnod+4),
     +                         node(kelem+(i-1)*nnod+3),
     +                         node(kelem+(i-1)*nnod+5),
     +                         node(kelem+(i-1)*nnod+6),
     +                         node(kelem+(i-1)*nnod+8),
     +                         node(kelem+(i-1)*nnod+7),
     +                         node(kelem+(i-1)*nnod+9)
        enddo
        write(21,'(a)') 'End elements'
        kelem = kelem+num*nnod
        if( numtyp .eq. 2) then
          write(21,'(a,a)') 'Mesh "surfq4" Dimension 3 Elemtype ',
     +                      'Quadrilateral Nnode 4'
          write(21,'(a)') 'Coordinates'
          write(21,'(a)') 'End coordinates'
          write(21,'(a)') 'Elements'
          num = mnode(2)
          nnod = nnode(2)
          do i=1, num
            write(21,13) i,
     +                   (node(kelem+(i-1)*nnod+j),j=1,nnod)
          enddo
          write(21,'(a)') 'End elements'
        endif
        close(21)

        open(21,file=resfile,form='formatted',status='unknown')
        write(21,'(a)') 'GID Post Results File 1.0'
        write(21,'(a)')
        write(21,'(a,a)') 'Result "nod_id" "Load Analysis" 1',
     +                    ' Scalar OnNodes'
        write(21,'(a)') 'ComponentNames "id"'
        write(21,'(a)') 'Values'
        do i=1,knode
          write(21,12) i, dble(nod_id(i))
        enddo
        write(21,*) 'End values'
        close(21)
      endif

12    format (1X,I8,10e17.8)
13    format (1X,15I8)

c      deallocate(coor)

9999  continue

      deallocate(part)
      deallocate(nod_gid)
      deallocate(nod_id)
c      deallocate(xadj_neigh_s)
c      deallocate(adjncy_neigh_s)
c      deallocate(xadj_neigh_r)
c      deallocate(adjncy_neigh_r)
      deallocate(nod_elm_bound)

      deallocate(mnode_t)
      deallocate(nnode_t)
c      deallocate(mnode)
c      deallocate(nnode)
      deallocate(num_elm_bound)
      deallocate(idx_elm_bound)
c      deallocate(node)

1111  continue

c      write(*,'(i5,a,i5)') mype,' Partcoor finished!!', mype

      call MPI_COMM_Free(mycomm, ierr)

      return
      end
