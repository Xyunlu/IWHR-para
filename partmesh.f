C============================================================C
C     Descript                                               C
C         将网格按对偶图的方式进行初步分区                   C
C     Input                                                  C
C         knode_t: 节点总数                                  C
C     Output                                                 C
C         nodg   : 局部节点对应的整体节点编号                C
C============================================================C
c     call partmesh(maxknode_i,nodg,mnode,nnode,node_glb,knode_i,knode_t)

      subroutine partmesh(maxknode_i, nodg,
c     +                    maxne, numg, nodeg, me,
     +                    mnodeg, nnodeg, node_g,
     +                    knode_i, knode_t)
      implicit none
c     Input data
c      character*24 prjname, field
      integer      mnodeg(2), nnodeg(2)
      integer      node_g(1)
c      integer      maxne, numg, nodeg(maxne, 8), me(maxne)
      integer      nodg(1), maxknode_i, knode_i

      include 'mpif.h'

      character*24 xyzfile,elmfile,iofile,grphfile
      character*8  ext

      integer, dimension(:), allocatable :: NODE, LM
      integer*8, dimension(:), allocatable :: NAJ, NAP

      integer*8, dimension(:), allocatable :: NUMCOL, NA
      integer*8, dimension(:), allocatable :: vtxdist
      integer*8, dimension(:), allocatable :: ipart

      integer*8, dimension(:), allocatable :: ipart_all, ipart_dist

      integer            nvtxs

      integer, dimension(:), allocatable :: ipool, nodnum

      integer, dimension(:), allocatable :: noddist, noddist0

      integer            MAXNODE, MAXT, NEDGS, JNA, ISIZE, MAXA
      integer            knode, MAXBAND
      integer            N1, N2

      integer            node_tmp(100)
      integer            numarg
c
      integer            ip,numtyp,nod_ist,nod_iend,ipe
      integer            nproc,knode_t,knode_dist,knode_dist0
      integer            num,nnode,nne,nod_g,nodipe,nodipe0,nodipe1
      integer            mmate,nmate
      integer            i,j,inod,icl,ne,l,mynod,nodi
      integer            ityp,kelem,kelemg,kemate,ie,num_ef,num_ef_max
c
      integer            ierr, mycomm, istatr(MPI_STATUS_SIZE), irequr
      integer            itag, istats(MPI_STATUS_SIZE), irequs, iflgs
      integer            iDBG, ichk, mype, npes, root
      integer            AllocateStatus, DeAllocateStatus
c
      character*1        material
      double precision   rtemp
      integer            itemp
C...  For test
      integer            i0
C.... Time
      double precision   T0,T1,T_READ,T_M2G,T_PART,T_OUT,T_Comm
      double precision   T10,T11,T20,T21,T110,T111

c      call MPI_Init(ierr)

      call MPI_COMM_DUP(MPI_COMM_WORLD, mycomm, ierr)
      call MPI_COMM_RANK(mycomm, mype, ierr)
      call MPI_COMM_SIZE(mycomm, npes, ierr)

      root=0
      t0=MPI_WTime()

C.............  TEST
c      if( mype .eq. 0) then
c        i0=1
c        call MPI_Send(i0,1,MPI_integer,1,10000,mycomm,ierr)
c      elseif( mype .eq. 1) then
c        call MPI_Recv(i0,1,MPI_integer,0,10000,mycomm,istat,ierr)
c        write(*,*) 'TEST;;; i0=', i0
c      endif
c      call MPI_COMM_RANK(mycomm, mype, ierr)
c      call MPI_COMM_SIZE(mycomm, npes, ierr)
c      print *,'mype=',mype
c
      T_READ=0.D0
      T_M2G=0.D0
      T_PART=0.D0
      T_Out=0.D0
      T_Comm=0.D0
c
      ichk = 0 
      iDBG = 0
c
c      numarg = 0
c      numarg = numarg + 1
c      call getarg(numarg, prjname)
c      numarg = numarg + 1
c      call getarg(numarg, field)
c      iofile = trim(prjname) // trim(field) // '.io'
c      xyzfile = trim(prjname) // '.xyz'
c      elmfile = trim(prjname) // '.elm'
c      xyzfile = 'coor0'
c      elmfile = 'elem0'
c      print *,'mype=',mype,',iofile=',iofile
c
C ..... OPEN IO FILE
c      if( mype .eq. 0) then
c        OPEN(21,FILE=iofile,FORM='FORMATTED',STATUS='OLD')
c        READ(21,*) material
c        READ(21,*) NUMTYP
c        CLOSE(21)
c        material = 'y'
c        NUMTYP = 2
c      endif
c      call MPI_Bcast(material, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
c      call MPI_Bcast(NUMTYP, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      material = 'y'
      NUMTYP = 3
 
C.... OPEN COOR File
c      if( mype .eq. 0) then
c        open(21,file=xyzfile,form='binary',status='old')
c        read(21) knode_t
c        close(21)
c      endif
      if(iDBG.eq.1 .and. mype .eq. ichk) 
     +  write(*,'(a,i12,a,i3,a,a)') 'PartMesh: knode_t = ', knode_t,
     +                      ',NUMTYP=',numtyp,
     +                      ',material =', material
c      print *,'mype, knode_t=',mype, knode_t

C.... OPEN NPART
c      open(21,file='nproc',form='formatted',status='old')
c      read(21,*) nproc
c      close(21)
c      if(nproc .ne. npes) then
c        write(*,'(a)') 'nproc not equal npes, set to npes!'
      nproc=npes
      if( mype .eq. 0) then
        open(21,file='nproc',form='formatted',status='unknown')
        write(21,*) nproc
        close(21)
      endif
c      endif
      T1=MPI_WTime()
c      write(*,'(i3,a,E12.5)') mype,',Prepare time :', T1-T0
c
      t0=MPI_WTime()
      allocate(vtxdist(npes+1), STAT = AllocateStatus)
      IF(AllocateStatus .NE. 0) STOP "* Not enough memory for vtxdist *"
c
      vtxdist(1) = 0
      knode_dist = knode_t / npes
      do ip=1, npes
        if(ip .gt. npes-mod(knode_t,npes) ) then
          vtxdist(ip+1) = knode_dist + 1
        else
          vtxdist(ip+1) = knode_dist
        endif
      enddo
c      write(*,'(100i10)') mype,npes,(vtxdist(i),i=1,npes+1)
c      
      do ip = 2, npes+1
        vtxdist(ip) = vtxdist(ip-1) + vtxdist(ip)
      enddo
c      write(*,'(100i10)') npes,(vtxdist(i),i=1,npes+1)
      if(vtxdist(npes+1) .ne. knode_t) then
        write(*,'(a,2i12)') 'VTXDIST Error!', vtxdist(npes+1), knode_t
        goto 9999
      endif
c      if(mype .eq. ichk) then
c        write(*,*) 'vtxdist:',(vtxdist(ip),ip=1,nproc+1)
c      endif
c
      nod_ist = vtxdist(mype+1)+1
      nod_iend = vtxdist(mype+2)
c      print *,'mype=',mype,'nod_ist,nod_iend=',nod_ist,nod_iend
      knode = nod_iend - nod_ist + 1
      MAXBAND = 50
      MAXNODE = knode+1
      MAXT = (knode)*MAXBAND
c      write(*,*) 'MAXNODE = ', MAXNODE
c      write(*,*) 'MAXT = ', MAXT
c
      allocate(LM(1000), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "*** Not enough memory ***"
      allocate(NUMCOL(MAXNODE), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "*** Not enough memory ***"
      allocate(NAP(MAXT), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "*** Not enough memory ***"
      allocate(NAJ(MAXT), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "*** Not enough memory ***"
      T1=MPI_WTime()
c      write(*,'(i3,a,E12.5)') mype,',Allocate time: ', t1-t0
c 
C.......OPEN ELEM0 FILE
c      if( mype .eq. 0) OPEN (23,FILE=elmfile,FORM='binary',STATUS='OLD')
c      OPEN (23,FILE=elmfile,FORM='binary',STATUS='OLD')
C.... Allocate ipool for element information
      kelem = 0
      num_ef_max = 0
      do ityp=1, numtyp
        num_ef = 0
c        read(23) num,nnode
        num = mnodeg(ityp)
        nnode = nnodeg(ityp)
        nne = nnode
        if( material .eq. 'y' .or. material .eq. 'Y') nne = nne-1
        isize=nnode
        allocate(ipool(isize), STAT = AllocateStatus)
        do ne=1,num
          t0=MPI_WTime()
c          if( mype .eq. 0) read(23) (ipool(j),j=1,nnode)
c          call MPI_Bcast(ipool(1), nnode, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
c          read(23) (ipool(j),j=1,nnode)
          do j=1, nne
c            ipool(j) = nodeg(ne,j)
            ipool(j) = node_g(kelem+(ne-1)*nnode+j)
          enddo
c          ipool(nnode) = me(ne)
          ipool(nnode) = node_g(kelem+ne*nnode)
c          if( ne .eq. 1) then
c            write(*,*) 'node_g:', (node_g(j), j=1,nnode)
c            write(*,*) kelem,',ipool=', (ipool(j),j=1,nnode)
c          endif
          t1=MPI_WTime()
          T_READ = T_READ+t1-t0
          do j=1,nne
            nodi = ipool(j)
            if(nodi .ge. nod_ist .and. nodi .le. nod_iend) goto 101
          enddo
          goto 111
101       num_ef = num_ef + 1
111       continue
          if(num_ef .gt. num_ef_max) num_ef_max = num_ef
        enddo
c        if( mype .eq. 0) read(23) mmate,nmate,(rtemp,i=1,mmate*nmate)
c        read(23) mmate,nmate,(rtemp,i=1,mmate*nmate)
        deallocate(ipool)
        kelem = kelem+num*nnode
      enddo
      if(iDBG.eq.1 .and. (mype.eq.ichk .or. mype.eq.ichk+1) ) then
        print *,mype, ',num_ef_max = ', num_ef_max
        write(*,'(i3,a,E12.5)') mype,',first read time: ', t_read
      endif
c      if( mype .eq. 0) close(23)
c      close(23)
c
c      if( mype .eq. 0) OPEN (23,FILE=elmfile,FORM='binary',STATUS='OLD')
c      OPEN (23,FILE=elmfile,FORM='binary',STATUS='OLD')
c 
      DO 350 I=1,KNODE+1
        NUMCOL(I)=0
350   CONTINUE
      JNA = KNODE
c 
      KELEM=0
      kelemg = 0
      KEMATE=0
      DO 2000 ITYP=1,NUMTYP
        num_ef = 0
C.......INPUT ENODE
c        if( mype .eq. 0) READ(23) NUM,NNODE
c        call MPI_Bcast(num, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
c        call MPI_Bcast(nnode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
c        READ(23) NUM,NNODE
        num = mnodeg(ityp)
        nnode = nnodeg(ityp)
        nne = nnode
        if( material .eq. 'y' .or. material .eq. 'Y') nne = nne-1
        allocate(NODE(num_ef_max*nnode), STAT = AllocateStatus)
        IF (AllocateStatus .NE. 0) then
          write(*,'(i4,a,3i10)') mype,',num_ef_max=',num_ef_max,nnode
          STOP "*** Not enough memory for NODE ***"
        endif
        isize=nnode
        allocate(ipool(isize), STAT = AllocateStatus)
        do i=1,num
          T0=MPI_WTime()
c          READ(23) (ipool(j),j=1,nnode)
          do j=1, nne
c            ipool(j) = nodeg(i,j)
            ipool(j) = node_g(kelemg+(i-1)*nnode+j)
          enddo
c          ipool(nnode) = me(i)
          ipool(nnode) = node_g(kelemg+i*nnode)
c          if( mype .eq. 0) READ(23) (ipool(j),j=1,nnode)
c          call MPI_Bcast(ipool(1), nnode, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          T1=MPI_WTime()
          T_READ=T_READ+t1-t0
c          if(num .lt. num_ef_max) then
c            allocate(NODE(num_ef_max*nnode), STAT = AllocateStatus)
c            IF (AllocateStatus .NE. 0) then
c              write(*,'(i4,a,3i10)') mype,',num_ef_max=',num_ef_max,nnode
c              STOP "*** Not enough memory for NODE ***"
c            endif
c          else
c            allocate(NODE(num*nnode), STAT = AllocateStatus)
c            IF (AllocateStatus .NE. 0) then
c              write(*,'(i4,a,3i10)') mype,',num=',num,nnode
c              STOP "*** Not enough memory for NODE ***"
c            endif
c          endif
c        WRITE(*,*) 'NUM =',NUM,' NNODE =',NNODE
c          num_ef = 0
          do j=1,nne
            nodi = ipool(j)
            if(nodi .ge. nod_ist .and. nodi .le. nod_iend) goto 100
          enddo
          goto 110
100       continue
c          if(mype.eq.ichk) then
c            print *,'num_ef =', num_ef
c          endif
          do j=1,nnode
            node(num_ef*nnode+j) = ipool(j)
          enddo
          num_ef = num_ef + 1
          if(num_ef .gt. num_ef_max) then
            print *,'num_ef gt num_ef_max!!!', num_ef, num_ef_max
            goto 9999
          endif
110       continue
        enddo
        deallocate(ipool, STAT = DeAllocateStatus)
        IF (DeAllocateStatus .NE. 0)
     +    STOP "** DeAllocate ipool Failed **"
        T0=MPI_WTime()
c        if( mype .eq. 0) read(23) mmate,nmate,(rtemp,i=1,mmate*nmate)
c        read(23) mmate,nmate,(rtemp,i=1,mmate*nmate)
        T1=MPI_WTime()
        T_READ=T_READ+t1-t0
        if(iDBG.eq.1 .and. mype.eq.ichk) then
          write(*,*) 'mype is', mype, ', num_ef =', num_ef
        endif
c
c        WRITE(*,*) 'NODE ='
c        WRITE(*,6) ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
c        IF (KELEM .LT. NUM*NNODE) KELEM = NUM*NNODE
        NNE = NNODE
        IF (MATERIAL.EQ.'Y' .OR. MATERIAL.EQ.'y') THEN
          NNE = NNE-1
        ENDIF
cc        WRITE(*,*) 'NUM =',NUM,' NNODE =',NNODE
CC        WRITE(*,*) 'NODE ='
CC        WRITE(*,6) ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
        T0=MPI_WTime()
        DO 1000 NE=1,NUM_ef
          L=0
          DO 700 INOD=1,NNE
            NODI=NODE((NE-1)*NNODE+INOD)
            L=L+1
            LM(L)=NODI
700       CONTINUE
c          if(mype.eq.ichk) then
c            WRITE (*,*) 'L,LM =',L
c            WRITE (*,'(1X,15I5)') (LM(I),I=1,L)
c          endif
c
          IF (L.GT.0) CALL GBANDWD(ityp,MAXBAND,JNA,NUMCOL,NAP,NAJ,
     +                             L,LM,nod_ist,nod_iend)
          IF (JNA.GT.MAXT) THEN
            WRITE(*,*) 'EXCEET ARRAY LENGTH MAXT ....',MAXT,' < ',JNA
            STOP 1111
          ENDIF
1000    CONTINUE
        T1=MPI_WTime()
        T_M2G = T_M2G + T1 - T0
        deallocate(NODE, STAT = DeAllocateStatus)
        IF (DeAllocateStatus .NE. 0) STOP "** DeAllocate NODE Failed **"
        kelemg = kelemg+num*nnode
2000  CONTINUE
c      CLOSE(23)
c      if( mype .eq. 0) CLOSE(23)
c      IF (DeAllocateStatus .NE. 0) STOP "** DeAllocate ipool Failed **"
c
c      if(mype.eq.ichk) then
c        print *,'JNA = ', JNA
c        print *,'NUMCOL = ', (NUMCOL(I),I=1,KNODE_t+1)
c      endif
      maxa = 0
      do i=1,knode+1
        maxa = maxa + numcol(i)
      enddo
      if(iDBG.eq.1 .and. mype.eq.ichk) then
        print *,'mype=', mype, ',maxa = ', maxa,',knode = ', knode
      endif
c
      T0=MPI_WTime()
      allocate(NA(MAXA), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "*** Not enough memory For NA***"
c
      CALL GRAPH(KNODE,NAP,NAJ,NUMCOL,NA,LM)
c      write(*,'(a)') 'GRAPH ok...........'
c
      deallocate(LM, STAT = DeAllocateStatus)
      IF (DeAllocateStatus .NE. 0) STOP "** DeAllocate LM Failed **"
      deallocate(NAP, STAT = DeAllocateStatus)
      IF (DeAllocateStatus .NE. 0) STOP "*** DeAllocate NAP Failed ***"
      deallocate(NAJ, STAT = DeAllocateStatus)
      IF (DeAllocateStatus .NE. 0) STOP "*** DeAllocate NAJ Failed ***"

      if(iDBG.eq.1 .and. mype .eq. ichk) then
        write(*,'(a)') 'Generateing Graph finished, begin part....'
        print *,'mype=', mype, ',maxa = ', maxa,',knode = ', knode
      endif
c
      T1 = MPI_WTime()
      T_M2G = T_M2G + T1 - T0
c
c      NEDGS = NUMCOL(knode+1)
      if(iDBG.eq.1 .and. mype.eq.ichk) then
        print *,'MAXA = ', maxa
c        print *,'numcol =', (numcol(i),i=1,knode+1)
c        print *,'NA = ', (NA(i), i=1, NEDGS)
      endif
c      write(*,'(a,E15.8,a)') 'Read elm file elasp time:', T_READ, ' s'
c      write(*,'(a,E15.8,a)') 'Mesh trans to Graph elasp time:', 
c     +                       T_M2G, ' s'
c
      do j=1, maxa
        na(j) = na(j) -1
      enddo
c      
      allocate(ipart(knode), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "*** Not enough memory ***"
c
      T0=MPI_WTime()
      call partsub0(vtxdist, numcol, na, ipart)
      T1=MPI_WTime()
      T_PART=T_PART+T1-T0
c
6     format(8i10)
c
      deallocate(NA, STAT = DeAllocateStatus)
      IF (DeAllocateStatus .NE. 0) STOP "*** DeAllocate NA Failed ***"
      deallocate(NUMCOL, STAT = DeAllocateStatus)
      IF (DeAllocateStatus .NE. 0) STOP "** DeAllocate NUMCOL Failed **"

      if( iDBG.eq.1 .and. mype.eq.ichk ) 
     +  write(*,*) 'Begin reduce and assemble partid...'
      T0=MPI_WTime()
C.... Reduce all ipart to ipart_all, and save to part.idx file
      T10=MPI_WTime()
      if ( mype .eq. 0) then
        allocate(ipart_all(knode_t), STAT = AllocateStatus)
        IF (AllocateStatus .NE. 0) STOP "*** Not enough memory ***"
        do i=vtxdist(mype+1)+1, vtxdist(mype+2)
          ipart_all(i) = ipart(i)
        enddo
        T110 = MPI_WTime()
        do ipe=2,npes
          nvtxs = vtxdist(ipe+1)-vtxdist(ipe)
          allocate(ipart_dist(nvtxs), STAT = AllocateStatus)
          IF (AllocateStatus .NE. 0) STOP "*** Not enough memory ***"
          call MPI_Recv(ipart_dist,nvtxs,MPI_Integer8,ipe-1,10000,
     +                  mycomm,istatr,ierr)
          do i=1,nvtxs
            ipart_all(i+vtxdist(ipe)) = ipart_dist(i)
          enddo
          deallocate(ipart_dist, STAT = DeAllocateStatus)
          IF (DeAllocateStatus .NE. 0)
     +      STOP "** DeAllocate ipart_dist Failed **"
        enddo
        T111 = MPI_WTime()
        if( iDBG.eq.1 ) 
     +    write(*,'(i4,a,E12.5,a)') mype,' Processor Recv ipart:',T111-T110,'s'
        T110 = MPI_WTime()
        open(21,file='part.idx',form='formatted',status='unknown')
        do i=1,knode_t
          write(21,'(i7)') ipart_all(i)
        enddo
        close(21)
        open(21,file='partb.idx',form='binary',status='unknown')
        write(21) (ipart_all(i),i=1,knode_t)
        close(21)
        T111 = MPI_WTime()
c        write(*,'(i4,a,E12.5,a)') mype,
c     +       ' processor output part.idx file:',T111-T110,'s'
        T110 = MPI_WTime()
        open(21,file='unodid',form='binary',status='unknown')
        write(21) (dble(ipart_all(i)),i=1,knode_t)
        close(21)
        T111 = MPI_WTime()
c        write(*,'(i4,a,E12.5,a)') mype,' processor output unodid file:',T111-T110,'s'
        deallocate(ipart_all, STAT = DeAllocateStatus)
        IF (DeAllocateStatus .NE. 0)
     +      STOP "** DeAllocate ipart_all Failed **"
      else
        nvtxs = vtxdist(mype+2)-vtxdist(mype+1)
        call MPI_Send(ipart,nvtxs,MPI_Integer8,0,10000,
     +                mycomm, ierr)
      endif
      T11=MPI_WTime()
      if(iDBG.eq.1 .and. mype.eq.0) 
     +   write(*,'(i4,a,E12.5,a)') mype,',Reduce and out:',T11-T10,'s'

      T10=MPI_WTime()
      allocate(Noddist(npes+1), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "*** Not enough memory ***"
      allocate(Noddist0(npes+1), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "*** Not enough memory ***"

      do i=1, npes+1
        Noddist(i) = 0
        Noddist0(i)= 0
      enddo
      do i=1,knode
        ip = ipart(i)+1
        Noddist0(ip) = Noddist0(ip) + 1
      enddo
      T20=MPI_WTime()
      call MPI_Reduce(Noddist0(1),Noddist(1),npes,MPI_Integer,MPI_SUM,
     +                0,mycomm,ierr)
      T21=MPI_WTime()
      if( iDBG.eq.1 .and. mype.eq. 0)
     +   write(*,'(i4,a,E12.5)') mype,',Reduce: ',T21-T20
      T_Comm=T_Comm+T21-T20
c
      T20=MPI_WTime()
      call MPI_Bcast(Noddist(1),npes,MPI_Integer,
     +               0,mycomm,ierr)
      T21=MPI_WTime()
      if( iDBG.eq.1 .and. mype.eq. 0)
     +   write(*,'(i4,a,E12.5)') mype,',Bcast: ',T21-T20
      T_Comm=T_Comm+T21-T20

C...
C.... Assemble nodg
C...
      do i=npes+1,2,-1
        Noddist0(i) = Noddist0(i-1)
      enddo
      Noddist0(1) = 0
      do i=2,npes+1
        Noddist0(i) = Noddist0(i)+Noddist0(i-1)
      enddo
C.... Allocate ipool array to save the global node number
      allocate(ipool(knode), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "* Not enough memory for ipool *"
      allocate(nodnum(npes), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "* Not enough memory for nodnum *"

      do ip=1,npes
        nodnum(ip) = 0
      enddo
      do i=1,knode
        ip = ipart(i)+1
        nod_g = vtxdist(mype+1)+i
        nodnum(ip) = nodnum(ip) + 1
        icl = nodnum(ip) + noddist0(ip)
        if( icl .gt. knode) then
        endif
        ipool(icl) = nod_g
      enddo
c
      deallocate(Nodnum, STAT = DeAllocateStatus)
      IF (DeAllocateStatus .NE. 0) STOP "* DeAllocate Nodnum Failed *"
C..
      knode_i = Noddist(mype+1)
      if( knode_i .gt. maxknode_i) then
        write(*,'(a)') 'Error! maxknode_i too small!!'
        write(*,'(i6,2(a,i7))') mype,',knode_i=',knode_i,'maxknode_i=',maxknode_i
      endif
c      print *,'mype=', mype, ', knode_i =', knode_i, ', npes=', npes
c      allocate(Nodg(knode_i), STAT = AllocateStatus)
c      IF (AllocateStatus .NE. 0) STOP "* Not enough memory for nodg *"

      inod = 0
      do ip=1, npes
        ipe = ip-1
        if( ipe .lt. mype ) then
C....   Mype proccessor Send nodipe to other proccessor
          nodipe = noddist0(ip+1)-noddist0(ip)
C..........................................................
          call MPI_Send(nodipe,1,MPI_Integer,ipe,1000,mycomm,ierr)
          if(nodipe .gt. 0 ) then
            nodipe0 = noddist0(ip)+1
c            call MPI_Send(ipool(nodipe0),nodipe,MPI_Integer,
c     +                    ipe,1000000,mycomm,ierr)
            T20=MPI_WTime()
            call My_sendai(ipe,mype,ipool(nodipe0),nodipe)
            T21=MPI_WTime()
            T_Comm=T_Comm+T21-T20
          endif
C...........................................................
C....   Mype proccessor Recv noddip0 from other proccessor
          call MPI_Recv(nodipe0,1,MPI_Integer,ipe,1000,mycomm,istatr,ierr)
          if( nodipe0 .gt. 0) then
            allocate(nodnum(nodipe0), STAT = AllocateStatus)
            IF (AllocateStatus .NE. 0) 
     +         STOP "* Not enough memory for nodnum *"
c            call MPI_Recv(nodnum(1),nodipe0,MPI_Integer,
c     +                    ipe,1000000,mycomm,istatr,ierr)
            T20=MPI_WTime()
            call My_Recvai(mype,ipe,nodnum(1),nodipe0)
            T21=MPI_WTime()
            T_Comm=T_Comm+T21-T20
            do j=1,nodipe0
              inod = inod+1
              Nodg(inod) = nodnum(j)
            enddo
            deallocate(Nodnum, STAT = DeAllocateStatus)
            IF (DeAllocateStatus .NE. 0) 
     +        STOP "* DeAllocate Nodnum Failed *"
          endif
        elseif(ipe .gt. mype) then
C....     Mype proccessor Recv noddip0 from other proccessor
          call MPI_Recv(nodipe0,1,MPI_Integer,ipe,1000,mycomm,istatr,ierr)
          if( nodipe0 .gt. 0) then
c            write(*,'(2(i2,a),i10)') mype, ' <==== ', ipe,' ..',nodipe0
            allocate(nodnum(nodipe0), STAT = AllocateStatus)
            IF (AllocateStatus .NE. 0) 
     +        STOP "* Not enough memory for nodnum *"
c            call MPI_Recv(nodnum(1),nodipe0,MPI_Integer,
c     +                    ipe,1000000,mycomm,istatr,ierr)
            T20=MPI_WTime()
            call My_Recvai(mype,ipe,nodnum(1),nodipe0)
            T21=MPI_WTime()
            T_Comm=T_Comm+T21-T20

            do j=1,nodipe0
              inod = inod+1
              Nodg(inod) = nodnum(j)
            enddo
            deallocate(Nodnum, STAT = DeAllocateStatus)
            IF (DeAllocateStatus .NE. 0) 
     +        STOP "* DeAllocate Nodnum Failed *"
          endif
C....     Mype proccessor Send nodipe to other proccessor
          nodipe = noddist0(ip+1)-noddist0(ip)
C..........................................................
          call MPI_Send(nodipe,1,MPI_Integer,ipe,1000,mycomm,ierr)
          if(nodipe .gt. 0 ) then
            nodipe0 = noddist0(ip)+1
c            call MPI_Send(ipool(nodipe0),nodipe,MPI_Integer,
c     +                    ipe,1000000,mycomm,ierr)
            T20=MPI_WTime()
            call My_Sendai(ipe,mype,ipool(nodipe0),nodipe)
            T21=MPI_WTime()
            T_Comm=T_Comm+T21-T20

c            write(*,'(2(i2,a),10i10)') Mype, ' ====> ', ipe,' ..', nodipe
          endif
C...........................................................
        else if(ipe .eq. mype) then
          nodipe = noddist0(ip+1)-noddist0(ip)
          do j=1,nodipe
            inod=inod+1
            Nodg(inod) = ipool(noddist0(ip)+j)
          enddo
        else
          write(*,'(a,3i10)') 'Error!',ipe,npes,mype
        endif
      enddo
c      print *,'mype=',mype,',inod=',inod,',knode_i=',knode_i
      T20=MPI_WTime()
      write(ext,'(i5)') mype
      xyzfile= 'nodg_' // trim(adjustl(ext))
c      open(21,file=xyzfile,form='formatted',status='unknown')
c      do i=1,knode_i
c        write(21,'(i15)') Nodg(i)
c      enddo
c      close(21)
      open(21,file=xyzfile,form='binary',status='unknown')
      write(21) (Nodg(i),i=1,knode_i)
      close(21)
      T21=MPI_WTime()
c      if(mype .eq. 0) write(*,'(i4,a,E12.5)') mype,'Write nodg file: ',T21-T20

      deallocate(ipool, STAT = DeAllocateStatus)
      IF (DeAllocateStatus .NE. 0) STOP "* DeAllocate ipool Failed *"
      deallocate(Noddist, STAT = DeAllocateStatus)
      IF (DeAllocateStatus .NE. 0) STOP "* DeAllocate Noddist Failed *"
      deallocate(Noddist0, STAT = DeAllocateStatus)
      IF (DeAllocateStatus .NE. 0) STOP "* DeAllocate Noddist0 Failed *"
      deallocate(vtxdist, STAT = DeAllocateStatus)
      IF (DeAllocateStatus .NE. 0) STOP "* DeAllocate vtxdist Failed *"
      deallocate(ipart, STAT = DeAllocateStatus)
      IF (DeAllocateStatus .NE. 0) STOP "* DeAllocate ipart Failed *"
c      deallocate(Nodg, STAT = DeAllocateStatus)
c      IF (DeAllocateStatus .NE. 0) STOP "* DeAllocate Nodg Failed *"
      T11=MPI_WTime()
      if(iDBG.eq.1 .and. mype.eq.0)
     +  write(*,'(i4,a,E12.5,a)') mype,',Assm and out nodg:',T11-T10,'s'
      if(iDBG.eq.1 .and. mype.eq.0)
     +  write(*,'(i4,a,E12.5)') mype,',Comm:',T_Comm
      T1=MPI_WTime()
      T_OUT=T_OUT+T1-T0
c
9999  continue 
c      print *,'Mesh2Graph Finished !', mype, ' proccessor...'
      call MPI_COMM_Free(mycomm, ierr)
C
      if(iDBG.eq.1 .and. mype.eq.0) then
        write(*,'(i4,5(a,E12.5),a)') mype,',Read file: ', T_READ,
     +                                's,M2G: ',T_M2G,'s,Part: ',T_PART,
     +                                's,OUT:',T_OUT,'s,Comm:',T_Comm,
     +                                's'
      endif

      return
      end

      SUBROUTINE GBANDWD(ityp,MAXBAND,JNA,NUMCOL,NAP,NAJ,ND,LM,ist,iend)
      IMPLICIT REAL*8 (A-H,O-Z)
c      IMPLICIT INTEGER*8(I-N)
      INTEGER JNA, JP
      INTEGER*8 NAP(*),NUMCOL(*),NAJ(*)
      INTEGER LM(*)
6     FORMAT(1X,5I15)
c      WRITE (*,*) 'ND= ',ND, (LM(I),I=1,ND)
      DO 400 I=1,ND
        NI1 = LM(I)
        if(ni1.lt.ist .or. ni1.gt.iend) goto 400
        NI = NI1 - ist + 1
        DO 300 J=1,ND
          NJ = LM(J)
          IF (NJ.EQ.NI1) GOTO 300
          NUMJ = NUMCOL(NI)
          IF (NUMJ.EQ.0) THEN
            JNA = JNA+1
            NAP(NI) = JNA
            NAJ(NI) = NJ
            NUMCOL(NI) = 1
          ELSE
            JP = NAP(NI)
            JV = NAJ(NI)
            IF (NJ.EQ.JV) GOTO 300
            DO K=1,NUMJ-1
              JV = NAJ(JP)
              JP = NAP(JP)
              IF (NJ.EQ.JV) GOTO 300
            ENDDO
            JNA = JNA+1
            NAP(JP) = JNA
            NAJ(JP) = NJ
            NUMCOL(NI) = NUMCOL(NI)+1
            if(NUMCOL(NI) .gt. MAXBAND) then
              write(*,*) 'Warring!!NUMCOL GT NBIND:',ityp,NI,NUMCOL(NI),MAXBAND
              WRITE (*,*) 'ND= ',ND, (LM(iI),iI=1,ND)
            endif
          ENDIF
300     CONTINUE
400   CONTINUE
      RETURN
      END
 
      SUBROUTINE GRAPH(NEQ,NAP,NAJ,NUMCOL,NA,LMI)
      IMPLICIT REAL*8 (A-H,O-Z)
c      IMPLICIT INTEGER*8(I-N)
      INTEGER*8 NUMCOL(*), NAP(*),NAJ(*),NA(*)
      INTEGER JP, NN
      INTEGER LMI(*)

      NN = 0
      DO 600 N=1,NEQ
        JP = NAP(N)
        JV = NAJ(N)
        LI = NUMCOL(N)
        DO 500 I=1,LI
          LMI(I) = JV
          JV = NAJ(JP)
          JP = NAP(JP)
500     CONTINUE
        CALL ORDER(LI,LMI)
        DO I=1,LI
          NN = NN+1
          NA(NN) = LMI(I)
        ENDDO
600   CONTINUE
      DO 800 N=1,NEQ-1
800   NUMCOL(N+1) = NUMCOL(N+1)+NUMCOL(N)
      DO 850 N=1,NEQ
850   NUMCOL(NEQ-N+2) = NUMCOL(NEQ-N+1)
      NUMCOL(1) = 0
c      WRITE(*,*) 'NUMCOL ='
c      WRITE(*,6) (NUMCOL(N),N=1,NEQ+1)
c      WRITE(*,*) 'NA ='
c      WRITE(*,6) (NA(I),I=1,NN)
1000  RETURN
6     FORMAT(1X,5I15)
      RETURN
      END

      SUBROUTINE ORDER(ND,LM)
      IMPLICIT REAL*8 (A-H,O-Z)
c      IMPLICIT INTEGER*8(I-N)
      INTEGER LM(*)
C       WRITE(*,*) '**** ORDER ****'
C       WRITE(*,*) (LM(I),I=1,ND)
      DO 200 I=1,ND
      LS=LM(I)+1
      DO 100 J=I,ND
      IF (LM(J).GT.LS) GOTO 100
      LS=LM(J)
      J0=J
100     CONTINUE
      LM(J0)=LM(I)
      LM(I)=LS
200     CONTINUE
C       WRITE(*,*) (LM(I),I=1,ND)
C       WRITE(*,*) '-----------------'
      RETURN
      END
