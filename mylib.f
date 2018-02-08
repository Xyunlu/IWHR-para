      subroutine My_init(myrank,numblk)
      implicit real*8 (a-h,o-z)
      common/messg/numcpus,myid,ndomains,msids(2,0:2048),msidm(2),
     +             iprint,mycomm,Message_length,MaxMsgTag
      logical filexist
      integer ierr
      include 'mpif.h'
c
c     initial message common block
c
      do j = 1,2
        msidm(j) = 0
        do i=0,2048
          msids(j,i) = 0
        end do
      end do
      call MPI_INIT( ierr )
      call MPI_COMM_DUP( MPI_COMM_WORLD, mycomm, ierr)
      call MPI_COMM_RANK( mycomm, myid, ierr )
      call MPI_COMM_SIZE( mycomm, numnodes, ierr )
      numcpus = numnodes
      ndomains = numcpus
c
      myrank = myid
      numblk = numnodes
c
      Message_length = 1024
      MaxMsgTag=536870000
      iprint = 0

      return
      end

      subroutine My_endjob(ierr)
      implicit real*8 (a-h,o-z)
      common/messg/numcpus,myid,ndomains,msids(2,0:2048),msidm(2),
     +             iprint,mycomm,Message_length,MaxMsgTag
      integer  ierr,numrank,mype
      include 'mpif.h'
c
      call MPI_COMM_RANK(mycomm,mype,ierr)
      call MPI_COMM_FREE(mycomm, ierr)
      call MPI_FINALIZE(ierr)
      if(mype.eq.0)
     +  print *,'The process has been successfully finalized !!'
      return
      end

      subroutine My_abort(ierr)
      implicit real*8 (a-h, o-z)
      common/messg/numcpus,myid,ndomains,msids(2,0:2048),msidm(2),
     +             iprint,mycomm,Message_length,MaxMsgTag
      integer ierr,ierrcode
      include 'mpif.h'
c
      call MPI_ABORT(mycomm, ierrcode, ierr)
      call MPI_COMM_FREE(mycomm, ierr)
      call MPI_FINALIZE(ierr)
      print *,"The process has been stopped with some error!!"
      return
      end

      subroutine My_sendar(idest,isrc,rarray_s,nal)
      implicit real*8 (a-h,o-z)
      common/messg/numcpus,myid,ndomains,msids(2,0:2048),msidm(2),
     +            iprint,mycomm,Message_length,MaxMsgTag
      dimension  rarray_s(1)
      include 'mpif.h'
      dimension nstat(MPI_STATUS_SIZE)
c
c      check parameters 
c
c       Message_length = 1000

      if(idest.lt.0) then
        write(*,*) 'Fatal error, destination id error when mpisend...'
        call My_endjob(ierr)
      end if
      if(isrc.lt.0) then
      write(*,*) 'Fatal error, source id error when mpisend...'
      call My_endjob(ierr)
      end if
      if(isrc.eq.idest) then
      write(*,*) 'Fatal error, source id is same with',
     +           ' dest id when mpisend...'
      call My_endjob(ierr)
      end if
      if((nal.lt.0).or.(nal.gt.100000000)) then
      write(*,*) 'Fatal error, array longth error when mpisend...'
      call My_endjob(ierr)
      end if

c
c      get message tag first
c
c      call timer(4,1)
      ns = isrc
      nd = idest
      nsource = isrc
      ndest = idest
      if(nd.eq.0) then
      mtimes = msidm(1)
      else
      mtimes = msids(1,ndest)
      endif
      nworkers = numcpus
      nsender = isrc
      nrecver = idest
      if(isrc.gt.nworkers) nsender = mod(nsource,nworkers)
      if(idest.gt.nworkers) nrecver = mod(ndest,nworkers)
c
      messagetag = isrc*1000000+idest*1000+mtimes
      if(messagetag .lt. 0) messagetag=-messagetag
      if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
c
      if( iprint .ge. 1 ) then
      write(*,'(a)') 'Start send one double precision message..........'
      write(*,'(a,i5,a,i5)') 'This real message is from ',isrc,' to ',idest
c      write(*,*) 'The source computer node is',isrc
c      write(*,*) 'The destination computer node is',nrecver
      write(*,*) 'Message tag is:',messagetag
      write(*,*) 'This longth of this real array is',nal
      if( iprint .gt. 1) write(*,9000) (rarray_s(i),i=1,nal)
      end if
c
      rnal = nal
      call MPI_SEND(rnal,1,MPI_DOUBLE_PRECISION,nrecver,messagetag,
     &              mycomm,ierr)

c
      if( nal .le. Message_length) then
        mtimes = mtimes + 1
        if(mtimes.gt.997) mtimes = mod(mtimes,997)
        messagetag = isrc*1000000+idest*1000+mtimes
        if(messagetag .lt. 0) messagetag=-messagetag
        if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
        call MPI_SEND(rarray_s(1),nal,MPI_DOUBLE_PRECISION,nrecver,
     &                messagetag,mycomm,ierr)
      else
        do i=1, nal, Message_length
          if( i+Message_length .gt. nal) exit
          mtimes = mtimes + 1
          if(mtimes.gt.997) mtimes = mod(mtimes,997)
          messagetag = isrc*1000000+idest*1000+mtimes
          if(messagetag .lt. 0) messagetag=-messagetag
          if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
          call MPI_SEND(rarray_s(i),Message_length,
     +                  MPI_DOUBLE_PRECISION,nrecver,
     +                  messagetag,mycomm,ierr)
        enddo
        mtimes = mtimes + 1
        if(mtimes.gt.997) mtimes = mod(mtimes,997)
        messagetag = isrc*1000000+idest*1000+mtimes
        if(messagetag .lt. 0) messagetag=-messagetag
        if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
        inal = nal -i +1
        call MPI_SEND(rarray_s(i),inal,MPI_DOUBLE_PRECISION,nrecver,
     +                messagetag,mycomm,ierr)
      endif
c      
      if( iprint .eq. 1 ) then
        write(*,*) 'sent from',isrc, ' to ', idest
        write(*,*) 'double precision one dimension array success!!!'
      end if
c
      mtimes = mtimes + 1
      if(mtimes.gt.997) mtimes = mod(mtimes,997)
      nd = ndest 
      if(nd.eq.0) then
        msidm(1) = mtimes
      else
        msids(1,ndest) = mtimes
      endif
9000  format(6e16.5)
c      call timer(4,2)

      return
      end

      subroutine My_recvar(idest,isrc,rarray_r,nal)
      implicit real*8 (a-h,o-z)
      common/messg/numcpus,myid,ndomains,msids(2,0:2048),msidm(2),
     +             iprint,mycomm,Message_length,MaxMsgTag
      dimension  rarray_r(1)
      include 'mpif.h'
      dimension nstat(MPI_STATUS_SIZE)
c
c      check parameters 
c
c       Message_length = 1000

c       call timer(4,1)
      if(idest.lt.0) then
      write(*,*) 'Fatal error, destination id error when recving...'
      call My_endjob(ierr)
      end if
      if(isrc.lt.0) then
      write(*,*) 'Fatal error, source id error when recving...'
      call My_endjob(ierr)
      end if
      if(isrc.eq.idest) then
      write(*,*) 'Fatal error, source id is same with',
     +           ' dest id when mpisend...'
      call My_endjob(ierr)
      end if
      if((nal.lt.0).or.(nal.gt.100000000)) then
      write(*,*) 'Fatal error, array longth error when recving...'
      call My_endjob(ierr)
      end if
c
c      cleaning buffer
c       
      do i =1,nal
      rarray_r(i) = 0.d0
      end do
c
c      get message tag first
c
      ns = isrc
      nd = idest
      nsource = isrc
      ndest = idest
      if(ns.eq.0) then
      mtimes = msidm(2)
      else
      mtimes = msids(2,ns)
      endif
      nworkers = numcpus
      nsender = nsource
      nrecver = idest
      if(nsource.gt.nworkers) nsender = mod(nsource,nworkers)
      if(idest.gt.nworkers)   nrecver = mod(ndest,nworkers)
c
      messagetag = nsource*1000000+ndest*1000+mtimes
      if(messagetag .lt. 0) messagetag=-messagetag
      if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
c
      if( iprint .eq. 1 ) then
       write(*,'(a)') 'Start receive one double precision message......'
       write(*,'(2(a,i5))') 'This real message is from ',ns,' to ',nd
c       write(*,*) 'The source computer node is',nsender
c       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
      end if
c
      call MPI_RECV(rnals,1,MPI_DOUBLE_PRECISION,nsender,messagetag,
     &              mycomm,nstat,ierr)
c
      nals = int(rnals)
      if(nals.ne.nal) then
       write(*,*) 'Different longth of double precision array received'
       write(*,*) 'Fatal error'
       call My_endjob(ierr)
      end if
c
c       mtimes = mtimes + 1
c       if(mtimes.gt.997) mtimes = mod(mtimes,997)
c       messagetag = nsource*1000000+ndest*1000+mtimes
c
      if( nal .le. Message_length) then
        mtimes = mtimes + 1
        if(mtimes.gt.997) mtimes = mod(mtimes,997)
        messagetag = nsource*1000000+ndest*1000+mtimes
        if(messagetag .lt. 0) messagetag=-messagetag
        if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
        call MPI_RECV(rarray_r(1),nal,MPI_DOUBLE_PRECISION,ns,
     &                messagetag,
     &                mycomm,nstat,ierr)
      else
        do i=1, nal, Message_length
          if( i+Message_length .gt. nal) exit
          mtimes = mtimes + 1
          if(mtimes.gt.997) mtimes = mod(mtimes,997)
          messagetag = nsource*1000000+ndest*1000+mtimes
          if(messagetag .lt. 0) messagetag=-messagetag
          if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
          call MPI_RECV(rarray_r(i),Message_length,
     +                  MPI_DOUBLE_PRECISION,ns,messagetag,
     &                  mycomm,nstat,ierr)
        enddo
        inal = nal-i+1
        mtimes = mtimes + 1
        if(mtimes.gt.997) mtimes = mod(mtimes,997)
        messagetag = nsource*1000000+ndest*1000+mtimes
        if(messagetag .lt. 0) messagetag=-messagetag
        if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
        call MPI_RECV(rarray_r(i),inal,
     +                MPI_DOUBLE_PRECISION,ns,messagetag,
     &                mycomm,nstat,ierr)
      endif
c
      if( iprint .ge. 1 ) then
      write(*,*) 'receive from',nsource, ' to ', ndest
      write(*,*) 'double precision one dimension array success!!!'
      write(*,*) 'This longth of this real array is',nal
      if( iprint .gt. 1) write(*,9000) (rarray_r(i),i=1,nal)
      end if
c
      mtimes = mtimes + 1
      if(mtimes.gt.997) mtimes = mod(mtimes,997)
      nd = ndest
      if(ns.eq.0) then
       msidm(2) = mtimes
      else
       msids(2,ns) = mtimes
      endif
9000  format(6e16.5)
c       call timer(4,2)

      return
      end

      subroutine My_sendai(idest,isrc,iarray_s,nal)
      implicit real*8 (a-h,o-z)
      common/messg/numcpus,myid,ndomains,msids(2,0:2048),msidm(2),
     +            iprint,mycomm,Message_length,MaxMsgTag
      dimension  iarray_s(1)
      include 'mpif.h'
      dimension nstat(MPI_STATUS_SIZE)
c
c      check parameters 
c
c       Message_length = 1000

      if(idest.lt.0) then
      write(*,*) 'Fatal error, destination id error when mpisend...'
      call My_endjob(ierr)
      end if
      if(isrc.lt.0) then
      write(*,*) 'Fatal error, source id error when mpisend...'
      call My_endjob(ierr)
      end if
      if(isrc.eq.idest) then
        write(*,*) 'Fatal error, source id is same with',
     +             ' dest id when mpisend...'
        call My_endjob(ierr)
      end if
      if((nal .le. 0).or.(nal .gt. 100000000)) then
        write(*,*) 'Fatal error, array longth error when mpisend...'
        call My_endjob(ierr)
      end if

c
c      get message tag first
c
c      call timer(4,1)
      ns = isrc
      nd = idest
      nsource = isrc
      ndest = idest
      if(nd.eq.0) then
      mtimes = msidm(1)
      else
      mtimes = msids(1,ndest)
      endif
      nworkers = numcpus
      nsender = isrc
      nrecver = idest
      if(isrc.gt.nworkers) nsender = mod(nsource,nworkers)
      if(idest.gt.nworkers) nrecver = mod(ndest,nworkers)
c
      messagetag = isrc*1000000+idest*1000+mtimes
      if(messagetag .lt. 0) messagetag=-messagetag
      if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
c
      if( iprint .ge. 1 ) then
      write(*,'(a)') 'Start send one integer array message..........'
      write(*,'(a,i5,a,i5)') 'This integer message is from ',isrc,' to ',idest
c      write(*,*) 'The source computer node is',isrc
c      write(*,*) 'The destination computer node is',nrecver
      write(*,*) 'Message tag is:',messagetag
      write(*,*) 'This longth of this integer array is',nal
      if( iprint .gt. 1) write(*,9000) (iarray_s(i),i=1,nal)
      end if
c
      call MPI_SEND(nal,1,MPI_INTEGER,nrecver,messagetag,
     &              mycomm,ierr)

c
      if( nal .le. Message_length) then
        mtimes = mtimes + 1
        if(mtimes.gt.997) mtimes = mod(mtimes,997)
        messagetag = isrc*1000000+idest*1000+mtimes
        if(messagetag .lt. 0) messagetag=-messagetag
        if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
        call MPI_SEND(iarray_s(1),nal,MPI_INTEGER,nrecver,
     &                messagetag,mycomm,ierr)
      else
        do i=1, nal, Message_length
          if( i+Message_length .gt. nal) exit
          mtimes = mtimes + 1
          if(mtimes.gt.997) mtimes = mod(mtimes,997)
          messagetag = isrc*1000000+idest*1000+mtimes
          if(messagetag .lt. 0) messagetag=-messagetag
          if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
          call MPI_SEND(iarray_s(i),Message_length,
     +                  MPI_INTEGER,nrecver,
     +                  messagetag,mycomm,ierr)
        enddo
        mtimes = mtimes + 1
        if(mtimes.gt.997) mtimes = mod(mtimes,997)
        messagetag = isrc*1000000+idest*1000+mtimes
        if(messagetag .lt. 0) messagetag=-messagetag
        if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
        inal = nal -i +1
        call MPI_SEND(iarray_s(i),inal,MPI_INTEGER,nrecver,
     +                messagetag,mycomm,ierr)
      endif
c      
      if( iprint .eq. 1 ) then
        write(*,*) 'sent from',isrc, ' to ', idest
        write(*,*) 'integer one dimension array success!!!'
      end if
c
      mtimes = mtimes + 1
      if(mtimes.gt.997) mtimes = mod(mtimes,997)
      nd = ndest 
      if(nd.eq.0) then
        msidm(1) = mtimes
      else
        msids(1,ndest) = mtimes
      endif
9000  format(6e16.5)
c      call timer(4,2)

      return
      end

      subroutine My_recvai(idest,isrc,iarray_r,nal)
      implicit real*8 (a-h,o-z)
      common/messg/numcpus,myid,ndomains,msids(2,0:2048),msidm(2),
     +             iprint,mycomm,Message_length,MaxMsgTag
      dimension  iarray_r(1)
      include 'mpif.h'
      dimension nstat(MPI_STATUS_SIZE)
c
c      check parameters 
c
c       Message_length = 1000

c       call timer(4,1)
      if(idest.lt.0) then
      write(*,*) 'Fatal error, destination id error when recving...'
      call My_endjob(ierr)
      end if
      if(isrc.lt.0) then
      write(*,*) 'Fatal error, source id error when recving...'
      call My_endjob(ierr)
      end if
      if(isrc.eq.idest) then
      write(*,*) 'Fatal error, source id is same with',
     +           ' dest id when mpisend...'
      call My_endjob(ierr)
      end if
      if((nal.le.0).or.(nal.gt.100000000)) then
      write(*,*) 'Fatal error, array longth error when recving...'
      call My_endjob(ierr)
      end if
c
c      cleaning buffer
c       
      do i =1,nal
      iarray_r(i) = 0.d0
      end do
c
c      get message tag first
c
      ns = isrc
      nd = idest
      nsource = isrc
      ndest = idest
      if(ns.eq.0) then
      mtimes = msidm(2)
      else
      mtimes = msids(2,ns)
      endif
      nworkers = numcpus
      nsender = nsource
      nrecver = idest
      if(nsource.gt.nworkers) nsender = mod(nsource,nworkers)
      if(idest.gt.nworkers)   nrecver = mod(ndest,nworkers)
c
      messagetag = nsource*1000000+ndest*1000+mtimes
      if(messagetag .lt. 0) messagetag=-messagetag
      if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
c
      if( iprint .eq. 1 ) then
       write(*,'(a)') 'Start receive one integer array message......'
       write(*,'(2(a,i5))') 'This real message is from ',ns,' to ',nd
c       write(*,*) 'The source computer node is',nsender
c       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
      end if
c
      call MPI_RECV(nals,1,MPI_INTEGER,nsender,messagetag,
     &              mycomm,nstat,ierr)
c
      if(nals.ne.nal) then
       write(*,*) 'Different longth of integer array received'
       write(*,*) 'Fatal error'
       call My_endjob(ierr)
      end if
c
c       mtimes = mtimes + 1
c       if(mtimes.gt.997) mtimes = mod(mtimes,997)
c       messagetag = nsource*1000000+ndest*1000+mtimes
c
      if( nal .le. Message_length) then
        mtimes = mtimes + 1
        if(mtimes.gt.997) mtimes = mod(mtimes,997)
        messagetag = nsource*1000000+ndest*1000+mtimes
        if(messagetag .lt. 0) messagetag=-messagetag
        if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
        call MPI_RECV(iarray_r(1),nal,MPI_INTEGER,ns,
     &                messagetag,
     &                mycomm,nstat,ierr)
      else
        do i=1, nal, Message_length
          if( i+Message_length .gt. nal) exit
          mtimes = mtimes + 1
          if(mtimes.gt.997) mtimes = mod(mtimes,997)
          messagetag = nsource*1000000+ndest*1000+mtimes
          if(messagetag .lt. 0) messagetag=-messagetag
          if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
          call MPI_RECV(iarray_r(i),Message_length,
     +                  MPI_INTEGER,ns,messagetag,
     &                  mycomm,nstat,ierr)
        enddo
        inal = nal-i+1
        mtimes = mtimes + 1
        if(mtimes.gt.997) mtimes = mod(mtimes,997)
        messagetag = nsource*1000000+ndest*1000+mtimes
        if(messagetag .lt. 0) messagetag=-messagetag
        if(messagetag .gt. MaxMsgTag) messagetag=mod(messagetag,MaxMsgTag)
        call MPI_RECV(iarray_r(i),inal,
     +                MPI_INTEGER,ns,messagetag,
     &                mycomm,nstat,ierr)
      endif
c
      if( iprint .ge. 1 ) then
      write(*,*) 'receive from',nsource, ' to ', ndest
      write(*,*) 'integer one dimension array success!!!'
      write(*,*) 'This longth of this integer array is',nal
      if( iprint .gt. 1) write(*,9000) (iarray_r(i),i=1,nal)
      end if
c
      mtimes = mtimes + 1
      if(mtimes.gt.997) mtimes = mod(mtimes,997)
      nd = ndest
      if(ns.eq.0) then
       msidm(2) = mtimes
      else
       msids(2,ns) = mtimes
      endif
9000  format(6e16.5)
c       call timer(4,2)

      return
      end
