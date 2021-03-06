C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C THREE DIMENSION FINITE ELEMENT PROGRAM OF CFRD EMBANKMENT !
C                  STATIC ANALYSIS                          !
C           BY LHJ(07/9/2016) in IWHR                       !
C                                                           !
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C!!!!!!!!!!!!!!!!!    PROGRAM MAIN   !!!!!!!!!!!!!!!!!!!!!!!!
C  THE FOLLOWING SUBROUTINES ARE CALLED IN THE MAIN PROGRAM !
C    SUBROUTINE NNAMES     DESCRIPTION                      !
C    INPUT                 输入数据                         !
C    EIVI                  计算初始弹模与泊松比             !
C    NHMA                  形成主对角线元素指示阵           !
C    WATER                 计算水荷载                       !
C    ZL25                  累加自重、分解、回代、中点增量   !
C    FAUVW                 缝的错动位移                     !
C    MBYL                  面板应力                         !
C 非零元存储，超松弛共轭梯度迭代法                          !
C 接触单元分为两类:                                         !
C   一类为面板下接触单元AE(9,IM)<1.0,                       !
C   其他的接触单元>1.0                                      !
C 混凝土单元分为两类:                                       !
C   一类为面板单元AE(7,IM)=0.0,AE(9,IM)=0.0,                !
C   其他的混凝凝土单元>1.0                                  !
C 对于接触单元和缝单元，位移为先短边(沉陷变形)，            !
C   后长边(剪切)，最后拉伸和压缩                            !
C 面板纵缝先沉陷(垂直于坝坡向内), 剪切(顺坝坡), 拉伸和压缩  !
C 面板应力先法向，顺坡和坝轴向                              !
C ssorg是超松弛共轭梯度法求解，收敛标准可自行调整在子程序中;!
C 对于土石坝工程而言，10-3一般即可满足要求                  !
C 是否考虑地基施工, 如有地基施工, zl25子程序中初始步应在    !
C   地基应力计算后开始计算位移                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C 防渗墙施工中作为覆盖层材料时的应力和变形要扣除，          !
C   在ZL25中需要更新防渗墙的材料号                          !
C 泊松比的下限对水平向位移影响很大,一般堆石在0.3~0.35,粘土0.4!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PROGRAM pEarthdamSTATIC
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      dimension uvw2(3)
      COMMON /A2/IPE(MAXN,8)                                          ! 单元节点信息
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN                               ! 单元数, 节点数, 节点数
      COMMON /A4/N,NH,MX,JR(3,MAXN)
      COMMON /A5/R(MAXN1)
      COMMON /A7/COP(MAXN,3),AE(11,40)                                ! 节点坐标, 材料矩阵
      COMMON /A8/ME(MAXN)                                             ! 单元材料编号
      COMMON /A11/NET1(100),NRW1(100),NRWD1(100)
      COMMON /A12/RR(MAXN)
      COMMON /A13/ET(MAXN),UT(MAXN),GAMT(MAXN)
      COMMON /A14/UVW(MAXN,3)
      COMMON /A15/STRZ(MAXN,6),EPGZ(MAXN,6)
      COMMON /A16/DZ(MAXN,3)
      COMMON /CMN4/NEE,NET,NRW,NDE
      COMMON /GAUSS/RSTG(3),H(3)
      COMMON /SQ/S0(MAXN),Q0(MAXN)
      COMMON /STRM/SSS(MAXN,3)
      COMMON /NJX/NJJ,MBJ,HDMAX,HDAM(100),HHY(100),II0(100),HHY2(100) !分步加载信息, II0-加载步单元分组
      COMMON /FA/AN,HDZ
      COMMON /EJOINT/IJKL(MAXN),VML(MAXN)
      COMMON /WATERNODE/NUH,IUH(MAXN,4),NDH,IDH(MAXN,4)               !水压单元
      COMMON /POSONGBI/PUMIN,PUMAX
      COMMON /INDEXSHIGONG/NUMBERSHIGONG
      
C=========================================================================C
C     node_glb: 整体网格单元节点编号，包括体单元和面单元
C     mmodeg, nnodeg: 各类型单元个数，各类型单元节点个数+材料编号
C     nodg: 子区域节点对应的整体节点编号 
C     coorg:  整体节点坐标
C=========================================================================C
      integer,dimension(:),allocatable :: nodg, node_glb, mnodeg, nnodeg
      integer            maxknode_i, knode_t, numc8, numq4

      real,dimension(:),allocatable :: coorg

      integer            AllocateStatus, DeAllocateStatus
      double precision   t0,t1,t_msh,t_cor,t_dat,t_sol,tt0,tt1

      double precision time_b, time_e
      
      call cpu_time(time_b)

      call My_init(mype,npes)

      kdgof = 3
      iblk = mype
      numblk = npes
      nblk = numblk

      allocate(idx_neigh(npes))

!     输入文件
      OPEN(15,FILE='INPUT.DAT',STATUS='OLD')
      OPEN(18,FILE='WATERENODE.DAT',STATUS='OLD')
!     输入文件
      
!     0号进程负责输出结果文件，其余进程跳过该段
      if( iblk .gt. 0) goto 1

!     输出结果文件包括位移和应力，竣工期和蓄水期
      OPEN(16,FILE='SOLUTION.DAT',STATUS='UNKNOWN')
      OPEN(17,FILE='EV1.dat',STATUS='UNKNOWN')
!     输出结果文件包括位移和应力，竣工期和蓄水期
      
!     输出结果文件包括应力水平，初始模量
      OPEN(20,FILE='LEVEL.DAT',STATUS='UNKNOWN')
      OPEN(30,FILE='INIMOLIANG.DAT',STATUS='UNKNOWN')
!     输出结果文件包括应力水平，初始模量
      
!     面板应力和面板稳定性
      OPEN(31,FILE='SLAB.dat',STATUS='UNKNOWN')
      OPEN(50,FILE='STABILITY.dat',STATUS='UNKNOWN')
!     面板应力和面板稳定性
      
!     ANSYS出图数据
      OPEN(51,FILE='UVWANSYS.DAT',STATUS='UNKNOWN')
      OPEN(52,FILE='STRESSANSYS.DAT',STATUS='UNKNOWN')
      OPEN(53,FILE='JGUVW.DAT',STATUS='UNKNOWN')
      OPEN(54,FILE='JGSTRESS.DAT',STATUS='UNKNOWN')
!     ANSYS出图数据
      
!     空库或竣工期坝体应力分布
      open(24,FILE='EMPTYSTRESS.DAT',STATUS='UNKNOWN')
!     空库或竣工期坝体应力分布
      
!     用于存储单元体积
      OPEN(55,FILE='TEMP.DAT')
!     用于存储单元体积
      
!     用于动力计算非零元指示矩阵
      OPEN(65,FILE='KMATRIXINDEX.DAT')
      OPEN(75,FILE='MA.DAT')
!     用于动力计算非零元指示矩阵

1     continue        !开始公共数据赋值

      ! 高斯积分权
      H(1)=0.5555555556
      H(2)=0.8888888889
      H(3)=H(1)
      RSTG(1)=-0.774596669241
      RSTG(2)=0.0
      RSTG(3)=-RSTG(1)
      !!!!!!!!!!!!!!!!!!
!     LIMIT THE POSSION TATIO OF SOFT SOIL
      PUMIN=0.3     ! 最小泊松比
      PUMAX=0.49    ! 最大泊松比
      ! LIMIT THE POSSION TATIO OF SOFT SOIL
      ! 施工期前几期坝体不计算变形，防渗墙不计算应力和变形
      NUMBERSHIGONG=0 
      ! 施工期前几期坝体不计算变形，防渗墙不计算应力和变形

      ! 读入数据

      ! ================================================= !
      ! 读入上下游水压力单元
      ! ================================================= !
      READ(18,*)NUH,NDH
      write(*,*)'UPSTREAM WATER ELEMENT',nuh   ! 上游水压面单元
      write(*,*)'DOWNSTREAM WATER ELEMENT',ndh ! 下游水压面单元
      IF(NUH.EQ.0)GOTO 80
      DO I=1,NUH
        READ(18,*) IK,NOELE,(IUH(I,J),J=1,4)
      ENDDO
80    CONTINUE
      
      IF(NDH.EQ.0)GOTO 81
      DO I=1,NDH
        READ(18,*)IK,NOELE,(IDH(I,J),J=1,4)
      ENDDO
81    CONTINUE
      ! ================================================= !
      
      ! ================================================= !
      ! 读入计算模型数据：单元、节点....
      ! ================================================= !
      CALL INPUT(GAMWA) ! GAMWA IS GRAVITY OF WATER 1.0
      ! ================================================= !
      WRITE(*,*) 'INPUT FINISH'
c      write(*,*) 'Mype, NJJ =', Mype, NJJ

      knode_t = NP                                 ! 节点总数
      maxknode_i = knode_t/nblk*2
      if (maxknode_i .lt. knode_t/nblk+1000) maxknode_i=knode_t/nblk+1000
      if (maxknode_i .lt. 20000 ) maxknode_i=20000

      allocate(Nodg(maxknode_i), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "* Not enough memory for nodg *"
c
      numc8 = ne
      nnodec8 = 9
c      numq4 = nuh + ndh
      nnodeq4 = 5
      kelemg = numc8*nnodec8 + nuh*nnodeq4 + ndh*nnodeq4
      allocate(mnodeg(3), nnodeg(3), STAT = AllocateStatus)
      mnodeg(1) = numc8
      nnodeg(1) = nnodec8
c      mnodeg(2) = numq4
      mnodeg(2) = nuh
      nnodeg(2) = nnodeq4
      mnodeg(3) = ndh
      nnodeg(3) = nnodeq4
      allocate(Node_glb(kelemg), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) 
     +  STOP "* Not enough memory for node_glb *"
      kelemg = 0
      do i=1, numc8
        do j=1,nnodec8-1
          node_glb(kelemg+(i-1)*nnodec8+j) = ipe(i,j)
        enddo
        node_glb(kelemg+i*nnodec8) = me(i)
      enddo
      kelemg = kelemg + numc8*nnodec8
      do i=1, nuh
        do j=1,nnodeq4-1
          node_glb(kelemg+(i-1)*nnodeq4+j) = iuh(i,j)
        enddo
        node_glb(kelemg+i*nnodeq4) = 1
      enddo
      kelemg = kelemg + nuh*nnodeq4
      do i=1, ndh
        do j=1, nnodeq4-1
          node_glb(kelemg+(i-1)*nnodeq4+j) = idh(i,j)
        enddo
        node_glb(kelemg+i*nnodeq4) = 1
      enddo
      kelemg = kelemg + ndh*nnodeq4

      ! 网格分区
      call partmesh(maxknode_i,nodg,
     +              mnodeg,nnodeg,node_glb,
     +              knode_i,knode_t)
      
      kcoor = 3
      allocate(coorg(knode_t*kcoor), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) 
     +  STOP "* Not enough memory for node_glb *"
      do i=1,knode_t
        do j=1,kcoor
          coorg((i-1)*kcoor+j) = cop(i,j)
        enddo
      enddo
      call partcoor(maxknode_i,knode_i,nodg,
     +              mnodeg,nnodeg,node_glb,
     +              knode_t,coorg)

      allocate(iN_lg(knode), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "* Not enough memory for iN_lg *"
      do i=1, knode_i
        iN_lg(i) = nodg(i)
      enddo
C     整理扩展节点的整体节点编号
      write(ext,'(i5)') mype
      filename='conn_' // trim(adjustl(ext))
      open(21,file=filename,form='formatted',status='old')
      read(21,*) neigh
      do i=1, neigh
        read(21,*) isend, nsend
        allocate(neighs(nsend))
        read(21,*) (neighs(j),j=1,nsend)
        ntmp = nsend
        allocate(ipool(ntmp))
        do j=1,nsend
          inod = neighs(j)
          if(inod .gt. knode_i .or. inod .le. 0) then
            print *,'Error!! inod, knode=',inod, knode
            call My_endjob(ierr)
          endif
          ipool(j) = iN_lg(inod)
        enddo
        call My_sendai(isend,mype,ipool,ntmp)
        deallocate(neighs, ipool)
c
        read(21,*) irecv, nrecv
        allocate(neighr(nrecv))
        read(21,*) (neighr(j),j=1,nrecv)
        ntmp = nrecv
        allocate(ipool(ntmp))
c        write(*,*) 'mype, irecv, ntmp_r=', mype, irecv, ntmp
        call My_recvai(mype,irecv,ipool,ntmp)
        do j=1, nrecv
          inod = neighr(j)
          if(inod .gt. knode .or. inod .lt. knode_i) then
            print *,'Error!! inod, knode=',inod, knode
            call My_endjob(ierr)
          endif
          if( j .gt. ntmp ) then
            print *,'j, ntmp=',j, ntmp
          endif
          iN_lg(inod) = ipool(j)
        enddo
        deallocate(neighr, ipool)
      enddo
      close(21)
c
      ! ================================================= !
      ! 计算单元体积....
      ! ================================================= !
      CALL SOLVEVOLUME           ! COMPUTE THE VOLUME OF SOLID ELEMENT 
      ! ================================================= !
      WRITE(*,*)'SOLVEVOLUME END'
      ! ================================================= !
      CALL EXAMCONTACTELEMENT    ! 核查接触单元和缝单元的节点排号是否正确
      ! ================================================= !
      WRITE(*,*)'EXAM FOR CONTACT ELEMENT END'
      CLOSE(15)
      
c      DO 30 I=1,NP
      DO 30 I=1,knode
      DO 30 J=1,3
30    UVW(I,J)=0.0
      
      ! 初始化 S0L,Q0L,STRZL,DZL
      num = mnode(1)
      if( .not. allocated(S0L) ) allocate( S0L(num), Q0L(num) )
c      if ( allocated(STRZL) ) deallocate(STRZL, DZL, EPGZL)
      if( .not. allocated(STRZL) ) allocate( STRZL(num,6), DZL(num,3), EPGZL(num, 6) )
c      DO 40 I=1,NE
      DO 40 I=1,num
      Q0L(I)=0.0
      S0L(I)=0.0
      DO 40 J=1,6
      STRZL(I,J)=0.0
      EPGZL(I,J)=0.0
40    CONTINUE

c      DO 50 I=1,NE
      DO 50 I=1,Num
      DO 50 J=1,3
50    DZL(I,J)=0.0
      
      ! ================================================= !
      ! 计算单元初始弹模和泊松比....
      ! ================================================= !
      WRITE(*,*)'INITIAL E G U'
      CALL EIVI  ! 计算初始弹模与泊松比 
      write(*,*)'INITIAL E G U FINISH'
      ! ================================================= !
      
110   FORMAT('============== THE LOAD STEP ============',I8)

      ! ================================================= !
      ! 加载步循环
      ! ================================================= !
c      NJJ = 20
      DO 90 KK=1,NJJ
        if( mype.eq.0 ) THEN
          WRITE(*,110) KK
          WRITE(16,110) KK
        ENDIF
        ! if(mype == 1) WRITE(16,110) KK
        if(mype == 0) WRITE(*,*) KK
        ! 形成结点自由度阵
        CALL MR(KK)  ! 形成节点方程编号，剔除无需加载的单元
        ! 
c        CALL NHMA(KK) ! 形成系数矩阵的存储结构
        CALL NHMA0(KK) ! 形成系数矩阵的DMSR存储结构
60      NET=NET1(KK)
        NERW=NRW1(KK)
        NERWDOWN=NRWD1(KK)
        if( mype.eq.0 ) WRITE(16,190) NET,NERW,NERWDOWN
        if( mype.eq.0 ) WRITE(*,190) NET,NERW,NERWDOWN
190     FORMAT(5X,' NET= ',I4,'       NERW=  ',I4,'       NRWDOWN=  ',I4)
        DO 70 I=1,N
          R(I)=0.0
70      CONTINUE
        IF(NERW.EQ.0 .AND. NERWDOWN.EQ.0) GOTO 999
        IF(NERW.EQ.1) H1=HHY(KK-1)
        IF(NERW.EQ.1) H2=HHY(KK)
        IF(NERW.EQ.1) CALL WATER(H1,H2,NUH,IUH,-1)
        
        IF(NERWDOWN.EQ.1) H3=HHY2(KK-1)
        IF(NERWDOWN.EQ.1) H4=HHY2(KK)
        IF(NERWDOWN.EQ.1) CALL WATER(H3,H4,NDH,IDH,-1)
        
        IF((NERWDOWN.EQ.1).OR.(NERW.EQ.1))THEN
c          DO 122 I=1,NP
          DO 122 I=1,knode
            DO 112 J=1,3
              IV=JRL(J,I)
              IF(IV.EQ.0) UVW2(J)=0.0
              IF(IV.NE.0) UVW2(J)=R(IV)
112         CONTINUE
            DO 16 J=1,3
              IF(UVW2(J).NE.0.0) GO TO 230
16          CONTINUE
            GO TO 122
230         CONTINUE
            if( mype.eq.0 ) WRITE(16,130) iN_lg(I),(UVW2(J),J=1,3)
122       CONTINUE
130       FORMAT(1X,'水荷载I,Rw(I,J)=',I6,3F12.1)
        ENDIF
        
        IF(NERW.EQ.1) WRITE(*,*)'THE LEVEL OF UP AND DOWN WATER',HHY(KK)
        IF(NERW.EQ.1) WRITE(*,100) ! NRW=1,计水荷载
100     FORMAT(1X,'RWATER END')
999     CONTINUE
        ! ZL25  累加自重、分解、回代、中点增量
        CALL ZL25(5,KK) 
        NET=0
        CALL ZL25(10,KK)

        IF(KK.EQ.MBJ)THEN
          ! 整理位移和SSS结果并输出
          if( mype .eq. 0 ) THEN
            write(*,*) '000整理位移并输出到文件'
            allocate(guvw(NP,3))
            do i=1, knode_i
              inod = iN_lg(i)
              do j=1, 3
                guvw(inod,j) = uvw(i,j)
              ENDDO
            ENDDO
            do jpe=1,npes-1
              call my_recvai(0,jpe,knod,1)
              if( allocated(ipool) ) print *,'Error!ipool is allocated!',0
              allocate(ipool(knod))
              call My_recvai(0,jpe,ipool,knod)
              ntmp = knod*3
              if( allocated(rpool) ) print *,'Error!rpool is allocated!',0
              allocate(rpool(ntmp))
              call My_recvar(0,jpe,rpool,ntmp)
              do i=1,knod
                inod = ipool(i)
                do j=1,3
                  guvw(inod,j) = rpool((i-1)*3+j)
                ENDDO
              ENDDO
              deallocate(ipool, rpool)
            ENDDO
            WRITE(53,6)((gUVW(I,J)*100,J=1,3),I=1,NP)
            deallocate(guvw)

            ! 整理单元应力并输出
            write(*,*) '000整理SSS并输出到文件'
            numi = mnode(1)
            do i=1,NumI
              IE = iE_lg(i)
              do j=1,3
                SSS(IE,J) = SSSL(i,j)
              ENDDO
            ENDDO
            do jpe = 1, npes-1
              call My_recvai(0,jpe,Numi,1)
              allocate(ipool(Numi))
              call my_recvai(0,jpe,ipool,Numi)
              allocate(rpool(NumI*3))
              call my_recvar(0,jpe,rpool,Numi*3)
              do i=1,NumI
                IE=ipool(i)
                do j=1,3
                  SSS(IE,j) = rpool((i-1)*3+j)
                enddo
              enddo
              deallocate(ipool, rpool)
            enddo
            DO K=1,NE
              write(24,*) K,(SSS(K,I),I=1,3)
            ENDDO
          ELSE
            !发送位移到0号节点
            call my_sendai(0,mype,knode_i,1)
            if( allocated(ipool) ) print *,'Error!ipool is allocated!',mype
            allocate(ipool(knode_i) )
            do i=1,knode_i
              ipool(i) = iN_lg(i)
            ENDDO
            call my_sendai(0,mype,ipool,knode_i)
            if( allocated(rpool) ) print *,'Error!rpool is allocated!',mype
            allocate(rpool(knode_i*3))
            do i=1, knode_i
              do j=1, 3
                rpool((i-1)*3+j) = uvw(i,j)
              ENDDO
            ENDDO
            call my_sendar(0,mype,rpool,knode_i*3)
            deallocate(ipool,rpool)

            !发送SSS结果到0号节点
            num = mnode(1)
            call my_sendai(0,mype,num,1)
            call my_sendai(0,mype,iE_lg,num)
            allocate(rpool(num*3))
            do i=1,NUM
              do j=1,3
                rpool((i-1)*3+j) = SSSL(i,j)
              ENDDO
            ENDDO
            call my_sendar(0,mype,rpool,num*3)
            deallocate(rpool)
          ENDIF
        
          if( mype .eq. 0 ) THEN  ! 0号进程接收整理和输出应力结果
            write(*,*) '0号进程整理接收和输出应力结果'
            NUM = mnode(1)
            do i=1,num
              IE = iE_lg(i)
              do j=1,6
                STRZ(IE,j) = STRZL(i,j)
              ENDDO
              do j=1,3
                SSS(IE,j) = SSSL(i,j)
              ENDDO
            ENDDO  
            do jpe=1,npes-1
              call My_recvai(0,jpe,Numi,1)
              print *,'0-- jpe,num=',jpe,numi
              allocate(ipool(Numi))
              call My_recvai(0,jpe,ipool,Numi)
              allocate(rpool(Numi*9))   !接收STRZL(numi,6), SSSL(numi,3)
              call My_recvar(0,jpe,rpool,Numi*9)
              DO i=1, NumI
                IE = ipool(i)
                do j=1,6
                  STRZ(IE,j) = rpool((i-1)*9+j)
                ENDDO
                do j=7,9
                  SSS(IE,j-6) = rpool((i-1)*9+j)
                ENDDO
              ENDDO
              deallocate(rpool)
              allocate(rpool(Numi*3))
              call my_recvar(0,jpe,rpool,NumI*3)
              do i=1, NumI
                IE = ipool(i)
                S0(IE) = rpool((i-1)*3+1)
                ET(IE) = rpool((i-1)*3+2)
                UT(IE) = rpool((i-1)*3+3)
              ENDDO
              deallocate(rpool,ipool)
            ENDDO
            DO 118 K=1,NE
              IM=ME(K)
              IF(AE(1,IM).GT.0.3)THEN
                SM1=STRZ(K,1)*0.01
                SM2=STRZ(K,2)*0.01
                SM3=STRZ(K,3)*0.01
                SM4=STRZ(K,4)*0.01
                SM5=STRZ(K,5)*0.01
                SM6=STRZ(K,6)*0.01
                SM7=SSS(K,1)*0.01
                SM8=SSS(K,2)*0.01
                SM9=SSS(K,3)*0.01
              ELSE
                SM1=STRZ(K,1)*0.01
                SM2=STRZ(K,2)*0.01
                SM3=STRZ(K,3)*0.01
                SM4=STRZ(K,4)*10.D0  ! mm -displacement of joint
                SM5=STRZ(K,5)*10.D0
                SM6=STRZ(K,6)*10.D0
                SM7=SSS(K,1)*0.01
                SM8=SSS(K,2)*0.01
                SM9=SSS(K,3)*0.01
              ENDIF
              WRITE(54,7)SM1,SM2,SM3,SM4,SM5,SM6,SM7,SM8,SM9,S0(K),ET(K),UT(K)
118         CONTINUE
          ELSE   !从进程发送应力结果到0号进程
            num = mnode(1)
            call My_sendai(0,mype,num,1)
            call My_sendai(0,mype,iE_lg,num)
            allocate(rpool(num*9))
            do i=1,Num
              do j=1,6
                rpool((i-1)*9+j) = STRZL(i,j)
              ENDDO
              do j=1,3
                rpool((i-1)*9+j+6) = SSSL(i,j)
              ENDDO
            ENDDO
            call My_sendar(0,mype,rpool,num*9)
            deallocate(rpool)
            allocate(rpool(num*3))
            do i=1,num
              rpool((i-1)*3+1) = S0L(i)
              rpool((i-1)*3+2) = ETL(i)
              rpool((i-1)*3+3) = UTL(i)
            ENDDO
            call My_sendar(0,mype,rpool,num*3)
            deallocate(rpool)
          ENDIF
        ENDIF ! 竣工期结束
        
        IF(KK.EQ.NJJ)THEN    !施工期结果整理输出
           ! 整理位移并输出
          if( mype .eq. 0 ) THEN
            write(*,*) '整理位移并输出到文件', knode_i
            allocate(guvw(NP,3))
            do i=1, knode_i
              inod = iN_lg(i)
              do j=1, 3
                guvw(inod,j) = uvw(i,j)
              ENDDO
            ENDDO
            do jpe=1,npes-1
              call my_recvai(0,jpe,knod,1)
              allocate(ipool(knod))
              call My_recvai(0,jpe,ipool,knod)
              ntmp = knod*3
              allocate(rpool(ntmp))
              call My_recvar(0,jpe,rpool,ntmp)
              do i=1,knod
                inod = ipool(i)
                do j=1,3
                  guvw(inod,j) = rpool((i-1)*3+j)
                ENDDO
              ENDDO
              deallocate(ipool, rpool)
            ENDDO
            write(*,*) '输出到文件UVWANSYS.DAT', NP
            WRITE(51,6)((gUVW(I,J)*100.D0,J=1,3),I=1,NP)
            deallocate(guvw)

            ! 整理单元应力并输出
            do K=1,mnode(1)
              IE = iE_lg(k)
              STRZ(IE,1) = STRZL(k,1)
              STRZ(IE,2) = STRZL(k,2)
              STRZ(IE,4) = STRZL(k,4)
              SSS(IE,1) = SSSL(k,1)
              SSS(IE,3) = SSSL(k,3)
            ENDDO
            do jpe=1,npes-1
              call My_recvai(0,jpe,NumI,1)
              allocate(ipool(NumI))
              call My_recvai(0,jpe,ipool,NumI)

              allocate(rpool(NumI))
              call My_recvar(0,jpe,rpool,NumI)
              do K=1,NumI
                IE = ipool(K)
                STRZ(IE,1) = rpool(K)
              ENDDO
              call My_recvar(0,jpe,rpool,NumI)
              do K=1,NumI
                IE = ipool(K)
                STRZ(IE,2) = rpool(K)
              ENDDO
              call My_recvar(0,jpe,rpool,NumI)
              do K=1,NumI
                IE = ipool(K)
                STRZ(IE,4) = rpool(K)
              ENDDO
              call My_recvar(0,jpe,rpool,NumI)
              do K=1,NumI
                IE = ipool(K)
                SSS(IE,1) = rpool(K)
              ENDDO
              call My_recvar(0,jpe,rpool,NumI)
              do K=1,NumI
                IE = ipool(K)
                SSS(IE,3) = rpool(K)
              ENDDO
              deallocate(ipool, rpool)
            ENDDO
            DO K=1,NE
              WRITE(50,91)K,STRZ(K,1),STRZ(K,2),STRZ(K,4),SSS(K,1),SSS(K,3)
            ENDDO
          ELSE
            !发送位移到0号进程
            call my_sendai(0,mype,knode_i,1)
            allocate( ipool(knode_i) )
            do i=1,knode_i
              ipool(i) = iN_lg(i)
            ENDDO
            call my_sendai(0,mype,ipool,knode_i)
            allocate(rpool(knode_i*3))
            do i=1, knode_i
              do j=1, 3
                rpool((i-1)*3+j) = uvw(i,j)
              ENDDO
            ENDDO
            call my_sendar(0,mype,rpool,knode_i*3)
            deallocate(ipool,rpool)

            !发送单元应力结果到0号进程
            num = mnode(1)
            call My_sendai(0,mype,num,1)
            call My_sendai(0,mype,iE_lg,num)
            allocate(rpool(num))
            do i=1,Num
              rpool(i) = STRZL(i,1)
            enddo
            call My_sendar(0,mype,rpool,Num)
            do i=1,Num
              rpool(i) = STRZL(i,2)
            ENDDO
            call My_sendar(0,mype,rpool,Num)
            do i=1,Num
              rpool(i) = STRZL(i,4)
            ENDDO
            call My_sendar(0,mype,rpool,Num)
            do i=1,Num
              rpool(i) = SSSL(i,1)
            ENDDO
            call My_sendar(0,mype,rpool,Num)
            do i=1,Num
              rpool(i) = SSSL(i,3)
            ENDDO
            call My_sendar(0,mype,rpool,Num)
            deallocate(rpool)
          ENDIF
6         FORMAT(3F15.5)
91        FORMAT(1X,I5,2X,5F10.3)
        
        ENDIF ! 正常运行期
      
90    CONTINUE

      if( mype.eq.0 ) THEN
        num=mnode(1)
        do i=1,NUM
          IE = iE_lg(i)
          ET(IE) = ETL(i)
          UT(IE) = UTL(i)
          GAMT(IE) = GAMTL(i)
          do j=1,6
            STRZ(IE,j) = STRZL(i,j)
          ENDDO
          do j=1,3
            SSS(IE,j) = SSSL(i,j)
          ENDDO
        ENDDO
        do jpe=1, npes-1
          call my_recvai(0,jpe,numi,1)
          allocate(ipool(numi))
          call my_recvai(0,jpe,ipool,numi)
          allocate(rpool(numi*3))
          call my_recvar(0,jpe,rpool,numi*3)
          do i=1, numi
            IE = ipool(i)
            ET(IE) = rpool((i-1)*3+1)
            UT(IE) = rpool((i-1)*3+2)
            GAMT(IE) = rpool((i-1)*3+3)
          ENDDO
          deallocate(rpool)
          allocate(rpool(numi*9))
          call My_recvar(0,jpe,rpool,numi*9)
          do i=1,numi
            IE = ipool(i)
            do j=1,6
              STRZ(IE,j) = rpool((i-1)*9+j)
            ENDDO
            do j=1,3
              SSS(IE,j) = rpool((i-1)*9+j+6)
            ENDDO
          ENDDO
          deallocate(ipool, rpool)
        ENDDO
        WRITE(17,220) (K,ET(K),UT(K),GAMT(K),K=1,NE)
        DO 10 K=1,NE
          WRITE(17,210) K,(STRZ(K,IJ),IJ=1,6),(SSS(K,J),J=1,3)
          IM=ME(K)
          IF(AE(1,IM).GT.0.3)THEN
            SM1=STRZ(K,1)*0.01
            SM2=STRZ(K,2)*0.01
            SM3=STRZ(K,3)*0.01
            SM4=STRZ(K,4)*0.01
            SM5=STRZ(K,5)*0.01
            SM6=STRZ(K,6)*0.01
            SM7=SSS(K,1)*0.01
            SM8=SSS(K,2)*0.01
            SM9=SSS(K,3)*0.01
          ELSE
            SM1=STRZ(K,1)*0.01
            SM2=STRZ(K,2)*0.01
            SM3=STRZ(K,3)*0.01
            SM4=STRZ(K,4)*10
            SM5=STRZ(K,5)*10
            SM6=STRZ(K,6)*10
            SM7=SSS(K,1)*0.01
            SM8=SSS(K,2)*0.01
            SM9=SSS(K,3)*0.01
          ENDIF
          WRITE(52,7)SM1,SM2,SM3,SM4,SM5,SM6,SM7,SM8,SM9,S0(K),ET(K),UT(K)
10      CONTINUE
      ELSE
        num = mnode(1)
        call My_sendai(0,mype,num,1)
        call My_sendai(0,mype,iE_lg,num)
        allocate(rpool(num*3))
        do i=1,NUM
          rpool((i-1)*3+1)=ETL(i)
          rpool((i-1)*3+2)=UTL(i)
          rpool((i-1)*3+3)=GAMTL(i)
        ENDDO
        call My_sendar(0,mype,rpool,num*3)
        deallocate(rpool)
        allocate(rpool(num*9))
        do i=1, NUM
          do j=1,6
            rpool((i-1)*9+j) = STRZL(i,j)
          ENDDO
          do j=1,3
            rpool((i-1)*9+j+6) = SSSL(i,j)
          ENDDO
        ENDDO
        call My_sendar(0,mype,rpool,num*9)
        deallocate(rpool)
      ENDIF

210   FORMAT(1X,I8,6F12.2,3F12.2)
220   FORMAT(2(3x,I8,F15.1,F15.3,F15.1))
7     FORMAT(12F15.2)

      if( mype .eq. 0) then
        num = mnode(1)
        do i=1,NUM
          IE=iE_lg(i)
          S0(IE) = S0L(i)
        ENDDO
        do jpe=1,npes-1
          call My_recvai(0,jpe,numi,1)
          allocate(ipool(numi))
          call my_recvai(0,jpe,ipool,numi)
          allocate(rpool(numi))
          call my_recvar(0,jpe,rpool,numi)
          do i=1,numi
            IE = ipool(i)
            S0(IE) = rpool(i)
          ENDDO
          deallocate(ipool,rpool)
        ENDDO
        DO I=1,NE
          WRITE(20,*)I,S0(I)
        ENDDO
        CLOSE(16)
      else
        num = mnode(1)
        call My_sendai(0,mype,num,1)
        call my_sendai(0,mype,iE_lg,Num)
        call My_sendar(0,mype,S0L,num)
      endif
      
99    continue

      call My_endjob(ierr)

      call cpu_time(time_e)
      print *,'Time of operation was ', time_e - time_b, ' seconds'
      
      STOP
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE WATER(H1,H2,NUH,IUH,IDR)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      COMMON /A4/N,NH,MX,JR(3,MAXN)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A5/R(MAXN1)
      COMMON /CMN4/NEE,NET,NRW,NDE
      COMMON /GAUSS/RSTG(3),H(3)
      DIMENSION  FF(4,3),IUH(MAXN,4),inod(4)
      DIMENSION UVW(3)

      !计算水压单元面压力
      write(*,*) 'In Water, num2,nnode2,H1,H2=',mnode(2),nnode(2),NUH,H1,H2
      nskp = mnode(1)*nnode(1)
      num2 = mnode(2)
      nnod2 = nnode(2)
      IF(H2.LE.H1) GOTO 400
c      DO I=1,NUH
      DO I=1,num2
c        X1=COP(IUH(I,1),1)
c        Y1=COP(IUH(I,1),2)
c        Z1=COP(IUH(I,1),3)
        nod1 = node(nskp+(i-1)*nnod2+1)
        inod(1) = NOD1
        X1=COOR((NOD1-1)*3+1)
        Y1=COOR((NOD1-1)*3+2)
        Z1=COOR((NOD1-1)*3+3)
c        X2=COP(IUH(I,2),1)
c        Y2=COP(IUH(I,2),2)
c        Z2=COP(IUH(I,2),3)
        nod2 = node(nskp+(i-1)*nnod2+2)
        inod(2) = nod2
        X2=COOR((NOD2-1)*3+1)
        Y2=COOR((NOD2-1)*3+2)
        Z2=COOR((NOD2-1)*3+3)
c        X3=COP(IUH(I,3),1)
c        Y3=COP(IUH(I,3),2)
c        Z3=COP(IUH(I,3),3)
        nod3 = node(nskp+(i-1)*nnod2+3)
        inod(3) = nod3
        X3=COOR((NOD3-1)*3+1)
        Y3=COOR((NOD3-1)*3+2)
        Z3=COOR((NOD3-1)*3+3)
c        X4=COP(IUH(I,4),1)
c        Y4=COP(IUH(I,4),2)
c        Z4=COP(IUH(I,4),3)
        nod4 = node(nskp+(i-1)*nnod2+4)
        inod(4) = nod4
        X4=COOR((NOD4-1)*3+1)
        Y4=COOR((NOD4-1)*3+2)
        Z4=COOR((NOD4-1)*3+3)
        WP1H1=(H1-Y1)
        IF(WP1H1.LT.0) WP1H1=0
        WP2H1=(H1-Y2)
        IF(WP2H1.LT.0) WP2H1=0
        WP3H1=(H1-Y3)
        IF(WP3H1.LT.0) WP3H1=0
        WP4H1=(H1-Y4)
        IF(WP4H1.LT.0) WP4H1=0
        WP1H2=(H2-Y1)
        IF(WP1H2.LT.0) WP1H2=0
        WP2H2=(H2-Y2)
        IF(WP2H2.LT.0) WP2H2=0
        WP3H2=(H2-Y3)
        IF(WP3H2.LT.0) WP3H2=0
        WP4H2=(H2-Y4)
        IF(WP4H2.LT.0) WP4H2=0
C
        AI=idr*((Y3-Y2)*(Z1-Z2)-(Y1-Y2)*(Z3-Z2)+
     &          (Y4-Y3)*(Z1-Z3)-(Y1-Y3)*(Z4-Z3))/2
        AJ=-idr*(((X3-X2)*(Z1-Z2)-(X1-X2)*(Z3-Z2))+
     &          ((X4-X3)*(Z1-Z3)-(X1-X3)*(Z4-Z3)))/2
        AK=idr*(((X3-X2)*(Y1-Y2)-(X1-X2)*(Y3-Y2))+
     &          ((X4-X3)*(Y1-Y3)-(X1-X3)*(Y4-Y3)))/2
        FF(1,1)=(WP1H2-WP1H1)*AI/4
        FF(1,2)=(WP1H2-WP1H1)*AJ/4
        FF(1,3)=(WP1H2-WP1H1)*AK/4
        FF(2,1)=(WP2H2-WP2H1)*AI/4
        FF(2,2)=(WP2H2-WP2H1)*AJ/4
        FF(2,3)=(WP2H2-WP2H1)*AK/4
        FF(3,1)=(WP3H2-WP3H1)*AI/4
        FF(3,2)=(WP3H2-WP3H1)*AJ/4
        FF(3,3)=(WP3H2-WP3H1)*AK/4
        FF(4,1)=(WP4H2-WP4H1)*AI/4
        FF(4,2)=(WP4H2-WP4H1)*AJ/4
        FF(4,3)=(WP4H2-WP4H1)*AK/4
        DO II=1,4
          DO JJ=1,3
c            IF(JR(JJ,IUH(I,II)).GT.0) THEN
c              R(JR(JJ,IUH(I,II)))=
c     &          R(JR(JJ,IUH(I,II)))+FF(II,JJ)
            IF(JRL(JJ,inod(II)).GT.0) THEN
              R(JRL(JJ,inod(II)))=
     &          R(JRL(JJ,inod(II)))+FF(II,JJ)
            ENDIF
          ENDDO
        ENDDO            
      ENDDO

400   RETURN
      end

!!!!!!!!!!!!!!!!!!!!!!   SUBRUTION ZL25   !!!!!!!!!!!!!!!!!!!!!
!     FUNCTION       累加自重、分解、回代、中点增量
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ZL25(IFO,KK)
      USE ComData
      USE FEMData
      USE solvmodule
      implicit real*8(a-h,o-z)
c      COMMON /A13/ET(MAXN),UT(MAXN),GAMT(MAXN)
      COMMON /A4/N,NH,MX,JR(3,MAXN)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A12/RR(MAXN1)
      COMMON /A5/R(MAXN1)
      COMMON /A1/SK(125930000),SM(MAXN1)
      COMMON /A14/UVW(MAXN,3)
      COMMON /ST/STR(MAXN,6),EPG(MAXN,6)
      COMMON /STRM/SSS(MAXN,3)
      COMMON /A15/STRZ(MAXN,6),EPGZ(MAXN,6)
      COMMON /NJX/NJJ,MBJ,HDMAX,hdam(100),hhy(100),ii0(100),HHY2(100)
      COMMON /A6/MA(MAXN1)
      COMMON /A16/DZ(MAXN,3)
      COMMON /A8/ME(MAXN)
      COMMON /INDEXSHIGONG/NUMBERSHIGONG
      COMMON /A7/COP(MAXN,3),AE(11,40)

c      if( mype.eq.0 ) then
c        write(*,*) '-----------------------------------'
c        write(*,*) 'IN ZL25, IFO,KK=', IFO,KK
c        write(*,*) '-----------------------------------'
c      endif
      CALL KRGAM(KK)    ! 计算刚度矩阵，形成总刚和右端荷载
      if( mype.eq.0 ) WRITE(*,110)
110   FORMAT(1X,'KRGAM END', i3)

      IF(IFO.EQ.10) GO TO 20
      DO 10 I=1,N
10      RR(I)=R(I)

20    DO 30 I=1,N
        IF(IFO.EQ.5)  R(I)=0.5*RR(I)
        IF(IFO.EQ.10) R(I)=RR(I)
30    CONTINUE

      if( mype.eq.0 ) WRITE(*,*)'EQUATION SOLVING STARTING'

      ! 输出分块系数矩阵，测试用
      write(ext,'(i5)') mype
      filename='matr_' // trim(adjustl(ext)) // '.dat'
      open(21,file=filename,form='formatted',status='unknown')
      write(21,*) N_update, na(N_update+1)
      write(21,*) (update(i),i=1,N_update)
      write(21,*) (na(i),i=1, N_update+1)
      do i=1, N_update
        n0 = na(i)+1
        n1 = na(i+1)
        do j=n0,n1
          write(21,'(2i6,E18.10)') update(i), ia(j), am(j)
        enddo
      enddo
      do i=1, N_update
        write(21,*) update(i), R(i)
      enddo
      close(21)

      call DCRS2DMSR(mype,N_update,N_external,
     +              na,ia,am,
     +              update,bindx,val,rhs,sol)

      do i=1, N_update
        rhs(i) = R(I)
      enddo
c      call SSORPCG
      call azsolv(N_update,N_external,update,
     +            bindx, val, rhs, sol)
      
      do i=1,N_update
        R(I) = sol(i)
      enddo
      do i=1, N
        update(i) = update(i) + 1
      enddo
      if( mype.eq.0 ) WRITE(*,*)'EQUATION SOLVING FINISH'
      if( mype == 1) then
        do i=1, 10
          write(*,*) 'sol:',sol(I)
        enddo
      endif
      
      if( allocated(disp) ) deallocate(disp)
      allocate(disp(3,knode))
      do i=1,knode
        do j=1,3
          disp(j,i) = 0.D0
        enddo
      enddo

      !整理位移到节点上
      do i=1,knode_i
        do j=1, 3
          iv = JRL(j,i)
c          if( JRG(j,i).eq.129 .or. JRG(j,i).eq.121 ) then
c            write(*,*) 'Nod_lg, 129/121---mype,inod,nodg=',mype,i,iN_lg(i)
c            write(*,*) 'iv,iw,sol=', iv,JRG(j,i), R(iv)
c          endif
          if( iv > 0 ) disp(j,i) = sol(iv)
        enddo
      enddo

      !交换边界节点位移
      write(ext,'(i5)') mype
      filename='conn_' // trim(adjustl(ext))
      open(21,file=filename,form='formatted',status='old')
      read(21,*) neigh
      do i=1, neigh
        read(21,*) isend, nsend
        allocate(neighs(nsend))
        read(21,*) (neighs(j),j=1,nsend)
        ntmp = nsend*kdgof
c        write(*,*) 'mype, isend, ntmp_s=', mype, isend, ntmp
        allocate(rpool(ntmp))
        do j=1,nsend
          inod = neighs(j)
          if(inod .gt. knode_i .or. inod .le. 0) then
            print *,'Error!! inod, knode=',inod, knode
            call My_endjob(ierr)
          endif
          do k=1,kdgof
            rpool((j-1)*kdgof+k) = disp(k,inod)
          enddo
        enddo
        call My_sendar(isend,mype,rpool,ntmp)
        deallocate(neighs, rpool)
c
        read(21,*) irecv, nrecv
        allocate(neighr(nrecv))
        read(21,*) (neighr(j),j=1,nrecv)
        ntmp = nrecv*kdgof
        allocate(rpool(ntmp))
c        write(*,*) 'mype, irecv, ntmp_r=', mype, irecv, ntmp
        call My_recvar(mype,irecv,rpool,ntmp)
        do j=1, nrecv
          inod = neighr(j)
          if(inod.gt.knode .or. inod.lt.knode_i) then
            print *,'Error!! inod, knode=',inod, knode
            call My_endjob(ierr)
          endif
          do k=1,kdgof
            if( (j-1)*kdgof+k .gt. ntmp ) then
              print *,'(j-1)*kdgof+k, ntmp=',(j-1)*kdgof+k, ntmp
            endif
            disp(k,inod) = rpool((j-1)*kdgof+k)
          enddo
        enddo
        deallocate(neighr, rpool)
      enddo
      close(21)

c      if( mype.eq.0 ) print *,'disp(1894)=',(disp(j,1894),j=1,3)

C...  检查位移是否过大，输出过大的节点编号
      DO 40 I=1,N
        IF(ABS(R(I)) .GT. 10.0)THEN
          WRITE(*,*) R(I)
c          DO J=1,NP
          DO J=1,knode_i
            DO K=1,3
c              IF(JR(K,J) .EQ. I)WRITE(*,*)J
              IF(JRL(K,J) .EQ. I) WRITE(*,*) iN_lg(j)
            ENDDO
          ENDDO
        ENDIF
        IF(ABS(R(I)).GT.10.0) WRITE(*,*)'DEFORMATION IS HUGE'
c        IF(ABS(R(I)).GT.10.0) STOP   
        IF(ABS(R(I)).GT.10.0) call my_endjob(ierr)
40    CONTINUE
      
      IF(IFO.EQ.5) GO TO 70
c      DO 65 I=1,NP
      DO 65 I=1,knode
        DO 60 J=1,3
c          IV=JR(J,I)
          IV=JRL(J,I)
          IF(IV.GT.N) GO TO 65
          IF(IV.EQ.0) UVW(I,J)=0.0
          IF(IV.NE.0)THEN
c            IF((NERWDOWN.EQ.1).OR.(NERW.EQ.1)) UVW(I,J)=UVW(I,J)+1.0*R(IV)
c            IF((NERWDOWN.EQ.0).AND.(NERW.EQ.0)) UVW(I,J)=UVW(I,J)+R(IV)
            IF((NERWDOWN.EQ.1).OR.(NERW.EQ.1)) UVW(I,J)=UVW(I,J)+1.0*disp(j,i)
            IF((NERWDOWN.EQ.0).AND.(NERW.EQ.0)) UVW(I,J)=UVW(I,J)+disp(j,i)
          ENDIF
60      CONTINUE
65    CONTINUE
c      if( mype.eq.1 ) write(*,*) 'JRG(1,1)=', JRG(1,1),',iN_lg(1)=', iN_lg(1)
          if( mype.eq.1 ) write(*,'(a,6E15.7)') 'UVW(29,:)=', (UVW(1,j),j=1,3)

70    continue
c70    CALL STRESS(KK)     ! 求应力
      CALL STRESS(KK)     ! 求应力
      if( mype.eq.0 ) WRITE(*,*)'STRESS END'

      ! 求应力时注意面板是否在施工期，应扣除施工期对面板应力和位移的影响
      Num = mnode(1)

      IF(KK.LE.NUMBERSHIGONG)THEN ! 地基施工不计算位移
c        DO 651 I=1,NP
        DO 651 I=1,knode
          DO 601 J=1,3
            IV=JRL(J,I)
            IF(IV.GT.N) GO TO 651
            IF(IV.EQ.0) UVW(I,J)=0.0
            IF(IV.NE.0) UVW(I,J)=0.0
601       CONTINUE
651     CONTINUE
c        DO 82 I=1,II0(KK)
        DO 82 I=1,num
          if( iE_lg(i) > II0(KK) ) cycle
          DO 82 J=1,6
            IF(IFO.EQ.10) STRZL(I,J)=STRZL(I,J)+STRL(I,J)
            IF(IFO.EQ.5)  STRL(I,J)=STRZL(I,J)+STRL(I,J)
            IF(IFO.EQ.10) STRL(I,J)=STRZL(I,J)
            IF(IFO.EQ.10) EPGZL(I,J)=EPGZL(I,J)+EPGL(I,J)
            IF(IFO.EQ.5)  EPGL(I,J)=EPGZL(I,J)+EPGL(I,J)
            IF(IFO.EQ.10) EPGL(I,J)=EPGZL(I,J)
82      CONTINUE
c        DO 84 I=1,II0(KK)
        DO 84 I=1,num
          IF( iE_lg(I) > II0(KK) ) cycle
          imate = node(I*nnode(1))
c          IF(AE(1,ME(I)).GE.2.4)THEN !FANGSHENQIANG CONSTRUCTION
          IF(AE(1,imate) .GE. 2.4)THEN !FANGSHENQIANG CONSTRUCTION
            DO 85 J=1,6
              IF(IFO.EQ.10) STRZL(I,J)=0.0
              IF(IFO.EQ.5)  STRL(I,J)=0.0
              IF(IFO.EQ.10) STRL(I,J)=0.0
              IF(IFO.EQ.10) EPGZL(I,J)=0.0
              IF(IFO.EQ.5)  EPGL(I,J)=0.0
              IF(IFO.EQ.10) EPGL(I,J)=0.0
85          CONTINUE
          ENDIF
84      CONTINUE
        ! 如果存在地基层，则该层仅计算应力，不将该层得到的位移计入最终位移
      ELSE !KK.NE.1
c        DO 80 I=1,II0(KK)
        DO 80 I=1,num
          if( iE_lg(I) > II0(KK) ) cycle
          DO 80 J=1,6
            IF(IFO.EQ.10) STRZL(I,J)=STRZL(I,J)+STRL(I,J)
            IF(IFO.EQ.5)  STRL(I,J)=STRZL(I,J)+STRL(I,J)
            IF(IFO.EQ.10) STRL(I,J)=STRZL(I,J)
            IF(IFO.EQ.10) EPGZL(I,J)=EPGZL(I,J)+EPGL(I,J)
            IF(IFO.EQ.5)  EPGL(I,J)=EPGZL(I,J)+EPGL(I,J)
            IF(IFO.EQ.10) EPGL(I,J)=EPGZL(I,J)
80      CONTINUE
      ENDIF 
c      if( mype.eq.1 ) then
c        write(*,*) 'IFO =', IFO
c        write(*,*) 'STRZL:', (STRZL(1,J),j=1,6)
c        write(*,*) 'STRL:', (STRL(1,J),j=1,6)
c        write(*,*) 'EPGZL:', (EPGZL(1,J),j=1,6)
c        write(*,*) 'EPGL:', (EPGL(1,J),j=1,6)
c      endif
      
      CALL MAIN(1) !MAIN STRESS
      if( mype.eq.0 ) WRITE(*,*)'MAINSTRESS FINISH'
!     CALL MAIN(0) !MAIN STRAIN AND SHEAR STRAIN
      CALL EVPK(IFO,KK)
      IF(IFO .NE. 10)GOTO 81
      IF(KK.EQ.MBJ.OR.KK.EQ.NJJ) CALL MBYL(KK,0)
!     IF(KK.EQ.MBJ.OR.KK.EQ.NJJ) CALL FAUVW(KK)
      IF(KK.EQ.MBJ.OR.KK.EQ.NJJ) then
        ! 0号进程负责输出结果文件
        if( mype .eq. 0) CALL OUTPUT(IFO,KK)
      ENDIF
!     IF(KK.GT.MBJ.AND.KK.LT.NJJ) CALL OUTPUT1(IFO,KK)
81    CONTINUE

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!    SUBRUTINRE INPUT   !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE INPUT(GAMWA)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /A2/IPE(MAXN,8)
      COMMON /A8/ME(MAXN)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /CMN4/NEE,NET,NRW,NDE
      COMMON /NJX/NJJ,MBJ,HDMAX,hdam(100),hhy(100),ii0(100),HHY2(100)
      COMMON /A10/JC(100001)
      COMMON /A11/NET1(100),NRW1(100),NRWD1(100)
      COMMON /FA/AN,HDZ
      common /iwater/iwate(10000)
      common /iwater2/iwateDOWN(10000)
      COMMON /JLMX/MX
      READ(15,*) NDE
      READ(15,*) NE,NP,NR                            ! 读入单元数,节点数,节点数
      READ(15,*) (IO,(IPE(I,J),J=1,8),ME(I),I=1,NE)  ! 读入c8单元节点号和材料号
      READ(15,*) (IO,(COP(I,J),J=1,3),I=1,NP)        ! 读入节点坐标
      READ(15,*) (JC(I),I=1,NR)                      ! 读入节点约束信息
      READ(15,*) NJJ,MBJ,HDMAX                       ! NJJ: 加载步
      READ(15,*) (HDAM(I),I=1,NJJ)
      READ(15,*) (HHY(I),I=1,NJJ)
      READ(15,*) (HHY2(I),I=1,NJJ)                   ! WATER LEVEL OF DOWNSTREAM SLOPE
      READ(15,*) (II0(I),I=1,NJJ)                    ! 加载步单元分组
      READ(15,*) AN,HDZ
      READ(15,*) GAMWA                               ! 水的重度？
      READ(15,*) NM,MX
      READ(15,*) ((AE(I,J),I=1,11),J=1,NM)
      !    READ(15,*) !AE(1，NM) is or no buoyant unit weight
      READ(15,*) (NET1(I),I=1,NJJ)
      READ(15,*) (NRW1(I),I=1,NJJ)
      READ(15,*) (NRWD1(I),I=1,NJJ) 
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!   SUBRUTINE EIVI  !!!!!!!!!!!!!!!!!!!!!!!
!     FUNCTION            初始弹模与泊松比
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE EIVI
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8)
      COMMON /FA/AN,HDZ
      COMMON /NJX/NJJ,MBJ,HDMAX,hdam(100),hhy(100),II0(100),HHY2(100)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /A2/IPE(MAXN,8)
      COMMON /A8/ME(MAXN)
      COMMON /A13/ET(MAXN),UT(MAXN),GAMT(MAXN)
      COMMON /EJOINT/IJKL(MAXN),VML(MAXN)
      COMMON /POSONGBI/PUMIN,PUMAX
      
      num = mnode(1)

      allocate(etl(num),utl(num), gamtl(num))
c      DO 100 K=1,NE
      DO 100 K=1,Num
c        NM=ME(K)
        NM = node(K*nnode(1))
        GAM=AE(1,NM)
        G=AE(7,NM) ! KB
        F=AE(8,NM) ! M
        D=AE(9,NM) ! DISTINGUISH THE GENERAL CONCRETE AND SLAB ELEMENT
        IF(GAM .ge. 2.4) GOTO 80 ! CONCRETE ELEMENT
      
        YY=0.0                
                    
        DO 20 I=1,8
c          IV=IPE(K,I)
          IV=node((k-1)*nnode(1)+i)
          IF(IV.EQ.0) GO TO 20
          DO 10 J=1,3
c10        XYZ(J,I)=COP(IV,J)
10        XYZ(J,I)=COOR((IV-1)*3+J)
c          YY=YY+COP(IV,2)         ! Y轴为竖向，向上为正
          YY=YY+COOR((IV-1)*3+2)   ! Y轴为竖向，向上为正
20      CONTINUE
      
        III=IJKL(K)
        YY=YY/FLOAT(III)           ! 单元平均高度
        YL=HDMAX                   ! YL最大坝高
        DO 30 KK=1,NJJ             ! 对各级循环
          Y0=HDAM(KK)              ! 各级坝高
          Y=Y0-YY                  ! 坝高程和单元高程差
          IF(Y .LE. 0.0) GO TO 30  ! 单元就在本层单元的最顶层
          IF(Y .LE. YL) YL=Y
30      CONTINUE
c        V=VML(K)
        V=VOL(K)
        IF(ABS(V).LE.1.0E-5) GO TO 40  ! 判断为实体元，接触元在40处理
        P=0.65*GAM*YL !K0=0.8 
        G0=AE(2,NM)
        CN=AE(3,NM)
        RF=AE(4,NM)
        C=AE(5,NM)
        FI=AE(6,NM)
        DELFI=AE(11,NM)
c        FI=FI-DELFI*ALOG10(P/10.0)
        FI=FI-DELFI*DLOG10(P/10.0)
        FI=FI*3.1415926/180.0
        QF=(2.0*C*COS(FI)+2.0*P*SIN(FI))/(1-SIN(FI))
        S=(GAM*YL-P)/QF                ! 初始应力水平
        IF(S.GE.1.0) S=1.0
        SLEVEL=S
        IF(S.GE.0.95) SLEVEL=0.95
        PE=P/10.0
        IF(PE.LE.1.0) PE=1.0
        E0=G0*10.0*(PE)**CN
        ETL(K)=E0*(1-RF*S)*(1-RF*S)
        IF(ETL(K).LT.10.0) ETL(K)=10.0
        BT=G*10.0*(PE)**F 
        UTL(K)=0.5-ETL(K)/6.0/BT 
        IF(UTL(K).LT.PUMIN) UTL(K)=PUMIN
        IF(UTL(K).GT.PUMAX) UTL(K)=PUMAX
        GAMTL(K)=GAM
        GO TO 100
40      IF(ABS(GAM-0.1).GT.1.0E-5) GO TO 50
        ! 以下为摩擦单元
        GAM=AE(9,NM)
        IF(GAM.GT.1.0) GO TO 45 
        ! 面板下的接触单元
        ALFA=AN/SQRT(1.0+AN*AN)
        P=0.5*2.4*ALFA
        GO TO 46
45      P=0.36*GAM*HDZ
46      G0=AE(2,NM)
        CN=AE(3,NM)
        RF=AE(4,NM)
        ETL(K)=G0*1.0*(P/10.0)**CN
        UTL(K)=G0*1.0*(P/10.0)**CN
        GAMTL(K)=AE(10,NM)
        IF(GAM.GT.1.0) GO TO 100
        FI=AE(6,NM)*3.1415926/180.0
        S=1.0/AN/TAN(FI)
        IF(S.GT.1.0)  S=1.0
        ETL(K)=ETL(K)*(1.0-RF*S)*(1.0-RF*S)
        UTL(K)=UTL(K)*(1.0-RF*S)*(1.0-RF*S)
        GO TO 100
        ! 以下为实体缝单元
50      X=AE(2,NM)
        Y=AE(3,NM)
        ETL(K)=22.50*X !SHEAR1
        UTL(K)=60.80*X+140.00*Y    ! SHEAR2
        GAMTL(K)=65.00*X+53.00*Y   ! TENSION AND COMPRESS
        GO TO 100
80      ETL(K)=AE(2,NM)
        UTL(K)=AE(3,NM)
        GAMTL(K)=GAM
100   CONTINUE

      do i=1, num
        ET(iE_lg(i)) = ETL(i)
        UT(iE_lg(i)) = UTL(i)
        GAMT(iE_lg(i)) = GAMTL(i)
      enddo
      
      if( mype .eq. 0) then
        do iprc=1, npes-1
          call My_recvai(mype,iprc,ndata,1)
          allocate(rdata(ndata),idata(ndata))
          call My_recvai(mype,iprc,idata,ndata)
          call My_recvar(mype,iprc,rdata,ndata)
          do i=1, ndata
            ET(idata(i)) = rdata(i)
          enddo
          call My_recvar(mype,iprc,rdata,ndata)
          do i=1, ndata
            UT(idata(i)) = rdata(i)
          enddo
          call My_recvar(mype,iprc,rdata,ndata)
          do i=1, ndata
            GAMT(idata(i)) = rdata(i)
          enddo
          deallocate(idata,rdata)
        enddo

        if( mype.eq.0 ) WRITE(16,110)
110     FORMAT('              K,EI,UI,GAM=')
        if( mype.eq.0 ) THEN
          DO 60 K=1,NE
            WRITE(30,120) K,ET(K),UT(K),GAMT(K)
60        CONTINUE
        ENDIF
120     FORMAT(2(1X,I6,1X,F11.1,F10.2,F11.1))
      else
        call My_sendai(0,mype,num,1)
        call My_sendai(0,mype,iE_lg,num)
        call My_sendar(0,mype,etl,num)
        call My_sendar(0,mype,utl,num)
        call My_sendar(0,mype,gamtl,num)
      endif

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!   SUBRUTINE WHDUS   !!!!!!!!!!!!!!!!!
!     FUNCTION        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE WHDUS(K,III,SS,UU)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8),U(3)
      DIMENSION BAL(3,3),SS(3),UVW(8,3),UU(3)
      COMMON /CMN51/FUN(8),P(2,8),XJR(2,3)
      COMMON /A2/IPE(MAXN,8)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /A5/R(MAXN1)
      COMMON /A4/N,NH,MX,JR(3,MAXN)
      COMMON /A13/ET(MAXN),UT(MAXN),GAMT(MAXN)

      DO 5 I=1,3
        UU(I)=0.0
        SS(I)=0.0
5     CONTINUE
      DO 10 I=1,III
c        IL=IPE(K,I)
        IL=node((K-1)*nnode(1)+I)
        DO 10 J=1,3
c          XYZ(J,I)=COP(IL,J)
          XYZ(J,I)=COOR((IL-1)*3+J)
          IV=JRL(J,IL)
          IF(IV.EQ.0) UVW(I,J)=0.0
          IF(IV.NE.0) UVW(I,J)=R(IV)
10    CONTINUE
      IF(III.EQ.8) CALL WHD4FX(XYZ,BAL)
      IF(III.EQ.6) CALL CL3(XYZ,BAL)
      IM=III/2
      DO 80 J=1,3
        U(J)=0.0
        DO 80 I=1,IM
          IC=I+IM
          U(J)=U(J)+UVW(I,J)-UVW(IC,J) !compress IS POSITIVE
80    CONTINUE
      DO 90 I=1,3
        U(I)=U(I)/FLOAT(IM)
90    CONTINUE
      DO 95 I=1,3
        UU(I)=0.0
        DO 95 J=1,3
          UU(I)=UU(I)+BAL(I,J)*U(J)
95    CONTINUE
      SS(1)=ETL(K)*UU(1)
      SS(2)=UTL(K)*UU(2)
      SS(3)=GAMTL(K)*UU(3)

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!   SUBRUTINE    PK1   !!!!!!!!!!!!!!!!
!     FUNCTION      摩擦单元及缝单元的弹性常数 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PK1(K,III,FF,DU)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      DIMENSION FF(3),DU(3),DU2(3),PK(3),F1(3)
      COMMON /NJX/NJJ,MBJ,HDMAX,hdam(100),hhy(100),ii0(100),HHY2(100)
      COMMON /STRM/SSS(MAXN,3)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /A8/ME(MAXN)
      COMMON /A13/ET(MAXN),UT(MAXN),GAMT(MAXN)
      COMMON /A16/DZ(MAXN,3)

      CALL WHDUS(K,III,FF,DU)
      DO 20 I=1,3
        F1(I)=SSS(K,I)+FF(I)
20    CONTINUE
c      MN=ME(K)
      MN=node(K*nnode(1))
      T=AE(1,MN)
      IF(ABS(T-0.1).GT.1.0E-5) GO TO 60 
      ! 以下对摩擦元处理
      IF(F1(3).LT.0.01) GOTO 140
      G0=AE(2,MN)
      CN=AE(3,MN)
      RF=AE(4,MN)
      FI=AE(6,MN)*3.1415926/180.0
      PK(3)=AE(10,MN)
      E0=G0*1.0*(F1(3)/10.0)**CN
      S=F1(3)*TAN(FI)
      S1=ABS(F1(1))/S
      S2=ABS(F1(2))/S
      IF(S1.GT.1.0) S1=1.0
      IF(S2.GT.1.0) S2=1.0
      T1=1.0-RF*S1
      T2=1.0-RF*S2
      PK(1)=E0*T1*T1
      PK(2)=E0*T2*T2
      DO 34 I=1,2
        IF(PK(I).LT.10.0) PK(I)=10.0
34    CONTINUE
      IF(PK(1).GT.G0) PK(1)=G0
      IF(PK(2).GT.G0) PK(2)=G0
      GO TO 120
140   CONTINUE
      PK(1)=ETL(K)
      PK(2)=utl(K)
      PK(3)=GAMTL(K)
      GO TO 120
60    X=AE(2,MN)  ! 以下为对缝元的处理 DU(1)沉陷 DU(2)剪切 DU(3)拉伸和压缩
      Y=AE(3,MN)
      DO 65 I=1,3
        DU2(I)=DZ(K,I)+DU(I)
65    CONTINUE
      IF(ABS(X).LT.1.0E-5) GO TO 90
      IF(DU2(3).GT.0.0) GO TO 70
      A=17.5*0.5
      B=47.6*0.5
      W=1.0-B*ABS(DU2(3))
      IF(W.LE.0.0)W=0.1
      PK(3)=A*X/W/W
      GO TO 80
70    A=65.0*0.5
      B=41.0*0.5
      W=1.0-B*DU2(3)
      IF(W.LE.0.0)W=0.1
      PK(3)=A*X/W/W
80    CONTINUE
      IF(DU2(3).GT.0.002) PK(3)=1000000.0
      IF(DU2(3).LE.-0.019) PK(3)=200.0
      A=22.5*0.5
      B=40.0*0.5
      W=1.0-B*ABS(DU2(1))
      IF(W.LE.0.0)W=0.1
      PK(1)=A*X/W/W
      U=ABS(DU2(2))
      IF(U.LT.0.0125) PK(2)=60.80*X*0.5
      IF(U.GE.0.0125) PK(2)=56.00*X*0.5
90    IF(ABS(Y).LT.1.0E-5) GO TO 115
      IF(DU2(3).GT.0.0) GO TO 100
      U=ABS(DU2(3))
      IF(U.LE.0.0115) PK(3)=400.00*Y*0.5+PK(3)
      IF(U.GT.0.0115) PK(3)=60.00*Y*0.5+PK(3)
      GO TO 110
100   U=ABS(DU2(3))
      IF(U.LE.0.0115) PK(3)=53.00*Y*0.5+PK(3)
      IF(U.GT.0.0115) PK(3)=19.60*Y*0.5+PK(3)
110   PK(1)=0.0+PK(1)
      PK(2)=140.00*Y*0.5+PK(2)
      IF(PK(3).LE.20000)PK(3)=20000
115   CONTINUE
120   CONTINUE
      ETL(K)=PK(1)
      UTL(K)=PK(2)
      GAMTL(K)=PK(3)

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!  SUBRUTINE EVPK   !!!!!!!!!!!!!!!!!!!!!!
!     FUNCTION              各级计算单元的非线性弹性常数        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE EVPK(IFO,KK)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8),FF(3),S(MAXN),FI(MAXN),DU1(3)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /STRM/SSS(MAXN,3)
      COMMON /A15/STRZ(MAXN,6),EPGZ(MAXN,6)
      COMMON /A2/IPE(MAXN,8)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /NJX/NJJ,MBJ,HDMAX,hdam(100),hhy(100),ii0(100),HHY2(100)
      COMMON /A13/ET(MAXN),UT(MAXN),GAM(MAXN)
      COMMON /A16/DZ(MAXN,3)
      COMMON /A8/ME(MAXN)
      COMMON /JLMX/MX
      COMMON /EJOINT/IJKL(MAXN),VML(MAXN)

      NUM = mnode(1)
      if( .not. allocated(DZL) ) allocate(DZL(num,3))

      KE=II0(KK)
c      DO 20 K=1,KE
      DO 20 K=1,NUM
        if( iE_lg(K) > KE ) cycle
        DO 30 I=1,8
c          IV=IPE(K,I)
          IV=node((K-1)*nnode(1)+I)
          IF(IV.EQ.0) GO TO 30
c          XYZ(1:3,I)=COP(IV,1:3)
          do j=1, 3
            XYZ(j,I)=COOR((IV-1)*3+j)
          enddo
30      CONTINUE
        III=IJKL(K)
c        V1=VML(K)
c        IM=ME(K)
        V1=VOL(K)
        IM=node(K*nnode(1))
        IF(AE(1,IM).GT.0.3) CALL EBMOD(K,IFO,S1,FI1)
!        S(K)=S1
!        FI(K)=FI1 !FRICTION AGLE
        IF(AE(1,IM).GT.0.2) GO TO 20
        CALL PK1(K,III,FF,DU1)
        IF(IFO.NE.10) GOTO 45
        DO 40 I=1,3
          IF(KK.GE.MBJ) SSSL(K,I)=FF(I)+SSSL(K,I)
          IF(KK.LT.MBJ.AND.FF(3).GT.1.0) SSSL(K,I)=FF(I)+SSSL(K,I)
!          IF(I.EQ.3)THEN
!          IF(DU1(I).GT.0)DU1(I)=0.0
!          ENDIF
          DZL(K,I)=DZL(K,I)+DU1(I) ! 竣工期没有缝位移
          JI=I+3
c          STRZ(K,JI)=DZL(K,I)*100 ! UNIT OF DISPLACEMENT OF JOINT IS CM
          STRZL(K,JI)=DZL(K,I)*100.D0 ! UNIT OF DISPLACEMENT OF JOINT IS CM
40      CONTINUE
45      CONTINUE
20    CONTINUE

150   FORMAT(' K,ET,VT,GAM,S=')
160   FORMAT(2(I5,F10.1,F9.2,F10.1,F5.2,i3))

110   RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!  SUBRUTINE  KRGAM   !!!!!!!!!!!!!!!!!
!     FUNCTION        形成总劲度阵及自重引起的荷载向量
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE  KRGAM(KK)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8)
      COMMON /A1/SK(125930000),SM(MAXN1)  ! 刚度矩阵
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /A6/MA(MAXN1)
      COMMON /A2/IPE(MAXN,8)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A4/N,NH,MX,JR(3,MAXN)
      COMMON /A5/R(MAXN1)
      COMMON /CMN4/K,NET,NRW,NDE
      COMMON /GAUSS/RSTG(3),H(3)
      COMMON /CMN3/SKE(24,24),RF(24),CE(24)
      COMMON /NJX/NJJ,MBJ,HDMAX,HD(100),HY(100),II0(100),HHY2(100)
      COMMON /A13/ET(MAXN),UT(MAXN),GAMT(MAXN)
      COMMON /A8/ME(MAXN)
      COMMON /EJOINT/IJKL(MAXN),VML(MAXN)
      COMMON/STRU3/ICOL(125930000)

      IF(KK.GT.1) JK=II0(KK-1)   ! 上步加载单元
      IF(KK.EQ.1) JK=0

      DO 10 I=1,NH
        SK(I)=0.0
10    CONTINUE

      DO 15 I=1,N
        SM(I)=0.0
15    CONTINUE

      ! 初始化刚度矩阵
      do i=1, nnz
        am(i) = 0.D0
      enddo

      KE=II0(KK)                ! 本步加载单元数
      if( mype.eq.0 ) write(*,*) '本步加载单元数：', KE
      NUM = mnode(1)
c      DO 20 K=1,KE
      DO 20 K=1,NUM
        ne_g = iE_lg(k)
        if( ne_g > KE ) cycle
c        NM=ME(K)
        NM = node(K*nnode(1))   ! 单元材料编号
        AIX=AE(1,NM)            ! 取材料参数
        ! 取单元节点坐标
        DO 30 I=1,8
c          IV=IPE(K,I)
          IV=node((k-1)*nnode(1)+I)  ! 节点编号
          if( iN_lg(IV) .eq. 29 ) print *,'29----,mype,K =', mype,K
          IF(IV.EQ.0) GO TO 30
          DO 28 J=1,3
c            XYZ(J,I)=COP(IV,J)
            XYZ(J,I)=COOR((IV-1)*3+J)       ! 节点坐标
28        CONTINUE
30      CONTINUE
        III=IJKL(K)         ! 单元节点个数?
c        V1=VML(K)
c        E=ET(K)
c        U=UT(K)
c        GAM=GAMT(K)
        V1=VOL(K)            ! 体积
        E=ETL(K)
        U=UTL(K)
        GAM=GAMTL(K)
        if( NE_G.eq.1 .and. mype.eq.1) then
          write(*,*) 'In KRGAM,1--K,V1,E,U,GAM=',K,V1,E,U,GAM
        endif

        IF(III.EQ.4) CALL KRS(K,JK,GAM)
!     FUNCTION  三棱柱单元的劲度阵
        IF(V1.LE.1.0E-5 .AND. III.EQ.6) CALL SJXK(E,U,GAM,XYZ)
!     FUNCTION  摩擦单元劲度阵
        IF(V1.LE.1.0E-5 .AND. III.EQ.8.AND.ABS(AIX-0.1).LT.1.0E-5)
     *     CALL WHDMK(E,U,GAM,XYZ)
!     FUNCTION  缝单元劲度阵
        IF(V1.LE.1.0E-5 .AND. III.EQ.8.AND.ABS(AIX-0.2).LT.1.0E-5)
     *     CALL WHDFK(E,U,GAM,XYZ)
!     FUNCTION  (8结点)单元劲度阵,及自重荷载列阵  
        IF(V1.GT.1.0E-5 .AND. III.EQ.8) 
     +     CALL STIF(K,JK,E,U,GAM,XYZ,III,ne_g)
!     FUNCTION  6结点单元劲度矩阵
        IF(V1.GT.1.0E-5 .AND. III.EQ.6)
     +     CALL STIF6(K,JK,E,U,GAM,XYZ,III)

c        if( ne_g .eq. 1 ) then
c          print *,'---------------------------------------'
c          print *,'iE_lg:', iE_lg(k),',k,v1,III=',k,v1,III
c          do i=1,24
c            write(*,*) (ske(i,j),j=1,24)
c          ENDDO
c          print *,'--------RF--------------------',mype,K,ne_g
c          write(*,'(24E12.4)') (RF(i),i=1,24)
c          print *,'---------------------------------------'
c        endif

        ! 组装总刚和右端项
        DO 45 II=1,8               !单刚中的行节点
c          IV=IPE(K,II)
          IV=node((K-1)*nnode(1)+II)
          IVG = iN_lg(IV)
          IF(IV.EQ.0) GO TO 45     !
          IQ1=(II-1)*3
          DO 40 IJ=1,3
c            IU=JR(IJ,IV)          !自由度编号
            IU=JRL(IJ,IV)          !局部自由度编号,只整理内部节点上的自由度
            IP1=IQ1+IJ             !本自由度对应的单刚行号
            IF(IU.EQ.0) GOTO 40
c            JV=MA(IU)
            JV = NA(IU)            !该行非零元在总刚中的起始位置
            DO 55 JI=1,8           !单刚中的列节点
c              IW=IPE(K,JI)
              IW=node((K-1)*nnode(1)+JI)  !列节点编号
              IF(IW.EQ.0) GO TO 55
              IQ2=(JI-1)*3
              DO 50 JJ=1,3         !列节点自由度循环
                IP2=IQ2+JJ         !列节点自由度对应单刚中的列号
c                IX=JR(JJ,IW)
                IX=JRG(JJ,IW)      !列节点自由度对应的整体列号(方程号)
                IF(IX.EQ.0)  GOTO 50
c                IF(IX.GT.IU) GOTO 50      !只存放下三角
c                if(iu.gt.1) il=MA(iu)-MA(iu-1)
c                if(iu.eq.1) il=MA(1)
                il = na(iu+1) - na(iu)   !本行非零元个数
                do 25 k3=1,il            !寻找单刚中元素对应在总刚中的位置
c                  if(iu.gt.1) k1=MA(iu-1)+k3
c                  if(iu.eq.1) k1=k3
                  k1 = JV+K3             !非零元在总刚中的位置
c                  k2=ICOL(k1)
                  k2=ia(k1)              !非零元对应的列号
                  if(k2.eq.ix) then      !累加单刚到总刚
c                    sk(k1)=sk(k1)+ske(ip1,ip2)
                     am(k1) = am(k1)+ske(ip1,ip2)
                    goto 26
                  endif
25              continue
26              continue
50            CONTINUE
55          CONTINUE
60          CONTINUE
            IF(NET.EQ.1) R(IU)=R(IU)+RF(IP1)  !Gravity
c            if( ne_g .eq. 333 ) then
c              write(*,'(a,12i5)') '---------------IV,IVG,K',IV,IVG,K
c              write(*,'(a,2i5,2E12.4)') 'IU,IP1,R(IU),RF(IP1)=',IU,IP1,R(IU),RF(IP1)
c              print *,'---------------------------------------'
c            ENDIF
40        CONTINUE
45      CONTINUE

20    CONTINUE

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!  SUBRUTINE WHDMK  !!!!!!!!!!!!!!!!
!     FUNCTION         摩擦单元劲度阵
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE WHDMK(E,U,GAM,XYZ)
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8),PK(3),BAK(3,3),BAL(3,3),BAD(3,3),A(8,8)
      COMMON /CMN3/SKE(24,24),RF(24),CE(24)

      PK(1)=E
      PK(2)=U
      PK(3)=GAM
      A1=XYZ(1,1)-XYZ(1,2)
      A2=XYZ(2,1)-XYZ(2,2)
      A3=XYZ(3,1)-XYZ(3,2)
      B=SQRT(A1*A1+A2*A2+A3*A3)
      A1=XYZ(1,1)-XYZ(1,3)
      A2=XYZ(2,1)-XYZ(2,3)
      A3=XYZ(3,1)-XYZ(3,3)
      AL=SQRT(A1*A1+A2*A2+A3*A3)
      S=B*AL
      DO 10 I=1,4
        A(I,I)=4.0
10    CONTINUE
      A(2,1)=2.0
      A(3,1)=2.0
      A(3,2)=1.0
      A(4,1)=1.0
      A(4,2)=2.0
      A(4,3)=2.0
      DO 20 I=1,3
      DO 20 J=I+1,4
        A(I,J)=A(J,I)
20    CONTINUE
C     IF(NEE.EQ.49) WRITE(16,150) B,AL,S,((A(I,J),J=1,4),I=1,4)
150   FORMAT(' B,L,S=',3F15.3/'  A(4,4)='/4F10.2/4F10.2/4F10.2/4F10.2)
      DO 30 I=1,4
      DO 30 J=1,4
        A(I,J)=A(I,J)*S/36.0
30    CONTINUE
      DO 40 I=5,8
        II=I-4
        DO 40 J=5,8
          JJ=J-4
          A(I,J)=A(II,JJ)
40    CONTINUE
      DO 50 I=1,4
      DO 50 J=5,8
        JJ=J-4
        A(I,J)=-A(I,JJ)
50    CONTINUE
      DO 60 I=5,8
        II=I-4
        DO 60 J=1,4
          A(I,J)=-A(II,J)
60    CONTINUE
      DO 70 I=1,24
        RF(I)=0.0
        CE(I)=0.0
        DO 70 J=1,24
          SKE(I,J)=0.0
70    CONTINUE
      DO 80 I=1,3
      DO 80 J=1,3
        BAK(I,J)=0.0
        IF(I.EQ.J) BAK(I,J)=PK(I)
80    CONTINUE
      CALL WHD4FX(XYZ,BAL)
      CALL CH29(BAL,BAK,BAD)   !!!!!!!!!!!!!!!!   !!!!!!!!!!!!!!!!
      
      DO 110 I=1,8
      DO 110 J=1,8
        DO 100 K=1,3
          M=(I-1)*3+K
          DO 90 L=1,3
            N=(J-1)*3+L
            SKE(M,N)=A(I,J)*BAD(K,L)
90        CONTINUE
100     CONTINUE
110   CONTINUE

      RETURN
      END
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!  SUBRUTINE WHDFK  !!!!!!!!!!!!!!
!     FUNCTION       缝单元劲度阵
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE WHDFK(E,U,GAM,XYZ)
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8),PK(3),BAK(3,3),BAL(3,3),BAD(3,3),A(8,8)
      COMMON /CMN3/SKE(24,24),RF(24),CE(24)
      PK(1)=E
      PK(2)=U
      PK(3)=GAM
      A1=XYZ(1,1)-XYZ(1,2)
      A2=XYZ(2,1)-XYZ(2,2)
      A3=XYZ(3,1)-XYZ(3,2)
      B=SQRT(A1*A1+A2*A2+A3*A3)
      A1=XYZ(1,1)-XYZ(1,3)
      A2=XYZ(2,1)-XYZ(2,3)
      A3=XYZ(3,1)-XYZ(3,3)
      AL=SQRT(A1*A1+A2*A2+A3*A3)
      IF(B.GT.AL) AL=B
      DO 10 I=1,4
      DO 10 J=1,4
      A(I,J)=0.25*AL
10      CONTINUE
      DO 40 I=5,8
      II=I-4
      DO 40 J=5,8
      JJ=J-4
      A(I,J)=A(II,JJ)
40      CONTINUE
      DO 50 I=1,4
      DO 50 J=5,8
      JJ=J-4
      A(I,J)=-A(I,JJ)
50      CONTINUE
      DO 60 I=5,8
      II=I-4
      DO 60 J=1,4
      A(I,J)=-A(II,J)
60      CONTINUE
      DO 70 I=1,24
      RF(I)=0.0
      CE(I)=0.0
      DO 70 J=1,24
      SKE(I,J)=0.0
70      CONTINUE
      DO 80 I=1,3
      DO 80 J=1,3
      BAK(I,J)=0.0
      IF(I.EQ.J) BAK(I,J)=PK(I)
80      CONTINUE
      CALL WHD4FX(XYZ,BAL)
      CALL CH29(BAL,BAK,BAD)
      DO 110 I=1,8
      DO 110 J=1,8
      DO 100 K=1,3
      M=(I-1)*3+K
      DO 90 L=1,3
      N=(J-1)*3+L
      SKE(M,N)=A(I,J)*BAD(K,L)
90    CONTINUE
100   CONTINUE
110   CONTINUE
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!! SUBRUTINE KRS   !!!!!!!!!!!!!!!!!!!!
!    FUNCTION           计算 Krs 子块
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE KRS(K,JK,GAM)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      COMMON /CMN4/NEE,NET,NRW,NDE
      COMMON /CMN3/SKE(24,24),RF(24),CE(24)
      COMMON /A13/ET(MAXN),UT(MAXN),GAMT(MAXN)
      COMMON/B11/BI(4),CI(4),DI(4)

      NE_G = iE_lg(K)    ! 单元整体编号
      DO 5 I=1,24
        RF(I)=0.0
        CE(I)=0.0
        DO 5 J=1,24
          SKE(I,J)=0.0
5     CONTINUE

      CALL FOUR(K,V)
      DO 30 L=1,4
        LQ=(L-1)*3
        DO 30 M=1,4
          MQ=(M-1)*3
          BR=BI(L)
          BS=BI(M)
          CR=CI(L)
          CS=CI(M)
          DR=DI(L)
          DS=DI(M)
          ITT=(-1)**(L+M)
c          V1=UT(K)
          V1=UTL(K)
          A1=V1/(1.0-V1)
          A2=(1.0-2.0*V1)/(1.0-V1)/2.0
c          T1=(1.0-V1)*ET(K)/(1.0-2.0*V1)/V/36.0/(1.0+V1)
          T1=(1.0-V1)*ETL(K)/(1.0-2.0*V1)/V/36.0/(1.0+V1)
          SKE(LQ+1,MQ+1)=ITT*T1*(BR*BS+A2*(CR*CS+DR*DS))
          SKE(LQ+1,MQ+2)=ITT*T1*(A1*BR*CS+A2*CR*BS)
          SKE(LQ+1,MQ+3)=ITT*T1*(A1*BR*DS+A2*DR*BS)
          SKE(LQ+2,MQ+1)=ITT*T1*(A1*CR*BS+A2*BR*CS)
          SKE(LQ+2,MQ+2)=ITT*T1*(CR*CS+A2*(BR*BS+DR*DS))
          SKE(LQ+2,MQ+3)=ITT*T1*(A1*CR*DS+A2*DR*CS)
          SKE(LQ+3,MQ+1)=ITT*T1*(A1*DR*BS+A2*BR*DS)
          SKE(LQ+3,MQ+2)=ITT*T1*(A1*DR*CS+A2*CR*DS)
          SKE(LQ+3,MQ+3)=ITT*T1*(DR*DS+A2*(BR*BS+CR*CS))
30    CONTINUE
      W=V*GAM
      DO 35 I=1,4
        LJ=(I-1)*3
        DO 35 J=1,3
          IJ=LJ+J
          CE(IJ)=0.25*W/9.80665
35    CONTINUE
      IF(NET.NE.1) GO TO 60
c      IF(K.LE.JK) GO TO 60
      IF(NE_G.LE.JK) GO TO 60
      DO 40 I=1,4
        J=I*3-1
        RF(J)=0.25*W
40    CONTINUE
60    CONTINUE

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE FOUR    !!!!!!!!!!!!
!     FUNCTION  !已改 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FOUR(K,V)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
C     IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(4),Y(4),Z(4)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /A2/IPE(MAXN,8)
      COMMON/B11/BI(4),CI(4),DI(4)
      COMMON /EJOINT/IJKL(MAXN),VML(MAXN)

      DO 10 I=1,4
c        KP=IPE(K,I)
c        X(I)=COP(KP,1)
c        Y(I)=COP(KP,2)
c        Z(I)=COP(KP,3)
        KP=node((K-1)*nnode(1)+I)
        X(I)=COOR((KP-1)*3+1)
        Y(I)=COOR((KP-1)*3+2)
        Z(I)=COOR((KP-1)*3+3)
10    CONTINUE
      DO 20 J=1,4
c        IF(J+1-4) 30, 30, 40
        if(j+1-4 <= 0) goto 30
        if(j+1-4 > 0) goto 40
30      L=J+1
        GO TO 50
40      L=1
50      CONTINUE
c        IF(J+2-4)60, 60, 70
        if(j+2-4 <= 0) goto 60
        if(j+2-4 > 0) goto 70
60      M=J+2
        GO TO 80
70      M=J-2
80      CONTINUE
c        IF(J+3-4)90, 90, 100
        if(j+3-4 <= 0) goto 90
        if(j+3-4 > 0) goto 100
90      I=J+3
        GO TO 200
100     I=J-1
200     BI(J)=Y(L)*(Z(I)-Z(M))+Y(M)*(Z(L)-Z(I))+Y(I)*(Z(M)-Z(L))
        CI(J)=X(L)*(Z(M)-Z(I))+X(M)*(Z(I)-Z(L))+X(I)*(Z(L)-Z(M))
        DI(J)=X(L)*(Y(I)-Y(M))+X(M)*(Y(L)-Y(I))+X(I)*(Y(M)-Y(L))
20    CONTINUE
c      V=VML(K)
      V=VOL(K)

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE STIF    !!!!!!!!!!!!!
!     FUNCTION        单元劲度阵(8结点),及自重荷载列阵  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE STIF(KA,JK,E,U,GAM,XYZ,III,NE_G)
      implicit real*8(a-h,o-z)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,8),RJX(3,3),Q(3,8),BV(24),D(9)
      COMMON /CMN5/FUN(8),P(3,8),XJR(3,3)
      COMMON /CMN4/NEE,NET,NRW,NDE
      COMMON /GAUSS/RSTG(3),H(3)            ! 高斯积分点坐标和权
      COMMON /CMN3/SKE(24,24),RF(24),CE(24)

      D1=E*(1.0-U)/((1.0+U)*(1.0-2.0*U))
      D2=E*U/((1.0+U)*(1.0-2.0*U))
      D3=E*0.5/(1.0+U)

      DO 5 I=1,24
        RF(I)=0.0
        CE(I)=0.0
        DO 5 J=1,24
          SKE(I,J)=0.0
5     CONTINUE

      DO 100 IT=1,3
        T=RSTG(IT)
        TH=H(IT)
        DO 90 IS=1,3
          S=RSTG(IS)
          SH=H(IS)
          DO 80 IR=1,3
            R=RSTG(IR)
            RH=H(IR)
!     FUNCTION--Jacobi矩阵的逆阵
            CALL RMSD(XYZ,Q,DET,R,S,T,RJX,III)
            DO 35 I=1,8
              LJ=(I-1)*3
              ALS=-FUN(I)*RH*SH*TH*DET*GAM
              DO 35 J=1,3
                IJ=LJ+J
                CE(IJ)=CE(IJ)+ALS/9.80665
35          CONTINUE
20          IF(NET.EQ.0) GO TO 40
            IF(NE_G.LE.JK) GO TO 40
            DO 30 I=1,8
              J=I*3-1
              RF(J)=RF(J)-FUN(I)*RH*SH*TH*DET*GAM
30          CONTINUE
c            if(NE_G.eq.106) THEN
c              write(*,*) 'In STIF,IT,IS,IR=',IT,IS,IR
c              write(*,*) '------------RF-----------',KA,NE_G
c              write(*,'(24E12.4)') (RF(j),j=1,24)
c            ENDIF
40          K3=0
            DO 10 I=1,8
              K3=K3+3
              K2=K3-1
              K1=K2-1
              BV(K1)=Q(1,I)
              BV(K2)=Q(2,I)
              BV(K3)=Q(3,I)
10          CONTINUE
            DO 70 I=1,24
              DO 60 J=I,24
                SKE(I,J)=SKE(I,J)+BV(I)*BV(J)*DET*RH*SH*TH
60            CONTINUE
70          CONTINUE
80        CONTINUE
90      CONTINUE
100   CONTINUE
      DO 110 I=2,24
        M=I-1
        DO 110 J=1,M
          SKE(I,J)=SKE(J,I)
110   CONTINUE
      KK=-3
      DO 130 I=1,8
        KK=KK+3
        K1=KK+1
        K2=KK+2
        K3=KK+3
        DO 125 J=1,8
          LL=(J-1)*3
          L1=LL+1
          L2=L1+1
          L3=L2+1
          IC=0
          DO 120 K=1,3
            M=KK+K
            DO 115 L=1,3
              NL=L+LL
              IC=IC+1
              D(IC)=SKE(M,NL)
115         CONTINUE
120       CONTINUE
          SKE(K1,L1)=D(1)*D1+(D(5)+D(9))*D3
          SKE(K2,L2)=D(5)*D1+(D(1)+D(9))*D3
          SKE(K3,L3)=D(9)*D1+(D(5)+D(1))*D3
          SKE(K1,L2)=D(2)*D2+D(4)*D3
          SKE(K2,L1)=D(4)*D2+D(2)*D3
          SKE(K2,L3)=D(6)*D2+D(8)*D3
          SKE(K3,L2)=D(8)*D2+D(6)*D3
          SKE(K3,L1)=D(7)*D2+D(3)*D3
          SKE(K1,L3)=D(3)*D2+D(7)*D3
125     CONTINUE
130   CONTINUE

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE EBMODEL   !!!!!!!!!!!!
!     FUNCTION                 E-B模型 !!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE EBMOD(K,IFO,S1,AF) 
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      DIMENSION FF(3)
      COMMON /JLMX/MX
      COMMON /ST/STR(MAXN,6),EPG(MAXN,6)
      COMMON /STRM/SSS(MAXN,3)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /A2/IPE(MAXN,8)
      COMMON /A8/ME(MAXN)
      COMMON /A13/E(MAXN),U(MAXN),GAMT(MAXN)
      COMMON /SQ/S0(MAXN),Q0(MAXN)
      COMMON /POSONGBI/PUMIN,PUMAX
c      NM=ME(K)

      NM = node(K*nnode(1))

      G0=AE(2,NM)
      CN=AE(3,NM)
      RF=AE(4,NM)
      G=AE(7,NM)
      F=AE(8,NM)
      DF=AE(9,NM)

      DO 10 I=1,3
        FF(I)=SSSL(K,I)
10    CONTINUE

c      if( K.eq.67 .and. mype.eq.1 ) THEN
c        write(*,*) '000 In ebmod, 67--,G0,CN,RF,G,F,DF='
c        write(*,*) G0,CN,RF,G,F,DF
c        write(*,*) 'FF:',(FF(i),i=1,3)
c      ENDIF

      IF(G0.GE.10000.AND.G0.LE.200000) GOTO 300
!     IF(G0.GT.200000.AND.G0.LE.1000000) GOTO 200
      IF(G0.GT.200000) GOTO 200
      IF(FF(1).EQ.0.AND.FF(2).EQ.0.AND.FF(3).EQ.0)GOTO 300
      P=(FF(1)+FF(3))/2.0
      Q=ABS(FF(1)-FF(3)) 
      IF(P.LT.10.0) P=10.0
      C=AE(5,NM)
      FI=AE(6,NM)
      DELFI=AE(11,NM)
      PE=P/10.0
      IF(PE.LE.1.0)PE=1.0
c      FI=FI-DELFI*ALOG10(PE)
      FI=FI-DELFI*DLOG10(PE)
      AF=FI
      FI=FI*3.1415926/180.0
      IF(MX.EQ.6) GO TO 40
      AMU=1.0-2.0*(FF(2)-FF(3))/(FF(1)-FF(2))
      ODE=(-AMU/1.7320508)
      CTA=ATAN(ODE)
      CM=3.0*SIN(FI)*P+3.0*C*SIN(FI)
      CP=1.7320508*COS(CTA)+SIN(CTA)*SIN(FI)
      QF=CM/CP
      GO TO 60
40    QF=(2.0*C*COS(FI)+2.0*P*SIN(FI))/(1-SIN(FI)) 
60    S=Q/QF

c      if( k.eq.67) THEN
c        write(*,*) '111 IN ebmod, 67--- S =', S
c      ENDIF

      IF(S.GE.1.0)S=1.0
      IF(S.LE.0.001)S=0.001
!     ET0=E0*(1-RF*S)*(1-RF*S)            
      IF(S.GE.0.95*S0L(K)) THEN
        S1=S
        E0=G0*10.0*(PE)**CN
        ET=E0*(1-RF*S)*(1-RF*S)
        BT=G*10.0*(PE)**F
        VT=0.5-ET/(6.0*BT)
        IF(VT.LT.PUMIN) VT=PUMIN
        IF(VT.GT.PUMAX) VT=PUMAX
        GOTO 20
      ENDIF
      IF(S.LE.0.75*S0L(K)) GO TO 35
      S1=S
      E0=G0*10.0*(PE)**CN
      ET=E0*(1-RF*S)*(1-RF*S)
      EEUR=G0*20.0*(PE)**CN
      ET=ET+(S0L(k)-S)/(0.2*S0L(k))*(EEUR-ET)
      BT=G*10.0*(PE)**F
      VT=0.5-ET/(6.0*BT)
      IF(VT.LT.PUMIN) VT=PUMIN
      IF(VT.GT.PUMAX) VT=PUMAX
      GO TO 20
      
35    S1=S
      EUR=AE(10,NM)
      ET=EUR*10.0*(PE)**CN
      VT=0.35
      
20    CONTINUE
      
      IF(IFO.EQ.5) GO TO 50
      IF(S.GT.S0L(K))S0L(K)=S
      IF(Q.GT.Q0L(K))Q0L(K)=Q
50    CONTINUE
      ETL(K)=ET
      UTL(K)=VT
100   CONTINUE
      GO TO 300

200   IF(FF(3).LT.0.0) GO TO 210
      ETL(K)=G0
      IF(FF(1).GT.2500.0) ETL(K)=G0*0.5
      GO TO 300
210   ETL(K)=640000.0
      IF(FF(3).LT.-8.0) ETL(K)=560000.0
      IF(FF(3).LT.-22.0) ETL(K)=460000.0
      IF(FF(3).LT.-40.0) ETL(K)=230000.0
      IF(FF(3).LT.-60.0) ETL(K)=100000.0
      IF(FF(3).LT.-80.0) ETL(K)=50000.0
      IF(FF(3).LT.-100.0) ETL(K)=10000.0
      IF(FF(3).LT.-120.0) ETL(K)=4000.0
      GOTO 300
300   CONTINUE

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE SJXK    !!!!!!!!!!!!!!
!     FUNCTION    三棱柱单元的劲度阵
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SJXK(E,U,GAM,XYZ)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      DIMENSION A(6,6),XYZ(3,8),XYL(3),PK(3)
      DIMENSION BAL(3,3),BAD(3,3),BAK(3,3)
      COMMON /CMN3/SKE(24,24),RF(24),CE(24)

      PK(1)=E
      PK(2)=U
      PK(3)=GAM
      DO 20 I=1,3
        XL1=0.0
        II=I+1
        IF(II.EQ.4)II=1
        DO 10 J=1,3
          XL=XYZ(J,I)-XYZ(J,II)
          XL1=XL1+XL*XL
10      CONTINUE
        XYL(I)=SQRT(XL1)
20    CONTINUE
      A1=XYL(1)
      A2=XYL(2)
      A3=XYL(3)
      S1=(A1+A2+A3)/2.0
      S=SQRT(S1*(S1-A1)*(S1-A2)*(S1-A3))
      DO 30 I=1,3
        DO 30 J=1,3
          IF(I.EQ.J) A(I,J)=S/6.0
          IF(I.NE.J) A(I,J)=S/12.0
30    CONTINUE
      DO 40  I=4,6
        DO 40  J=4,6
          IF(I.EQ.J) A(I,J)=S/6.0
          IF(I.NE.J) A(I,J)=S/12.0
40    CONTINUE
      DO 50 I=1,3
        DO 50 J=4,6
          JJ=J-3
          IF(I.EQ.JJ) A(I,J)=-S/6.0
          IF(I.NE.JJ) A(I,J)=-S/12.0
50    CONTINUE
      DO 60 I=4,6
        DO 60 J=1,3
          JJ=I+3
          IF(I.EQ.JJ) A(I,J)=-S/6.0
          IF(I.NE.JJ) A(I,J)=-S/12.0
60    CONTINUE
      DO 70 I=1,24
        RF(I)=0.0
        CE(I)=0.0
        DO 70 J=1,24
          SKE(I,J)=0.0
70    CONTINUE
      CALL CL3(XYZ,BAL)
      DO 80 I=1,3
        DO 80 J=1,3
          BAK(I,J)=0.0
          IF(I.EQ.J) BAK(I,J)=PK(I)
80    CONTINUE
      CALL CH29(BAL,BAK,BAD)
      DO 110 I=1,6
        DO 110 J=1,6
          DO 100 K=1,3
            M=(I-1)*3+K
            DO 90 L=1,3
              N=(J-1)*3+L
              SKE(M,N)=A(I,J)*BAD(K,L)
90          CONTINUE
100       CONTINUE
110   CONTINUE

      RETURN
      END
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE CL3    !!!!!!!!!!!!!!
!     FUNCTION            三角形无厚度单元的局部坐标系
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CL3(XYZ,BAL)
c      USE ComData
c      USE FEMData
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8),BAL(3,3)

      B1=XYZ(1,2)-XYZ(1,1)
      B2=XYZ(2,2)-XYZ(2,1)
      B3=XYZ(3,2)-XYZ(3,1)
      B4=SQRT(B1*B1+B2*B2+B3*B3)
      BL=B1/B4
      BM=B2/B4
      BN=B3/B4
      B1=XYZ(1,3)-XYZ(1,1)
      B2=XYZ(2,3)-XYZ(2,1)
      B3=XYZ(3,3)-XYZ(3,1)
      B4=SQRT(B1*B1+B2*B2+B3*B3)
      EL=B1/B4
      EM=B2/B4
      EN=B3/B4
      B1=EN*BM-EM*BN
      B2=EL*BN-BL*EN
      B3=BL*EM-BM*EL
      B4=SQRT(B1*B1+B2*B2+B3*B3)
      BAL(3,1)=B1/B4
      BAL(3,2)=B2/B4
      BAL(3,3)=B3/B4
      XY=BAL(3,1)-BAL(3,2)
      B4=SQRT(BAL(3,3)*BAL(3,3)*2.0+XY*XY)
      BAL(2,1)=-BAL(3,3)/B4
      BAL(2,2)=BAL(3,3)/B4
      BAL(2,3)=XY/B4
      B1=BAL(2,2)*BAL(3,3)-BAL(3,2)*BAL(2,3)
      B2=BAL(3,1)*BAL(2,3)-BAL(2,1)*BAL(3,3)
      B3=-BAL(3,1)*BAL(2,2)+BAL(2,1)*BAL(3,2)
      B4=SQRT(B1*B1+B2*B2+B3*B3)
      BAL(1,1)=B1/B4
      BAL(1,2)=B2/B4
      BAL(1,3)=B3/B4

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE CH29    !!!!!!!!!!!!!!
!     FUNCTION        矩阵运算LKL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE  CH29(BAL,BAK,BAD)
      implicit real*8(a-h,o-z)
      DIMENSION BAL(3,3),BLA(3,3),BAD(3,3),BAK(3,3),BAC(3,3)

      DO 10 I=1,3
      DO 10 J=1,3
        BLA(I,J)=BAL(J,I)
10    CONTINUE

      DO 30 I=1,3
        DO 30 J=1,3
          BAC(I,J)=0.0
          DO 35 II=1,3
            BAC(I,J)=BAC(I,J)+BLA(I,II)*BAK(II,J)
35        CONTINUE
30    CONTINUE
      DO 40 I=1,3
      DO 40 J=1,3
        BAD(I,J)=0.0
        DO 42 II=1,3
          BAD(I,J)=BAD(I,J)+BAC(I,II)*BAL(II,J)
42      CONTINUE
40    CONTINUE
      RETURN
      END
      
      SUBROUTINE  CH28(BAL,SZ,BAD)
      implicit real*8(a-h,o-z)
      DIMENSION BAL(3,3),BAD(3),SZ(6)
 
      DO 40 I=1,3
      BAD(I)=0.0
      DO 42 II=1,3
      BAD(I)=BAD(I)+BAL(I,II)*SZ(II)
42    CONTINUE
40    CONTINUE
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE PDF    !!!!!!!!!!!!
!     FUNCTION           无厚度四边形单元形函数及其偏导数
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PFD(R,S,AAB,XYZ)
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8),XI(4),ETA(4),E(3)
      COMMON /CMN51/FUN(8),P(2,8),XJR(2,3)
      DATA XI/-1.0,1.0,-1.0,1.0/
      DATA ETA/-1.0,-1.0,1.0,1.0/
      DO 10 I=1,4
      FUN(I)=0.25*(1.0+XI(I)*R)*(1.0+ETA(I)*S)
10    CONTINUE
      DO 20 I=5,8
      J=I-4
      FUN(I)=-FUN(J)
20    CONTINUE
      DO 30 I=1,4
      P(1,I)=0.25*XI(I)*(1.0+ETA(I)*S)
      P(2,I)=0.25*ETA(I)*(1.0+XI(I)*R)
30    CONTINUE
      DO 40 I=5,8
      J=I-4
      P(1,I)=-P(1,J)
      P(2,I)=-P(2,J)
40    CONTINUE
      DO 70 I=1,2
      DO 60 J=1,3
      DE1=0.0
      DO 50 K=1,4
      DE1=DE1+P(I,K)*XYZ(J,K)
50    CONTINUE
      XJR(I,J)=DE1
60    CONTINUE
70    CONTINUE
      DO 80 I=1,3
80    E(I)=0.0
      DO 90 I=1,3
      E(1)=E(1)+XJR(1,I)*XJR(1,I)
      E(2)=E(2)+XJR(2,I)*XJR(2,I)
      E(3)=E(3)+XJR(1,I)*XJR(2,I)
90    CONTINUE
      DET=E(1)*E(2)-E(3)*E(3)
      AAB=SQRT(DET)
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE VOLM    !!!!!!!!!!!!!!
!     FUNCTION                    计算四面体的体积
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION  VOLM(X,Y,Z)
      implicit real*8(a-h,o-z)
C     IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(4),Y(4),Z(4)
      VOLM=((X(2)-X(1))*((Y(3)
     1  -Y(1))*(Z(4)-Z(1))-(Y(4)-Y(1))*(Z(3)-Z(1)))
     1  -(Y(2)-Y(1))*((X(3)-X(1))*(Z(4)-Z(1))-(X(4)-X(1))*(Z(3)-Z(1)))+
     1   (Z(2)-Z(1))*((X(3)-X(1))*(Y(4)-Y(1))-(X(4)-X(1))*(Y(3)-Y(1))))
     1 /6.0
C     WRITE(16,6) (X(IO),Y(IO),Z(IO),IO=1,4)
6     FORMAT(2X,'X=',3E12.5)
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE  V68    !!!!!!!!!!!!!!!
!     FUNCTION           (6或8)面体的体积
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE V68(XYZ,III,V)
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8),X(4),Y(4),Z(4),NOD(4,6),NOD1(4,3)
      DATA NOD/1,4,2,6,
     *     1,5,4,6,
     *     1,3,4,5,
     *     1,4,8,2,
     *     1,8,4,7,
     *     1,4,3,7/
      DATA NOD1/1,4,5,6,
     *      1,5,2,6,
     *      1,2,3,6/
      V=0.0
      IF(III.EQ.4) GO TO 60
      IF(III.EQ.6) GO TO 30
      DO 20  I=1,6
      DO 10 J=1,4
      K=NOD(J,I)
      X(J)=XYZ(1,K)
      Y(J)=XYZ(2,K)
      Z(J)=XYZ(3,K)
10      CONTINUE
      V=V+ABS(VOLM(X,Y,Z))
20      CONTINUE
      GO TO 100
30      DO 50 I=1,3
      DO 40 J=1,4
      K=NOD1(J,I)
      X(J)=XYZ(1,K)
      Y(J)=XYZ(2,K)
      Z(J)=XYZ(3,K)
40      CONTINUE
      V=V+ABS(VOLM(X,Y,Z))
50      CONTINUE
      GO TO 100
60      DO 70 I=1,4
      X(I)=XYZ(1,I)
      Y(I)=XYZ(2,I)
      Z(I)=XYZ(3,I)
70      CONTINUE
      V=ABS(VOLM(X,Y,Z))
100   RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE WHD4FX    !!!!!!!!!!!!!!
!     FUNCTION      无厚度单元的方向(8结点)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE WHD4FX(XYZ,BAL)
      implicit real*8(a-h,o-z)
C     IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XYZ(3,8),BAL(3,3)

      B1=XYZ(1,2)-XYZ(1,1)
      B2=XYZ(2,2)-XYZ(2,1)
      B3=XYZ(3,2)-XYZ(3,1)
      B4=SQRT(B1*B1+B2*B2+B3*B3)
      BL=B1/B4
      BM=B2/B4
      BN=B3/B4
      B1=XYZ(1,3)-XYZ(1,1)
      B2=XYZ(2,3)-XYZ(2,1)
      B3=XYZ(3,3)-XYZ(3,1)
      C4=SQRT(B1*B1+B2*B2+B3*B3)
      EL=B1/C4
      EM=B2/C4
      EN=B3/C4
      IF(B4.GE.C4) GO TO 10
      BAL(2,1)=EL
      BAL(2,2)=EM
      BAL(2,3)=EN
      GO TO 20
10    BAL(2,1)=BL
      BAL(2,2)=BM
      BAL(2,3)=BN
20    B1=EN*BM-EM*BN
      B2=EL*BN-BL*EN
      B3=BL*EM-BM*EL
      B4=SQRT(B1*B1+B2*B2+B3*B3)
      BAL(3,1)=B1/B4
      BAL(3,2)=B2/B4
      BAL(3,3)=B3/B4
      XX=BAL(2,2)*BAL(3,3)-BAL(2,3)*BAL(3,2)
      Y=BAL(3,1)*BAL(2,3)-BAL(2,1)*BAL(3,3)
      Z=BAL(2,1)*BAL(3,2)-BAL(3,1)*BAL(2,2)
      A=SQRT(XX*XX+Y*Y+Z*Z)
      BAL(1,1)=XX/A
      BAL(1,2)=Y/A
      BAL(1,3)=Z/A

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE MBYL    !!!!!!!!!!!!!!!
!     FUNCTION                   计算面板应力
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MBYL(KK,ID)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      DIMENSION BAL(3,3),BAK(3,3),BAD(3,3),BLA(3,3),SZ(6)
      COMMON /A15/STRZ(MAXN,6),EPGZ(MAXN,6)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /A8/ME(MAXN)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /ST/STR(MAXN,6),EPG(MAXN,6)
      COMMON /NJX/NJJ,MBJ,HDMAX,hdam(100),hhy(100),ii0(100),HHY2(100)
      COMMON /STRM/SSS(MAXN,3)

      IF(ID.EQ.0 .and. mype.eq.0 )
     +  WRITE(31,*) ' K     STR(n,/,z)   (unit:t/m/m)   NM='
      IF(ID.EQ.0 .and. mype.eq.0 )
     +  WRITE(31,*)'法向  顺坡向    坝轴向'
      CALL MBFXYL(BAL)

      DO 30 I=1,3
      DO 30 J=1,3
        BLA(I,J)=BAL(J,I)
30    CONTINUE
      DO 40 I=1,3
      DO 40 J=1,3
        BAL(I,J)=BLA(I,J)
40    CONTINUE

      NUM = mnode(1)
      KE=II0(KK)
      DO 60 K=1,NUM
        if( iE_lg(K) > KE ) cycle
c        NM=ME(K)
        NM=node(K*nnode(1))
        GAM=AE(1,NM)
        DF=AE(9,NM)
        G=AE(7,NM)
        IF(GAM.LT.2.35) GO TO 50
        IF(G.GT.1.0E-5.OR.DF.GT.1.0E-5) GO TO 50
        DO 10 J=1,6
          IF(ID.EQ.0) SZ(J)=STRZL(K,J)
          IF(ID.EQ.1) SZ(J)=STRL(K,J)
10      CONTINUE
        BAK(1,1)=SZ(1)
        BAK(2,2)=SZ(2)
        BAK(3,3)=SZ(3)
        BAK(1,2)=SZ(4)
        BAK(2,3)=SZ(5)
        BAK(1,3)=SZ(6)
        BAK(2,1)=SZ(4)
        BAK(3,2)=SZ(5)
        BAK(3,1)=SZ(6)
        CALL CH29MB(BAL,BAK,BAD)
        IF(ID.EQ.0) WRITE(31,100) iE_lg(K),SZ(1:6),NM
        IF(ID.EQ.0) WRITE(31,100) iE_lg(K),(BAD(I,I),I=1,3),(SSSL(K,J),J=1,3), NM
        sssl(k,1)=BAD(1,1)
        sssl(k,2)=BAD(2,2)
        sssl(k,3)=BAD(3,3)

50      CONTINUE

        IF(ID.EQ.0) cycle
        STRL(K,1)=BAD(1,1)
        STRL(K,2)=BAD(2,2)
        STRL(K,3)=BAD(3,3)
        STRL(K,4)=BAD(1,2)
        STRL(K,5)=BAD(2,3)
        STRL(K,6)=BAD(1,3)
60    CONTINUE

100   FORMAT(I6,6f10.3,I6)

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE MBFX    !!!!!!!!!!!!!!!!
!     FUNCTION   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MBFXYL(BAL)
      implicit real*8(a-h,o-z)
      DIMENSION BAL(3,3)
      COMMON /FA/AN,HDZ

      DO 10 I=1,3
      DO 10 J=1,3
        BAL(I,J)=0.0
10    CONTINUE
      BAL(3,3)=1.0
      T=SQRT(AN*AN+1.0)
      BAL(1,1)=1.0/T
      BAL(1,2)=-AN/T
      BAL(2,1)=AN/T
      BAL(2,2)=1.0/T

      RETURN
      END
      
      SUBROUTINE  CH29MB(BAL,BAK,BAD)
      implicit real*8(a-h,o-z)
      DIMENSION BAL(3,3),BLA(3,3),BAD(3,3),BAK(3,3),BAC(3,3)

      DO 10 I=1,3
      DO 10 J=1,3
        BLA(I,J)=BAL(J,I)
10    CONTINUE
      DO 30 I=1,3
      DO 30 J=1,3
        BAC(I,J)=0.0
        DO 35 II=1,3
          BAC(I,J)=BAC(I,J)+BLA(I,II)*BAK(II,J)
35      CONTINUE
30    CONTINUE
      DO 40 I=1,3
      DO 40 J=1,3
        BAD(I,J)=0.0
        DO 42 II=1,3
          BAD(I,J)=BAD(I,J)+BAC(I,II)*BAL(II,J)
42      CONTINUE
40    CONTINUE

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE STIF6    !!!!!!!!!!!!!!!
!     FUNCTION       6结点单元劲度矩阵
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE STIF6(KA,JK,E,U,GAM,XYZ,III)
      implicit real*8(a-h,o-z)
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,8),RJX(3,3),Q(3,8),BV(24),D(9),RST(2)
      COMMON /CMN5/FUN(8),P(3,8),XJR(3,3)
      COMMON /CMN4/NEE,NET,NRW,NDE
      COMMON /CMN3/SKE(24,24),RF(24),CE(24)
      D1=E*(1.0-U)/((1.0+U)*(1.0-2.0*U))
      D2=E*U/((1.0+U)*(1.0-2.0*U))
      D3=E*0.5/(1.0+U)
      DO 5 I=1,24
      RF(I)=0.0
      CE(I)=0.0
      DO 5 J=1,24
      SKE(I,J)=0.0
5     CONTINUE
      R=1.0/3.0
      S=R
      RST(1)=-0.5773503
      RST(2)=-RST(1)
      RSTH=0.5
      DO 100 IT=1,2
      T=RST(IT)
      CALL RMSD(XYZ,Q,DET,R,S,T,RJX,III)
      DO 35 I=1,6
      LJ=(I-1)*3
      ALS=-FUN(I)*RSTH*DET*GAM
      DO 35 J=1,3
      IJ=LJ+J
      CE(IJ)=CE(IJ)+ALS/9.80665
35      CONTINUE
20      IF(NET.EQ.0) GO TO 40
      IF(KA.LE.JK) GO TO 40
      DO 30 I=1,6
      J=I*3-1
      RF(J)=RF(J)-FUN(I)*RSTH*DET*GAM
30      CONTINUE
40      K3=0
      DO 10 I=1,6
      K3=K3+3
      K2=K3-1
      K1=K2-1
      BV(K1)=Q(1,I)
      BV(K2)=Q(2,I)
      BV(K3)=Q(3,I)
10      CONTINUE
      DO 70 I=1,18
      DO 60 J=I,18
      SKE(I,J)=SKE(I,J)+BV(I)*BV(J)*DET*RSTH
60    CONTINUE
70    CONTINUE
80    CONTINUE
90    CONTINUE
100   CONTINUE
      DO 110 I=2,18
      M=I-1
      DO 110 J=1,M
      SKE(I,J)=SKE(J,I)
110   CONTINUE
      KK=-3
      DO 130 I=1,6
      KK=KK+3
      K1=KK+1
      K2=KK+2
      K3=KK+3
      DO 125 J=1,6
      LL=(J-1)*3
      L1=LL+1
      L2=L1+1
      L3=L2+1
      IC=0
      DO 120 K=1,3
      M=KK+K
      DO 115 L=1,3
      NL=L+LL
      IC=IC+1
      D(IC)=SKE(M,NL)
115   CONTINUE
120   CONTINUE
      SKE(K1,L1)=D(1)*D1+(D(5)+D(9))*D3
      SKE(K2,L2)=D(5)*D1+(D(1)+D(9))*D3
      SKE(K3,L3)=D(9)*D1+(D(5)+D(1))*D3
      SKE(K1,L2)=D(2)*D2+D(4)*D3
      SKE(K2,L1)=D(4)*D2+D(2)*D3
      SKE(K2,L3)=D(6)*D2+D(8)*D3
      SKE(K3,L2)=D(8)*D2+D(6)*D3
      SKE(K3,L1)=D(7)*D2+D(3)*D3
      SKE(K1,L3)=D(3)*D2+D(7)*D3
125   CONTINUE
130   CONTINUE
      RETURN
      END
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE FAUVW    !!!!!!!!!!!!!!
!     FUNCTION      缝的错动位移 
!!!   警告!!!!!!!!!!!!!!!!!!!!
!     本子函数没有调用，所以没有进行并行化处理
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FAUVW(KK)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8),IO(MAXN),Q(MAXN,3)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /A2/IPE(MAXN,8)
      COMMON /A8/ME(MAXN)
      COMMON /A14/UVW(MAXN,3)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A16/DZ(MAXN,3)
      COMMON /NJX/NJJ,MBJ,HDMAX,hdam(100),hhy(100),ii0(100),HHY2(100)
      COMMON /EJOINT/IJKL(MAXN),VML(MAXN)
      WRITE(16,*) ' K   FANGuvw(short,loang,n) (unit:cm)  NM' 
      IQ=0
      KE=II0(KK)
      DO 100 K=1,KE
      NM=ME(K)
      T=AE(1,NM)
      IF(ABS(T-0.1).LT.1.0E-5) GO TO 100
      III=0
      DO 20 I=1,8
      IL=IPE(K,I)
      IF(IL.EQ.0) GO TO 20
      III=III+1
      DO 10 J=1,3
      XYZ(J,I)=COP(IL,J)
10      CONTINUE
20      CONTINUE
      V=VML(K)
      IF(ABS(V).GE.1.0E-5) GO TO 100
      IQ=IQ+1
      IO(IQ)=K
      DO 45 I=1,3
      Q(IQ,I)=DZ(K,I)*100.0
45      CONTINUE
100   CONTINUE
      WRITE(16,120) (IO(I),(Q(I,J),J=1,3),ME(IO(I)),I=1,IQ)
120   FORMAT(I6,3E12.4,I3)
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE OUTPUT    !!!!!!!!!!!!!!
!     FUNCTION            输出计算结果
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OUTPUT(IFO,KK)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION UVW(MAXN,3),INP(MAXN),C1(3),C2(3),L1(3),L2(3)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /A8/ME(MAXN)
      COMMON /A5/R(MAXN1)
      COMMON /A15/STRZ(MAXN,6),EPGZ(MAXN,6)
      COMMON /STRM/SSS(MAXN,3)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A4/N,NH,MX,JR(3,MAXN)
      COMMON /NJX/NJJ,MBJ,HDMAX,hdam(100),hhy(100),ii0(100),HHY2(100)
      COMMON /A14/UVW1(MAXN,3)

      IF(IFO.EQ.5) GO TO 50
      DO 705 J=1,3
        C1(J)=0.0
        C2(J)=0.0
705   CONTINUE
c      DO 710 I=1,NP
      DO 710 I=1,knode
      DO 710 J=1,3
        IF(UVW1(I,J).GT.C1(J)) L1(J)=I
        IF(UVW1(I,J).GT.C1(J)) C1(J)=UVW1(I,J)
        IF(UVW1(I,J).LT.C2(J)) L2(J)=I
        IF(UVW1(I,J).LT.C2(J)) C2(J)=UVW1(I,J)
710   CONTINUE
      WRITE(16,*) ' NP, UVWMAX,UVWMIN (unit:  m)'
      WRITE(16,563) (L1(I),C1(I),I=1,3)
      WRITE(16,563) (L2(I),C2(I),I=1,3)
563   FORMAT(3(I8,F13.3))

      IK=0
c      DO 240 I=1,NP
      DO 240 I=1,knode
        DO 16 J=1,3
          IF(ABS(UVW1(I,J)).GE.0.0) GO TO 230 
16      CONTINUE
        GO TO 240
230     CONTINUE
        IK=IK+1
        INP(IK)=I
        DO 10 J=1,3
          UVW(IK,J)=UVW1(I,J)*100.0
10      CONTINUE
240   CONTINUE
      WRITE(16,90)
90    FORMAT('  I=','   Ux     ','     Uy     ','    Uz (unit: cm)')
      WRITE(16,100) (INP(I),(UVW(I,J),J=1,3),I=1,IK)
100   FORMAT(I5,3F15.5)
      WRITE(16,80)
80    FORMAT(1X,'  K      STR(K,6)         SSS(K,3) (unit: t/m/m)')
      KE=II0(KK)
      DO 30 K=1,KE
        MN=ME(K)
        T=AE(1,MN)
        WRITE(16,110) K,(STRZ(K,I),I=1,6),(SSS(K,J),J=1,3)
30    CONTINUE
110   FORMAT(1X,I6,6F15.2,3F15.2)
50    CONTINUE

      RETURN
      END

      SUBROUTINE OUTPUT1(IFO,KK)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION UVW(MAXN,3),INP(MAXN),C1(3),C2(3),L1(3),L2(3)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /A8/ME(MAXN)
      COMMON /A5/R(MAXN1)
      COMMON /A15/STRZ(MAXN,6),EPGZ(MAXN,6)
      COMMON /STRM/SSS(MAXN,3)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A4/N,NH,MX,JR(3,MAXN)
      COMMON /NJX/NJJ,MBJ,HDMAX,hdam(100),hhy(100),ii0(100),HHY2(100)
      COMMON /A14/UVW1(MAXN,3)
      IF(IFO.EQ.5) GO TO 50
      DO 705 J=1,3
      C1(J)=0.0
      C2(J)=0.0
705   CONTINUE
      DO 710 I=1,NP
      DO 710 J=1,3
      IF(UVW1(I,J).GT.C1(J)) L1(J)=I
      IF(UVW1(I,J).GT.C1(J)) C1(J)=UVW1(I,J)
      IF(UVW1(I,J).LT.C2(J)) L2(J)=I
      IF(UVW1(I,J).LT.C2(J)) C2(J)=UVW1(I,J)
710   CONTINUE
      WRITE(55,*) ' THE COMPUTING STEP',KK     
      WRITE(55,*) ' NP, UVWMAX,UVWMIN (unit:  m)'
      WRITE(55,563) (L1(I),C1(I),I=1,3)
      WRITE(55,563) (L2(I),C2(I),I=1,3)
563   FORMAT(3(I8,F13.3))
      IK=0
      DO 240 I=1,NP
      DO 16 J=1,3
      IF(ABS(UVW1(I,J)).GE.0.0) GO TO 230 
16    CONTINUE
      GO TO 240
230   CONTINUE
      IK=IK+1
      INP(IK)=I
      DO 10 J=1,3
      UVW(IK,J)=UVW1(I,J)*100.0
10      CONTINUE
240   CONTINUE
      WRITE(55,90)
90    FORMAT('  I=','   Ux     ','     Uy     ','    Uz (unit: cm)')
      WRITE(55,100) (INP(I),(UVW(I,J),J=1,3),I=1,IK)
100   FORMAT(I5,3F15.5)
      WRITE(55,80)
80    FORMAT(1X,'  K      STR(K,6)         SSS(K,3) (unit: t/m/m)')
      KE=II0(KK)
      DO 30 K=1,KE
      MN=ME(K)
      T=AE(1,MN)
      WRITE(55,110) K,(STRZ(K,I),I=1,6),(SSS(K,J),J=1,3)
30      CONTINUE
110   FORMAT(1X,I6,6F15.2,3F15.2)
50    CONTINUE
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE STRESS    !!!!!!!!!!!!!!
!     FUNCTION         计算单元的应力和应变
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE STRESS(KK)
      use ComData
      use FEMData
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8),RJX(3,3),Q(3,8)
      DIMENSION UVW(3,8),X(9),C(6),MM(4),BP(6,3)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /ST/STR(MAXN,6),EPG(MAXN,6)
      COMMON /A2/IPE(MAXN,8)
      COMMON /A8/ME(MAXN)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A5/R(MAXN1)
      COMMON /CMN4/NEE,NET,NRW,NDE
      COMMON /A4/N,NH,MX,JR(3,MAXN)
      COMMON /NJX/NJJ,MBJ,HDMAX,hdam(100),hhy(100),ii0(100),HHY2(100)
      COMMON /A13/ET(MAXN),UT(MAXN),GAMT(MAXN)
      COMMON/B11/BI(4),CI(4),DI(4)
      COMMON /EJOINT/IJKL(MAXN),VML(MAXN)
      DATA  MM/1.0,-1.0,1.0,-1.0/

      num = mnode(1)
      !初始化单元应力STRL和EPGL
      if( .not. allocated(strl) ) allocate(strl(num,6), epgl(num,6))

c      DO 10 I=1,NE
      DO 10 I=1,NUM
      DO 10 J=1,6
        STRL(I,J)=0.0D0
        EPGL(I,J)=0.0D0
10    CONTINUE

      KE=II0(KK)
c      DO 150 K=1,KE
      DO 150 K=1,num
        if( iE_lg(K) .gt. KE) goto 150
        NEE=K           !单元编号
        DO 20 I=1,8
c          IV=IPE(K,I)
          IV=node((k-1)*nnode(1)+I)
          IF(IV.EQ.0) GO TO 20
          DO 22 J=1,3   
c            XYZ(J,I)=COP(IV,J)
            XYZ(J,I)=COOR((IV-1)*3+J)
22        CONTINUE
20      CONTINUE
        III=IJKL(K)     !单元有效节点个数
c        V1=VML(K)
        V1=VOL(K)       !单元体积
        IF(V1.LT.1.0E-8) GO TO 100
c        NM=ME(K)
        NM=node(K*nnode(1))       !单元材料编号
        E=ETL(K)
        U=UTL(K)
        D1=E*(1.0-U)/((1.0+U)*(1.0-2.0*U))
        D2=E*U/((1.0+U)*(1.0-2.0*U))
        D3=E*0.5/(1.0+U)
        DO 30 I=1,8
c          IV=IPE(K,I)
          IV=node((K-1)*nnode(1)+I)
          IF(IV.EQ.0) GO TO 30
          IF(IV.gt.knode .or. IV.lt.0 ) then
            write(*,*) 'Error! IV, knode=', IV, knode
            call my_endjob(ierr)
          ENDIF
          DO 32 J=1,3
c            IW=JR(J,IV)
            IW=JRG(J,IV)
            IF(IW.EQ.0) UVW(J,I)=0.0
c            IF(IW.NE.0) UVW(J,I)=R(IW)
            IF(IW.NE.0) UVW(J,I)=disp(j,iv)
32        CONTINUE
30      CONTINUE
        IF(III.EQ.4) GO TO 200
        T=0.0
        IF(III.EQ.6) GO TO 36
        R1=0.0
        S=0.0
        GO TO 45
36      R1=0.3333333
        S=R1
45      CALL RMSD(XYZ,Q,DET,R1,S,T,RJX,III)
        DO 35 I=1,9
35      X(I)=0.0
        DO 40 I=1,8
c          IV=IPE(K,I)
          IV=node((K-1)*nnode(1)+I)
          IF(IV.EQ.0) GO TO 40
          X(1)=X(1)+Q(1,I)*UVW(1,I)
          X(2)=X(2)+Q(2,I)*UVW(2,I)
          X(3)=X(3)+Q(3,I)*UVW(3,I)
          X(4)=X(4)+Q(1,I)*UVW(2,I)
          X(5)=X(5)+Q(1,I)*UVW(3,I)
          X(6)=X(6)+Q(2,I)*UVW(3,I)
          X(7)=X(7)+Q(2,I)*UVW(1,I)
          X(8)=X(8)+Q(3,I)*UVW(1,I)
          X(9)=X(9)+Q(3,I)*UVW(2,I)
c          if( mype.eq.1 .and. k.eq.1 ) then
c            write(*,*) 'X8:', i,x(8),Q(3,i),uvw(1,i)
c            write(*,*) 'X9:', i,x(9),Q(3,i),uvw(2,i)
c          endif
40      CONTINUE
c        if( mype.eq.1 .and. k.eq.1 ) then
c          write(*,*) "X:", (x(i),i=1,9)
c        endif
        C(1)=D1*X(1)+D2*(X(2)+X(3))
        C(2)=D1*X(2)+D2*(X(3)+X(1))
        C(3)=D1*X(3)+D2*(X(1)+X(2))
        C(4)=D3*(X(7)+X(4))
        C(5)=D3*(X(9)+X(6))
        C(6)=D3*(X(8)+X(5))
        DO 80 I=1,6
          STRL(K,I)=-C(I)
80      CONTINUE
        DO 90 I=1,3
90      EPGL(K,I)=-X(I)
        EPGL(K,4)=-(X(7)+X(4))
        EPGL(K,5)=-(X(9)+X(6))
        EPGL(K,6)=-(X(8)+X(5))
        GO TO 100
200     CALL FOUR(K,V)
        DO 210 J=1,6
210     X(J)=0.0
        DO 220 I=1,3
        DO 220 J=1,6
220     BP(J,I)=0.0
        DO 230 J=1,4
          BP(1,1)=MM(J)*BI(J)
          BP(2,2)=MM(J)*CI(J)
          BP(3,3)=MM(J)*DI(J)
          BP(4,2)=BP(1,1)
          BP(6,3)=BP(1,1)
          BP(4,1)=BP(2,2)
          BP(5,3)=BP(2,2)
          BP(5,2)=BP(3,3)
          BP(6,1)=BP(3,3)
          DO 230 L=1,6
            DO 230 M=1,3
              X(L)=X(L)+BP(L,M)*UVW(M,J)/V/6.0
230     CONTINUE
        C(1)=D1*X(1)+D2*(X(2)+X(3))
        C(2)=D1*X(2)+D2*(X(3)+X(1))
        C(3)=D1*X(3)+D2*(X(1)+X(2))
        C(4)=D3*X(4)
        C(5)=D3*X(5)
        C(6)=D3*X(6)
        DO 240 I=1,6
          STRL(K,I)=-C(I)
240       EPGL(K,I)=-X(I)
100     CONTINUE
150   CONTINUE

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE  MAIN   !!!!!!!!!!!!!!!!
!     FUNCTION      计算主应力或主应变   0－应变   1－应力
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MAIN(IM)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8),SW(6),PS(3)
      COMMON /ST/STR(MAXN,6),EPG(MAXN,6)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /STRM/SSS(MAXN,3)
      COMMON /A2/IPE(MAXN,8)
      COMMON /D13/GMV(MAXN)
      COMMON /EJOINT/IJKL(MAXN),VML(MAXN)
      
      NUM = mnode(1)
      if( .not. allocated(SSSL)) allocate(SSSL(num,3), GMVL(num))

c      DO 100 K=1,NE
      DO 100 K=1,NUM         !单元循环
        ION=0
        III=0
        DO 20 I=1,8
c          IV=IPE(K,I)
          IV=node((K-1)*nnode(1) + I)
          IF(IV.EQ.0) GO TO 20
          DO 10 J=1,3
c10        XYZ(J,I)=COP(IV,J)
10        XYZ(J,I)=COOR((IV-1)*3+J)
          III=III+1
20      CONTINUE
c        V=VML(K)
        V=VOL(K)
        IF(ABS(V).LE.1.0E-5) GO TO 100
        DO 30 I=1,3
          PS(I)=0.0
30      CONTINUE
        IF(IM.EQ.1) GO TO 60
        DO 40 I=1,3
c40      SW(I)=EPG(K,I)*1.0E6
40      SW(I)=EPGL(K,I)*1.0D6
        DO 50 I=4,6
c50      SW(I)=EPG(K,I)*0.5*1.0E6
50      SW(I)=EPGL(K,I)*0.5*1.0D6
        GO TO 80
60      DO 70 I=1,6
c          SW(I)=STR(K,I)
          SW(I)=STRL(K,I)
          IF(SW(I).GT.1.0E5) ION=1
70      CONTINUE
c        if( k.eq.67 .and. mype.eq.1 ) then
c          write(*,*) 'In MAIN, 67---SW=',(SW(i),i=1,6)
c        endif

80      B=-SW(1)-SW(2)-SW(3)
        C=SW(1)*SW(2)+SW(3)*SW(2)+SW(1)*SW(3)-SW(4)*SW(4)-SW(5)*SW(5)
     #     -SW(6)*SW(6)
        D=SW(1)*SW(5)*SW(5)+SW(2)*SW(6)*SW(6)+SW(3)*SW(4)*SW(4)
     #    -2.0*SW(4)*SW(5)*SW(6)-SW(1)*SW(2)*SW(3)
        P=C-B*B/3.0
        Q=2.0*B*B*B/27.0-B*C/3.0+D
        R=0.25*Q*Q+P*P*P/27.0
        IF(ABS(R).GT.0.5E-11) GO TO 110
        FACE=1.0
        IF(Q.LT.0.0) FACE=-1.0
        Y=(ABS(Q)*0.5)**(1.0/3.0)
        Y=Y*FACE
        PS(1)=2.0*Y-B/3.0
        PS(2)=-Y-B/3.0
        PS(3)=PS(2)
        GO TO 120
110     S=SQRT(-P/3.0)
        T=-0.5*Q/S/S/S
        IF(ABS(T).GT.1.0) T=SIGN(1.0,-Q)
        T=1.5707963267949-ASIN(T)
        Y=2.0*S*COS(T/3.0)
        PS(1)=Y-B/3.0
        Y=2.0*S*COS((T+6.2831853071796)/3.0)
        PS(2)=Y-B/3.0
        Y=2.0*S*COS((T+2.0*6.2831853071796)/3.0)
        PS(3)=Y-B/3.0
        IF(PS(3).LT.PS(2)) GO TO 120
        SIII=PS(2)
        PS(2)=PS(3)
        PS(3)=SIII
        IF(PS(3).LT.PS(1)) GO TO 120
        SII=PS(1)
        PS(1)=PS(3)
        PS(3)=SII
120     CONTINUE
        IF(IM.NE.1) GO TO 150
        DO 140 I=1,3
c          SSS(K,I)=PS(I)
          SSSL(K,I)=PS(I)
140     CONTINUE
c        if( K.eq.67 .and. mype.eq.1 ) then
c          write(*,*) 'IN main, 67---SSSL=',(SSSL(k,i),i=1,3)
c        endif
        GO TO 100
c150     GMV(K)=ABS(PS(1)-PS(3))*1.0E-6
150     GMVL(K)=ABS(PS(1)-PS(3))*1.0E-6
100   CONTINUE
170   FORMAT(I6,E15.6)
190   FORMAT(I6,6E12.4)

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE MJZDDH    !!!!!!!!!!!!!!
!     FUNCTION          将不属于本级荷载的编号剔除
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MJZDDH(KK,NO)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      COMMON /NJX/NJJ,MBJ,HDMAX,hdam(100),hhy(100),ii0(100),HHY2(100)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A4/N,NH,MX,JR(3,MAXN)
      COMMON /A2/IPE(MAXN,8)
      ! IPE为全局单元节点信息，子区域用node数组取代
      integer, dimension(:), allocatable :: node0


      NEJK=II0(KK)                   !本加载步需要加载的单元数
      num = mnode(1)

      allocate(node0(knode))
      do i=1,knode
        node0(i) = 0
      enddo

      ! 找出本加载步需要激活的节点
      NO = 0
      knode0 = 0
      num0 = 0
      do i=1,num
        if( ie_lg(i) .gt. NEJK ) cycle
        num0 = num0 + 1
        do j=1, 8
          iv = node((i-1)*nnode(1)+j)
          if( iv .gt. knode_i )  cycle
          if( node0(iv) .eq. 0) then
            node0(iv) = 1
            knode0 = knode0 + 1
            if( iv .gt. NO ) NO=IV
          endif
        enddo
      enddo

      ! 剔除本加载步不加载的单元包含的节点
      do i=1, knode
        if( node0(i) .eq. 1) cycle
        JRL(1,I)=0
        JRL(2,I)=0
        JRL(3,I)=0
      enddo

      deallocate(node0)

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE MR    !!!!!!!!!!!!!!!
!     FUNCTION              形成结点自由度列表
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MR(KK)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      COMMON /A10/JC(100001)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A4/N,NH,MX,JR(3,MAXN)
      COMMON /NJX/NJJ,MBJ,HDMAX,hdam(100),hhy(100),ii0(100),HHY2(100)

      integer :: JRR(3,5000)
      integer :: Ns, Nt
      integer :: AllocateStatus

      if( .not. allocated(JRL) ) then
        allocate(JRL(3,knode), STAT = AllocateStatus)
        IF (AllocateStatus .NE. 0) STOP "* Not enough memory for JRL*"
      endif
      if( .not. allocated(JRG) ) then
        allocate(JRG(3,knode), STAT = AllocateStatus)
        IF (AllocateStatus .NE. 0) STOP "* Not enough memory for JRG*"
      endif

c      DO 10 I=1,NP
      DO 10 I=1,knode
      DO 10 J=1,3
        JRL(J,I) = 1
        JRG(J,I) = 0
10    CONTINUE

!     形成内部节点自由度列表
c      DO 20 I=1,NR
      DO 20 I=1,knode_i
c        K=JC(I)
        K=JC(iN_lg(I)) ! 节点自由度约束信息,iN_lg-整体节点编号
        J=K/1000       ! 节点编号
        L=(K-J*1000)/100
        M=(K-J*1000-L*100)/10
        N=K-J*1000-L*100-M*10
        do idx=1,knode_i
          if( iN_lg(idx) .eq. J ) exit
        enddo
        JRL(1,idx)=L
        JRL(2,idx)=M
        JRL(3,idx)=N
20    CONTINUE
      CALL MJZDDH(KK,NO)
      N=0
      DO 70 I=1,NO
      DO 60 J=1,3
c        IF(JR(J,I)) 60,60,50
c        if(jr(j,i) <= 0) goto 60
c        if(jr(j,i) > 0) goto 50
        if(JRL(j,i) > 0) THEN
50        N=N+1
          JRL(J,I)=N
        ENDIF
60    CONTINUE
70    CONTINUE

      if( mype .eq. 0) then
        nupdate(0) = 0
        nupdate(1) = N
        do jpe=1,npes-1
          call My_recvai(0,jpe,nupdate(jpe+1),1)
          nupdate(jpe+1) = nupdate(jpe+1)+nupdate(jpe)
        enddo
        do jpe=1,npes-1
          call My_sendai(jpe,0,nupdate(jpe),1)
        enddo
        Ns = nupdate(mype)
      else
        call My_sendai(0,mype,N,1)
        call My_recvai(mype,0,Ns,1)
      endif
      
      N_update = N
      if( allocated(update) ) deallocate(update)
      allocate(update(N), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "* allocate update error *"

      if( allocated(update_index) ) deallocate(update_index)
      allocate(update_index(N), STAT = AllocateStatus)
      IF (AllocateStatus .NE. 0) STOP "* allocate update_index error *"
     
      ! 计算内部节点的局部变量号(JRL)和整体变量号(JRG)
      Nt = 0
      do i=1, NO
        do j=1, 3
          if( JRL(j,i) > 0 ) then
            Nt = Nt + 1
            if( Nt <=0 .or. Nt > N ) then
               print *,'Error!!i,j,nt=', i,j,nt
               call My_endjob(ierr)
            endif
            JRG(j,i) = JRL(j,i) + Ns
            update(Nt) = JRG(j,i)
            update_index(Nt) = JRL(j,i)
          endif
        enddo
      enddo
      if(Nt .ne. N_update) then
        write(*,*) "Error,Nt, N_update=", Nt, N_update
        call My_endjob(ierr)
      endif

      ! 从文件读入每个节点的整体变量编号(仅测试，JRG文件由串行程序生成给出)
      JRG4 = 0
      if( JRG4 .eq. 1) then
        write(*,*)'mype,NP=', mype,NP
        open(21,file='JRG',form='formatted',status='old')
        do i=1,NP
          read(21,*) itmp,(JRR(j,i),j=1,3)
        enddo
        close(21)
        write(*,*) 'READ JRG file ok.....', knode
        do i=1, knode_i
          if( iN_lg(i) .eq. 29 ) print *,'29---mype,i=',mype,i
          do j=1, 3
            JRG(j,i) = JRR(j,iN_lg(i))
          enddo
        enddo
      endif

      ! 从相邻子区域得到和整理边界扩展节点方程号
      N_external = 0
      write(ext,'(i5)') mype
      filename='conn_' // trim(adjustl(ext))
      open(21,file=filename,form='formatted',status='old')
      read(21,*) neigh
      do i=1, neigh
        read(21,*) isend, nsend
        allocate(neighs(nsend))
        read(21,*) (neighs(j),j=1,nsend)
        ntmp = nsend*kdgof
c        write(*,*) 'mype, isend, ntmp_s=', mype, isend, ntmp
        allocate(ipool(ntmp))
        do j=1,nsend
          inod = neighs(j)
          if(inod .gt. knode_i .and. inod .le. 0) then
            print *,'Error!! inod, knode=',inod, knode
            call My_endjob(ierr)
          endif
          do k=1,kdgof
            ipool((j-1)*kdgof+k) = JRG(k,inod)
          enddo
        enddo
        call My_sendai(isend,mype,ipool,ntmp)
        deallocate(neighs, ipool)
c
        read(21,*) irecv, nrecv
        allocate(neighr(nrecv))
        read(21,*) (neighr(j),j=1,nrecv)
        ntmp = nrecv*kdgof
        allocate(ipool(ntmp))
c        write(*,*) 'mype, irecv, ntmp_r=', mype, irecv, ntmp
        call My_recvai(mype,irecv,ipool,ntmp)
        do j=1, nrecv
          inod = neighr(j)
          if(inod .gt. knode .and. inod .le. 0) then
            print *,'Error!! inod, knode=',inod, knode
            call My_endjob(ierr)
          endif
          do k=1,kdgof
            if( (j-1)*kdgof+k .gt. ntmp ) then
              print *,'(j-1)*kdgof+k, ntmp=',(j-1)*kdgof+k, ntmp
            endif
            JRG(k,inod) = ipool((j-1)*kdgof+k)
            if( JRG(k,inod) > 0) then
              N_external = N_external + 1
            endif
          enddo
        enddo
        deallocate(neighr, ipool)
      enddo
      close(21)

      Nt = 0
      do i=1, NO
        do j=1, 3
          if( JRL(j,i) > 0 ) then
            Nt = Nt + 1
            if( Nt <=0 .or. Nt > N ) then
               print *,'Error!!i,j,nt=', i,j,nt
               call My_endjob(ierr)
            endif
c            JRG(j,i) = JRL(j,i) + Ns
            update(Nt) = JRG(j,i)
c            update_index(Nt) = JRL(j,i)
          endif
        enddo
      enddo
      if(Nt .ne. N_update) then
        write(*,*) "Error,Nt, N_update=", Nt, N_update
        call My_endjob(ierr)
      endif
c      write(*,*) mype,'--MR--N_update,N_external=',N_update,N_external

100   FORMAT('JRL=',4(4I4,2X))

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE NHMA0    !!!!!!!!!!!!!
!     FUNCTION            形成系数矩阵DCRS格式存储空间(na, ia, a)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE NHMA0(kk)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      COMMON /A6/MA(MAXN1)
      COMMON /A2/IPE(MAXN,8)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A4/N,NH,MX,JR(3,MAXN)
      COMMON /NJX/NJJ,MBJ,HDMAX,HDAM(100),HHY(100),II0(100),HHY2(100)
      COMMON /STRU3/ICOL(125930000)
      DIMENSION ISEQ(160,MAXN1),NSORT(159)

      integer, dimension(:), allocatable :: nai, naj
      integer ise

      ! 根据单元信息计算非零元个数

      NEJK = II0(kk)
      num = mnode(1)
      if( allocated(na) ) deallocate(na)
      allocate(na(N+1))

C     !每个自由度最多81个相关自由度（含自身）.
      maxbnd = 81
      maxt = maxbnd*N
      allocate(nai(maxt))

      do i=1, N+1
        na(i) = 0
      enddo
      do i=1, maxt
        nai(i) = 0
      enddo

      ! 初步形成非零元结构
      nextern = 0
      nnz = 0
      do IE=1,num
        if( ie_lg(IE) .gt. NEJK ) cycle
        do i=1, 8
          iv = node((IE-1)*nnode(1)+i)
          if( iv > knode) then
            write(*,*) 'Error!mype,iv,ie,i,knode,nnode(1)=',mype,iv,ie,i,knode,nnode(1)
            call My_endjob(ierr)
          endif
          if( iv == 0)  cycle
          do j=1, kdgof
            ijn = JRL(j,iv)
c            if( mype==0 .and. iv == 743 .and. j==1 ) 
c     +         write(*,'(a,5i4)') '0--743,ie,i,j,JRL,JRG=',ie,i,j,JRL(j,iv),JRG(j,iv)
            if( ijn == 0) cycle
            if( iv > knode_i) then
              write(*,*) 'Error!mype,iv,ie,knode,knode_i=',mype,iv,ie,knode,knode_i
              call My_endjob(ierr)
            endif
C           !找出每个自由度的关联自由度，并存放在NAI内。这里找出的是整体编号
            do i1=1, 8
              iw = node((IE-1)*nnode(1)+i1)
              do j1=1, kdgof
                ijn1 = JRG(j1,iw)
                if( ijn1 == 0 ) cycle
                iexist = 0
                do i2=1,na(ijn)
                  if( nai((ijn-1)*maxbnd+i2) == ijn1 ) then
                    iexist = 1
                    exit
                  endif
                enddo
                if( iexist == 1 ) cycle
                na(ijn) = na(ijn) + 1
                if( na(ijn) > maxbnd ) then
                   write(*,*) "Error! na .gt. Maxbnd", na(ijn), maxbnd
                   call My_endjob()
                endif
                nai((ijn-1)*maxbnd + na(ijn)) = ijn1
                nnz = nnz + 1
              enddo
            enddo
          enddo
        enddo
      enddo

      ! 整理非零元结构
      do i=N_update+1,2,-1
        na(i) = na(i-1)
      enddo
      if( allocated(ia) ) deallocate(ia, am)
      allocate(ia(nnz), am(nnz))

      nnz = 0
      na(1) = 0
      do i=1, N_update
        na(i+1) = na(i+1)+na(i)
        do j=1, maxbnd
          if( nai((i-1)*maxbnd+j) == 0 ) cycle
          nnz = nnz+1
          ia(nnz) = nai((i-1)*maxbnd+j)
        enddo
      enddo
c      write(*,*) 'mype,nnz,na =', mype, nnz, na(N_update+1)
      if( allocated(nai) ) deallocate(nai)

      do i=1, N_update
        n0 = na(i)+1
        n1 = na(i+1)
        allocate(nai(n1-n0+1), naj(n1-n0+1))
        do j=n0, n1
          nai(j-n0+1) = ia(j)
          naj(j-n0+1) = j-n0+1
        ENDDO
        iv = 0
        iu = n1-n0
        call qsort(nai, naj, iv, iu)
        do j=n0,n1
          ia(j) = nai(j-n0+1)
        enddo
        deallocate(nai, naj)
      ENDDO

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE NHMA    !!!!!!!!!!!!!
!     FUNCTION             主对角线元素的指示矩阵
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE NHMA(kk)
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      COMMON /A6/MA(MAXN1)
      COMMON /A2/IPE(MAXN,8)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A4/N,NH,MX,JR(3,MAXN)
      COMMON /NJX/NJJ,MBJ,HDMAX,HDAM(100),HHY(100),II0(100),HHY2(100)
      COMMON /STRU3/ICOL(125930000)
      DIMENSION ISEQ(160,MAXN1),NSORT(159)
            
      DO 10 I=1,N
        MA(I)=0
10    CONTINUE       
C     !每个自由度最多81个相关自由度（含自身）.
C     !ISEQ（82，）存放该自由度相关自由度总数，以方便下面进行排序。
      DO 15 I=1,N
      DO 15 J=1,160
        ISEQ(J,I)=0
15    CONTINUE       
       
      DO 20 IE=1,II0(kk)
      DO 35 I=1,8
        IV=IPE(IE,I)
        IF(IV.EQ.0) GOTO 35
        DO 30 J=1,3
          IJN=JR(J,IV)
          IF(IJN.EQ.0) GOTO 30
C         !找出每个自由度的相关自由度，并存放在ISEQ内。
          DO 55 I1=1,8
            IW=IPE(IE,I1)
            IF(IW.EQ.0) GOTO 55   
            DO 50 J1=1,3
              IJN1=JR(J1,IW)
              IF(IJN1.EQ.0) GOTO 50
              IF(IJN1.LT.IJN) THEN
                L=ISEQ(160,IJN)
                IF(L.EQ.0) THEN    ! 第一个相关自由度
                  ISEQ(1,IJN)=IJN1
                  ISEQ(160,IJN)=1
                  GOTO 50
                ENDIF
                DO 25 I2=1,L      ! 第2，3。。。个自由度
                  L1=ISEQ(I2,IJN)
                  IF(IJN1.EQ.L1) GOTO 55  ! IW节点已在之前的一个单元循环过,此时
25              CONTINUE
                L=L+1
                ISEQ(L,IJN)=IJN1
                ISEQ(160,IJN)=L
!               WRITE(100,*)IJN,IE,ISEQ(82,IJN)
              ENDIF
50          CONTINUE
55        CONTINUE
30      CONTINUE
35    CONTINUE
20    CONTINUE
C     !所有单元循环已经结束，下面对每个自由度的相关自由度由小到大进行排序。
      
      MA(1)=1            ! 存放该行及前面行元素总数
      ICOL(1)=1
      DO 70 I=2,N        ! 对每个自由度循环
        L=ISEQ(160,I)    ! 对于每个自由度来说，比他小的相关自由度的个数（不包括自身）
        IF(L.EQ.0) THEN   !!只有自己
          MA(I)=MA(I-1)+1                  ! +1指对角元
          ICOL(MA(I))=I
          GOTO 70
        ENDIF
        ! 由小到大排序并重新赋值
        DO 71 J=1,L
          NSORT(J)=ISEQ(J,I)
71      CONTINUE
        CALL SORT(NSORT,L)
        DO 72 J=1,L
          ISEQ(J,I)=NSORT(J)
72      CONTINUE
C        !获取指示矩阵
        MA(I)=MA(I-1)+L+1        !存放该行及前面行元素总数
        DO 73 J=1,L
          J1=MA(I-1)+J
          ICOL(J1)=ISEQ(J,I)
73      CONTINUE
        ICOL(MA(I))=I
70    CONTINUE
      
      NH=MA(N)      
      WRITE(*,*)'NH=',MA(N),'N=',N,'NE=',II0(KK)
      IF(MA(N).NE.NH) STOP
      IF(KK.EQ.NJJ)THEN
        WRITE(65,*)NH   
        WRITE(65,100)(ICOL(I),I=1,NH)
        WRITE(75,100)(MA(I),I=1,N)
      ENDIF
100   FORMAT(10(2X,I10))
      RETURN
      END

C-------------------  ------------------------------
      SUBROUTINE SORT(NA,N)
      implicit real*8(a-h,o-z)
      DIMENSION NA(159)
      DO 5 I=1,N-1
      DO 5 J=I+1,N
        IF(NA(I).GT.NA(J)) THEN
          T=NA(I)
          NA(I)=NA(J)
          NA(J)=T
        ENDIF
5     ENDDO
C         
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE FPJD    !!!!!!!!!!!!!!!
!     FUNCTION   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE   FPJD(R,S,T,DET,XYZ)
      implicit real*8(a-h,o-z)
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION    XYZ(3,8),XI(8),ETA(8),ZETA(8)
      COMMON /CMN5/FUN(8),P(3,8),XJR(3,3)
      COMMON /CMN4/NEE,NET,NRW,NDE
      DATA XI/-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0/
      DATA ETA/-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0/
      DATA ZETA/-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0/

      DO 10 I=1,8
        XI0=XI(I)*R
        ETA0=ETA(I)*S
        ZETA0=ZETA(I)*T
        FUN(I)=0.125*(1.0+XI0)*(1.0+ETA0)*(1.0+ZETA0)
        P(1,I)=0.125*XI(I)*(1.0+ETA0)*(1.0+ZETA0)
        P(2,I)=0.125*ETA(I)*(1.0+XI0)*(1.0+ZETA0)
        P(3,I)=0.125*ZETA(I)*(1.0+XI0)*(1.0+ETA0)
 10   CONTINUE
      DO 40 I=1,3
        DO 30 J=1,3
          DET=0.0
          DO 20 K=1,8
            DET=DET+P(I,K)*XYZ(J,K)
 20       CONTINUE
          XJR(I,J)=DET
 30     CONTINUE
 40   CONTINUE
      DET=XJR(1,1)*XJR(2,2)*XJR(3,3)
     *     +XJR(2,1)*XJR(3,2)*XJR(1,3)
     *     +XJR(1,2)*XJR(2,3)*XJR(3,1)
     *     -XJR(1,3)*XJR(2,2)*XJR(3,1)
     *     -XJR(3,2)*XJR(2,3)*XJR(1,1)
     *     -XJR(1,2)*XJR(2,1)*XJR(3,3)
      IF(DET .LT. 1.0E-8) GO TO 100
Cc    WRITE(16,500) DET
      RETURN
 500  FORMAT(10X,'DET=***',F15.3)
 100  WRITE(16,600) NEE,R,S,T
      WRITE(16,500) DET
 600  FORMAT(1X,'FPJD ERROR***  ELEMENT(',I5,')',4X,'R=',F10.5,4X,'S=',
     *    F10.3,' T=',F10.5)
c      STOP
      call My_endjob(ierr)
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE RMSD    !!!!!!!!!!!!!!
!     FUNCTION          Jacobi矩阵的逆阵
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE RMSD(XYZ,Q,DET,R,S,T,RJX,III)
      implicit real*8(a-h,o-z)
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,8),RJX(3,3),Q(3,8)
      COMMON /CMN5/FUN(8),P(3,8),XJR(3,3)
      COMMON /CMN4/NEE,NET,NRW,NDE

      IF(III.EQ.8) CALL FPJD(R,S,T,DET,XYZ)
      IF(III.EQ.6) CALL FPJD6(R,S,T,DET,XYZ)
      REC=1.0/DET
      RJX(1,1)=REC*( XJR(2,2)*XJR(3,3)-XJR(2,3)*XJR(3,2))
      RJX(2,1)=REC*(-XJR(2,1)*XJR(3,3)+XJR(2,3)*XJR(3,1))
      RJX(3,1)=REC*( XJR(2,1)*XJR(3,2)-XJR(2,2)*XJR(3,1))
      RJX(1,2)=REC*(-XJR(1,2)*XJR(3,3)+XJR(1,3)*XJR(3,2))
      RJX(2,2)=REC*( XJR(1,1)*XJR(3,3)-XJR(1,3)*XJR(3,1))
      RJX(3,2)=REC*(-XJR(1,1)*XJR(3,2)+XJR(1,2)*XJR(3,1))
      RJX(1,3)=REC*( XJR(1,2)*XJR(2,3)-XJR(1,3)*XJR(2,2))
      RJX(2,3)=REC*(-XJR(1,1)*XJR(2,3)+XJR(1,3)*XJR(2,1))
      RJX(3,3)=REC*( XJR(1,1)*XJR(2,2)-XJR(1,2)*XJR(2,1))
      DO 70 I=1,III
        DO 60 J=1,3
          X=0.0
          DO 50 K=1,3
            X=X+RJX(J,K)*P(K,I)
 50       CONTINUE
          Q(J,I)=X
 60     CONTINUE
 70   CONTINUE
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBRUTINE FPJD6    !!!!!!!!!!!!!
!     FUNCTION       积分点上的
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE   FPJD6(R,S,T,DET,XYZ)
      implicit real*8(a-h,o-z)
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION    XYZ(3,8),ZETA(6)
      COMMON /CMN5/FUN(8),P(3,8),XJR(3,3)
      COMMON /CMN4/NEE,NET,NRW,NDE
      DATA ZETA/-1.0,-1.0,-1.0,1.0,1.0,1.0/

      DO 5 I=1,8
        FUN(I)=0.0
        DO 5 J=1,3
          P(J,I)=0.0
5     CONTINUE  
      DO 10 J=1,2
        K=(J-1)*3
        DO 10 L=1,3
          I=K+L
          ZETA0=ZETA(I)*T
          IF(L.EQ.1) FUN(I)=0.5*R*(1.0+ZETA0)
          IF(L.EQ.2) FUN(I)=0.5*S*(1.0+ZETA0)
          IF(L.EQ.3) FUN(I)=0.5*(1.0-R-S)*(1.0+ZETA0)
10    CONTINUE
      DO 15 J=1,2
        DO 15 K=1,3
          I=(J-1)*3+K
          ZETA0=ZETA(I)*T
          IF(K.EQ.1) P(1,I)=0.5*(1.0+ZETA0)
          IF(K.EQ.1) P(2,I)=0.0
          IF(K.EQ.1) P(3,I)=0.5*ZETA(I)*R
          IF(K.EQ.2) P(1,I)=0.0
          IF(K.EQ.2) P(2,I)=0.5*(1.0+ZETA0)
          IF(K.EQ.2) P(3,I)=0.5*ZETA(I)*S
          IF(K.EQ.3) P(1,I)=-0.5*(1.0+ZETA0)
          IF(K.EQ.3) P(2,I)=-0.5*(1.0+ZETA0)
          IF(K.EQ.3) P(3,I)=0.5*ZETA(I)*(1.0-R-S)
15    CONTINUE
      DO 40 I=1,3
        DO 30 J=1,3
          DET=0.0
          DO 20 K=1,6
            DET=DET+P(I,K)*XYZ(J,K)
 20       CONTINUE
          XJR(I,J)=DET
 30     CONTINUE
 40   CONTINUE
      DET=XJR(1,1)*XJR(2,2)*XJR(3,3)
     *     +XJR(2,1)*XJR(3,2)*XJR(1,3)
     *     +XJR(1,2)*XJR(2,3)*XJR(3,1)
     *     -XJR(1,3)*XJR(2,2)*XJR(3,1)
     *     -XJR(3,2)*XJR(2,3)*XJR(1,1)
     *     -XJR(1,2)*XJR(2,1)*XJR(3,3)
      IF(DET .LT. 1.0E-5) GO TO 100
      RETURN
 500  FORMAT(10X,'DET=***',F15.3)
 100  WRITE(16,600) NEE,R,S,T
      WRITE(16,500) DET
 600  FORMAT(1X,'FPJD6 ERROR***  ELEMENT(',I5,')',4X,'R=',F10.5,4X,'S=',
     *    F10.3,' T=',F10.5)
c      STOP
      call My_endjob(ierr)
      END

      ! 计算单元体积
      SUBROUTINE SOLVEVOLUME
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8)
      COMMON /EJOINT/IJKL(MAXN),VML(MAXN)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /A2/IPE(MAXN,8)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A8/ME(MAXN)

      num = mnode(1)
      allocate(vol(num))
c      DO 100 K=1,NE
      DO 100 K=1,NUM
c        IM=ME(K)
        IM = node(K*nnode(1))
        VML(K)=0.0
50      III=0
        DO 20 I=1,8
c          IV=IPE(K,I)
          IV = node((k-1)*nnode(1)+I)
          IF(IV.EQ.0) GO TO 20
          III=III+1
c          XYZ(1:3,I)=COP(IV,1:3)
          do j=1,3
            xyz(j,I) = coor((iv-1)*3+j)
          enddo
20      CONTINUE
40      IJKL(K)=III
        IF(AE(1,IM).LE.0.3)GOTO 101
         ! VOLUME OF CONTACT ELEMENT AND JOINT ELEMENT EQUAL TO 0
        CALL V68(XYZ,III,V)
        V1=ABS(V)
c        VML(K)=V1
        VoL(K)=V1
101     CONTINUE
        IF(AE(1,IM).gE.0.3 .AND. V1.LE.0.0)THEN
          WRITE(*,*)'WRONG ELEMENT',K
          STOP
        ENDIF
c        write(55,*)k,im,vml(k)
100   CONTINUE

      do i=1, num
        VML(iE_lg(i)) = Vol(i)
      enddo
      
      if( mype .eq. 0) then
        do iprc=1, npes-1
          call My_recvai(mype,iprc,ndata,1)
          allocate(rdata(ndata),idata(ndata))
          call My_recvai(mype,iprc,idata,ndata)
          call My_recvar(mype,iprc,rdata,ndata)
          do i=1, ndata
            VML(idata(i)) = rdata(i)
          enddo
          deallocate(idata,rdata)
        enddo
        do i=1, ne
          write(55,*) i, me(i), VML(i)
        enddo
      else
        call My_sendai(0,mype,num,1)
        call My_sendai(0,mype,iE_lg,num)
        call My_sendar(0,mype,vol,num)
      endif

      RETURN
      END

      SUBROUTINE EXAMCONTACTELEMENT
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
      DIMENSION XYZ(3,8),VELEMEN(MAXN),NUMNODE(MAXN)
      COMMON /A7/COP(MAXN,3),AE(11,40)
      COMMON /EJOINT/IJKL(MAXN),VML(MAXN)
      COMMON /A2/IPE(MAXN,8)
      COMMON /A3/NE,NP,NR,NERW,NERWDOWN
      COMMON /A8/ME(MAXN)
      REAL LL1,LL2,LL3

      num = mnode(1)
c      DO 100 K=1,NE
      DO 100 K=1,Num
c        IM=ME(K)
        IM = node(K*nnode(1))
        IF(AE(1,IM) .LE. 0.2)THEN
          DO 20 I=1,8
c            IV=IPE(K,I)                                           
            IV = node((k-1)*nnode(1)+I)
            IF(IV.EQ.0) GO TO 20
c            XYZ(1:3,I)=COP(IV,1:3)
            do j=1,3
              xyz(j,I) = coor((iv-1)*3+j)
            enddo
20        CONTINUE
          NUMNODE(K)=IJKL(K)
c          VELEMEN(K)=VML(K)
          VELEMEN(K)=Vol(K)
      
          IF(NUMNODE(K).EQ.4)THEN
            STOP 'ERROR'
          ENDIF
      
          IF(NUMNODE(K).EQ.8)THEN
c            L1=IPE(K,1)
c            L2=IPE(K,2)
c            L3=IPE(K,3)
c            L4=IPE(K,4)
c            L5=IPE(K,5)
c            L6=IPE(K,6)
c            L7=IPE(K,7)
c            L8=IPE(K,8)
c            LL1=SQRT((COP(L1,1)-COP(L2,1))**2+(COP(L1,2)-COP(L2,2))**2+
c     *               (COP(L1,3)-COP(L2,3))**2)
c            LL2=SQRT((COP(L1,1)-COP(L3,1))**2+(COP(L1,2)-COP(L3,2))**2+
c     *               (COP(L1,3)-COP(L3,3))**2)
c            LL3=SQRT((COP(L1,1)-COP(L5,1))**2+(COP(L1,2)-COP(L5,2))**2+
c     *               (COP(L1,3)-COP(L5,3))**2)
c   
            L1=node((k-1)*nnode(1)+1)
            L2=node((k-1)*nnode(1)+2)
            L3=node((k-1)*nnode(1)+3)
            L4=node((k-1)*nnode(1)+4)
            L5=node((k-1)*nnode(1)+5)
            L6=node((k-1)*nnode(1)+6)
            L7=node((k-1)*nnode(1)+7)
            L8=node((k-1)*nnode(1)+8)
            LL1=SQRT((Coor((L1-1)*kcoor+1)-COOR((L2-1)*kcoor+1))**2
     *              +(COOR((L1-1)*kcoor+2)-COOR((L2-1)*kcoor+2))**2
     *              +(COOR((L1-1)*kcoor+3)-COOR((L2-1)*kcoor+3)))**2
            LL2=SQRT((Coor((L1-1)*kcoor+1)-COOR((L3-1)*kcoor+1))**2
     *              +(COOR((L1-1)*kcoor+2)-COOR((L3-1)*kcoor+2))**2
     *              +(COOR((L1-1)*kcoor+3)-COOR((L3-1)*kcoor+3)))**2
            LL3=SQRT((Coor((L1-1)*kcoor+1)-COOR((L5-1)*kcoor+1))**2
     *              +(COOR((L1-1)*kcoor+2)-COOR((L5-1)*kcoor+2))**2
     *              +(COOR((L1-1)*kcoor+3)-COOR((L5-1)*kcoor+3)))**2
   
            IF(((LL3-LL1).GT.0.00001).OR.((LL3-LL2).GT.0.00001))THEN
              WRITE(*,*)'CONTACT ELEMENT IS FAULT',k 
              STOP 'DISTANCE ERROR'
            ENDIF

          ENDIF

          IF(NUMNODE(K).EQ.6)THEN
c            L1=IPE(K,1)
c            L2=IPE(K,2)
c            L3=IPE(K,3)
c            L4=IPE(K,4)
c            L5=IPE(K,5)
c            L6=IPE(K,6)
c            L7=IPE(K,7)
c            L8=IPE(K,8)
c            LL1=SQRT((COP(L1,1)-COP(L2,1))**2+(COP(L1,2)-COP(L2,2))**2+
c     *               (COP(L1,3)-COP(L2,3))**2)
c            LL2=SQRT((COP(L1,1)-COP(L3,1))**2+(COP(L1,2)-COP(L3,2))**2+
c     *               (COP(L1,3)-COP(L3,3))**2)
c            LL3=SQRT((COP(L1,1)-COP(L4,1))**2+(COP(L1,2)-COP(L4,2))**2+
c     *               (COP(L1,3)-COP(L4,3))**2)

            L1=node((k-1)*nnode(1)+1)
            L2=node((k-1)*nnode(1)+2)
            L3=node((k-1)*nnode(1)+3)
            L4=node((k-1)*nnode(1)+4)
            L5=node((k-1)*nnode(1)+5)
            L6=node((k-1)*nnode(1)+6)
            L7=node((k-1)*nnode(1)+7)
            L8=node((k-1)*nnode(1)+8)
            LL1=SQRT((Coor((L1-1)*kcoor+1)-COOR((L2-1)*kcoor+1))**2
     *              +(COOR((L1-1)*kcoor+2)-COOR((L2-1)*kcoor+2))**2
     *              +(COOR((L1-1)*kcoor+3)-COOR((L2-1)*kcoor+3)))**2
            LL2=SQRT((Coor((L1-1)*kcoor+1)-COOR((L3-1)*kcoor+1))**2
     *              +(COOR((L1-1)*kcoor+2)-COOR((L3-1)*kcoor+2))**2
     *              +(COOR((L1-1)*kcoor+3)-COOR((L3-1)*kcoor+3)))**2
            LL3=SQRT((Coor((L1-1)*kcoor+1)-COOR((L4-1)*kcoor+1))**2
     *              +(COOR((L1-1)*kcoor+2)-COOR((L4-1)*kcoor+2))**2
     *              +(COOR((L1-1)*kcoor+3)-COOR((L4-1)*kcoor+3)))**2
   
            IF(((LL3-LL1).GT.0.00001).OR.((LL3-LL2).GT.0.00001))THEN
              WRITE(*,*)'CONTACT ELEMENT IS FAULT' ! 对面距离比其它面距离大
              STOP 'DISTANCE ERROR'
            ENDIF
          ENDIF
        ENDIF !AE(1,IM).GT.0.2
100   CONTINUE

      RETURN
      END
      
      SUBROUTINE SSORPCG
      USE ComData
      USE FEMData
      implicit real*8(a-h,o-z)
C** IN THIS SUBROUTINE, THE STORAGE MODE OF THE MATRIX A IS AS FOLLOWS:
C** 1. THE LOW TRIANGLE MATRIX IS STORED;
C** 2. THE NONZERO ELEMENTS ARE STORED;
C** 3. THE ELEMENTS ARE STORED ACCORDING THE ROW;
      COMMON /A1/SK(125930000),DM(MAXN1)
      COMMON /A6/MA(MAXN1)
      COMMON /A4/N,NH,MX,JR(3,MAXN)
      COMMON /STRU3/ICOL(125930000)
      COMMON /A5/R(MAXN1)
      
      DIMENSION C(MAXN1),D(MAXN1),E(MAXN1),X(MAXN1)
      DIMENSION G(MAXN1),Y(MAXN1),Z(MAXN1) 
            
      OMIGA=1.2
      OMIGA1=(2-OMIGA)/OMIGA
      EPS=1E-2
      RLAMDA=1E3
      
      DO 1 I=1,N
        R(I)=R(I)*RLAMDA
1     CONTINUE
      
      KIter=0
      DO I=1,N
       Y(I)=0.0
       X(I)=0.0
C** G=AX-B
       G(I)=-R(I)
       D(I)=0.0
      ENDDO
      DO 10 I=1,N
C** Y(I)=W(-1)G
        IF(I.EQ.1) THEN
C** for the first diagonal, the coresponding 
C**   nonzero non-diagonal element
c** whose dof is less than 1 is zero
          L1=1
        ELSE
          L1=MA(I-1)+1
        ENDIF
        L2=MA(I)-1
        DO 15 J=L1,L2
          IF(ICOL(J).GE.I) THEN
            WRITE(*,*)'ERROR1'
          ENDIF
          G(I)=G(I)-SK(J)*Y(ICOL(J))
15      CONTINUE
        Y(I)=G(I)*OMIGA/SK(MA(I))
C** Z(I)=-VY
        Z(I)=-OMIGA1*(SK(MA(I)))*Y(I)
10    CONTINUE  
C** COMPUTE THE VECTOR D=W(-T)Z
      DO 20 I=N,1,-1
        IF(I.EQ.1) THEN
C** for the first diagonal, the coresponding
C**  nonzero non-diagonal element 
c** whose dof is less than 1 is zero
          L1=1
        ELSE
          L1=MA(I-1)+1
        ENDIF
        L2=MA(I)-1
C** FOR I equals to N, D(I) equals to 0
        D(I)=OMIGA*(Z(I)-D(I))/(SK(MA(I)))
C** D(I) HAS BEEN SET TO ZERO IN LOOP 10
        DO 25 J=L1,L2
          D(ICOL(J))=D(ICOL(J))+SK(J)*D(I)
25      CONTINUE
20    CONTINUE
C** CHECK THE CONVERGENCE
60    DELTA=0.0
        
      DO I=1,N
        DELTA=DELTA+Y(I)**2*OMIGA1*(SK(MA(I)))
      ENDDO
      IF(DELTA.LE.EPS) THEN
        DO I=1,N
          R(I)=X(I)/RLAMDA
        ENDDO
        WRITE(*,*)'END ITERATION STEP  ', KIter,'TORLERANCE  ',DELTA1
        RETURN
      ENDIF
C*    * SET THE TEMPERORY VARIBLE ZER
      DELTA1=0.0
      DO I=1,N
        C(I)=Z(I)-OMIGA1*SK(MA(I))*D(I)
        DELTA1=DELTA1+D(I)*(Z(I)+C(I))
      ENDDO
      
      T=DELTA/DELTA1
      DELTA1=0.0
      DO I=1,N
        E(I)=0.0
      ENDDO
      DO 40 I=1,N
        X(I)=X(I)+T*D(I)
C** E(I)=W(-1)C
        
        IF(I.EQ.1) THEN
C** for the first diagonal, the coresponding
C**  nonzero non-diagonal element 
c** whose dof is less than 1 is zero
          L1=1
        ELSE
          L1=MA(I-1)+1
        ENDIF
        L2=MA(I)-1
        DO 45 J=L1,L2
          IF(ICOL(J).GE.I) THEN
            WRITE(*,*)'ERROR1'
          ENDIF
          C(I)=C(I)-SK(J)*E(ICOL(J))
45      CONTINUE
        E(I)=C(I)*OMIGA/(SK(MA(I)))
        Y(I)=Y(I)+T*(D(I)+E(I))
        DELTA1=DELTA1+Y(I)**2*OMIGA1*(SK(MA(I)))
40    CONTINUE
            
      BETA0=DELTA1/DELTA
      DO I=1,N
        Z(I)=-OMIGA1*(SK(MA(I)))*Y(I)+BETA0*Z(I)
        D(I)=0.0
      ENDDO
C** D(K+1)=W(-T)Z(K+1)
      DO 50 I=N,1,-1
        IF(I.EQ.1) THEN
C** for the first diagonal, the coresponding
C**  nonzero non-diagonal element 
c** whose dof is less than 1 is zero
          L1=1
        ELSE
          L1=MA(I-1)+1
        ENDIF
        L2=MA(I)-1
C** FOR I equals to N, D(I) equals to 0
        D(I)=OMIGA*(Z(I)-D(I))/(SK(MA(I)))

C** D(I) HAS BEEN SET TO ZERO IN LOOP 10
        DO 55 J=L1,L2
          D(ICOL(J))=D(ICOL(J))+SK(J)*D(I)
55      CONTINUE
50    CONTINUE
      KIter=KIter+1
      IF(MOD(KIter,50).EQ.0) THEN
        WRITE(*,*)'ITERATION STEP IS ', KIter,'TORLERANCE IS ',DELTA1
      ENDIF
      GOTO 60
C           
      RETURN
      END

