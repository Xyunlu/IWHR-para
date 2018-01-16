!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    THREE DIMENSION FINITE ELEMENT PROGRAM OF earth EMBANKMENT            !
!                        DYNAMIC ANALYSIS                                    !
!                             BY LHJ(26/5/2015) in IWHR                                    !
!                         覆盖层材料容重赋予0.0!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PROGRAM earthdamDYNAMIC
      COMMON /A3/NE,NP,NR
      COMMON /GAUSS/RSTG(3),H(3)
      COMMON /A15/STRZ(500000,6)
      COMMON /STRM/SSS(500000,3)
      COMMON /STM/STRMM(500000,3)
      COMMON /A13/ET(500000),UT(500000),GAMT(500000) 
      COMMON /ECON/ETC(500000),UTC(500000)
      COMMON /STRU3/ICOL(125930000)
      COMMON /A6/MA(990000)
      COMMON /A8/ME(500000)
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /A4/N,NH,MX,JR(3,500000)
      COMMON /ELE/IPEYW(500000,8)
      COMMON /A2/IPE(500000,8)
      COMMON /FSELE/NUMFS,IPLX(20000)
      COMMON /FACTOR11/ACCEFACTOR,STRAINFACTOR,EPS2
      COMMON /MBJOINT/MB,IP(50000)
      COMMON /SLABADDMASS/IPSLAB(50000),NUH,MP(50000)

      OPEN(7,file='DYNAMICP.DAT',status='old')
      OPEN(9,file='tc1.DAT',status='old')
      OPEN(11,FILE='EV1.DAT',STATUS='OLD')          ! static analysis
      OPEN(15,FILE='INPUT.DAT',STATUS='OLD')        !static analysis
      OPEN(21,FILE='EWATER.DAT',STATUS='OLD')       !slab additional mass
      OPEN(22,FILE='FSELEMENT.DAT',STATUS='OLD')    !dynamic stability
      OPEN(30,FILE='KMATRIXINDEX.DAT',STATUS='OLD') !static analysis
      OPEN(40,FILE='MA.DAT',STATUS='OLD')           !static analysis
      OPEN(48,FILE='WATERENODE.DAT',STATUS='OLD')   !static analysis
      !!!!!!!!!index of nonzero k matrix

      OPEN(8,file='A1.dat')

      OPEN(16,FILE='DA.DAT')
      OPEN(18,file='DIS.dat')
      OPEN(19,file='AXYZ.dat')

      !TYPICAL NODE AND ELEMENT TIME HISTORY!
      OPEN(50,FILE='AXYZ6037.DAT')
      OPEN(51,FILE='AXYZ5986.DAT')
      OPEN(52,FILE='AXYZ5921.DAT')
      OPEN(53,FILE='TMAX6929.DAT')
      OPEN(54,FILE='TMAX4059.DAT')
      OPEN(55,FILE='TMAX952.DAT')
      !TYPICAL NODE AND ELEMENT TIME HISTORY!

      !TIMEHISTORY STABILITY ANALYSIS!
      OPEN(31,FILE='STRESSDL.DAT')
      !TIMEHISTORY STABILITY ANALYSIS!


      !ANSYS POST RESULTS!
      OPEN(61,FILE='AGANSYS.DAT')
      OPEN(62,FILE='STANSYS.DAT')
      OPEN(63,FILE='UDANSYS.DAT')
      !ANSYS POST RESULTS!

      !PERMANENT DEFORMATION COMPUTING!
      open(100,file='A2SHEN.dat')
      open(200,file='A2IWHR.dat')
      !PERMANENT DEFORMATION COMPUTING!

      H(1)=0.5555555556
      H(2)=0.8888888889
      H(3)=H(1)
      RSTG(1)=-0.774596669241
      RSTG(2)=0.0
      RSTG(3)=-RSTG(1)
      !加速度放大倍数!
      ACCEFACTOR=1
      !加速度放大倍数!
      !剪应变折减系数!
      STRAINFACTOR=1.0
      !剪应变折减系数!
      EPS2=0.1
      !收敛控制标准
      WRITE(*,*)'CFRD OR COREDAM DYNAMIC COMPUTING'
      WRITE(*,*)'面板材料静力参数中第九项为0'
!      WRITE(*,*)'用于确定面板附加质量'
      READ(48,*)NUH
      write(*,*)'UPSTREAM WATER ELEMENT',nuh
      IF(NUH.EQ.0)GOTO 80
      DO I=1,NUH
      READ(48,*)MP(I),IPSLAB(I)
      ENDDO
80    CONTINUE
      READ(22,*)NUMFS
      DO I=1,NUMFS
      READ(22,*)IPLX(I)
      ENDDO
      CALL INPUT(NNM)
      !IN ORDER TO DEFINE THE FJZL
      DO I=1,NE
      READ(21,*)IK,(IPEYW(IK,J),J=1,8),IJ
      ENDDO
      !IN ORDER TO DEFINE THE FJZL
      CALL MR
!!!!!!!!replace the old NHMA subroutine!!!!!!!!!!!!!!!
      READ(30,*)NH
      READ(30,*)(ICOL(I),I=1,NH)
      READ(40,*)(MA(I),I=1,N)
!!!!!!!!replace the old NHMA subroutine!!!!!!!!!!!!!!!
      CALL SOLVEVOLUME
      WRITE(*,*)'PRE IS FINISHED'

      READ(11,*) (IK,ET(K),UT(K),GAMT(K),K=1,NE)
      DO 10 IK=1,NE
      ETC(K)=ET(K)
      UTC(K)=UT(K)
      READ(11,*) K,(STRZ(K,IJ),IJ=1,6),(SSS(K,J),J=1,3)
10      CONTINUE
      WRITE(*,*) 'STATIC  STRESS READ END'
       CLOSE(11)
       CLOSE(30)
       CLOSE(40)
      CALL AMOUNT
      IF(MB.NE.0)THEN
      CALL MWI
      WRITE(*,*) 'slab FJZL END'
      ENDIF
!     存在防渗墙是需注意，静力参数中第9项不能为0，将面板单元静力参数第9项设为0

      DO 210 K=1,NE
      IF(ET(K).LT.250000.0) ET(K)=250000.0 !initial dynamic modulus
      DO 210 J=1,3
      STRMM(K,J)=SSS(K,J)
210      CONTINUE
      
      CALL GDFIRS
      NET=0
      WRITE(*,*)'DYNAMIC START'
      CALL  WILJG(NNM)
       CLOSE(16)
      STOP
      END

C      AMOUNT OF THE CONCRETE FACED AND JIONTS
      SUBROUTINE AMOUNT
      DIMENSION XYZ(3,8)
      COMMON /EJOINT/IJKL(500000),VML(500000)
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /A2/IPE(500000,8)
      COMMON /MBJOINT/MB,IP(50000)
      COMMON /A3/NE,NP,NR
      COMMON /D15/IG(50000),DU1(500000,3)
      COMMON /A8/ME(500000)/IFANG/IHJ
      MB=0
      DO 10 K=1,NE
      NM=ME(K)
      GAM=AE(1,NM)
      DF=AE(9,NM)
      G=AE(7,NM)
      IF(AE(4,NM).GE.0.1) GOTO 10
      IF(GAM.LT.2.35) GOTO 10
      IF(G.GT.1.0E-5.OR.DF.GT.1.0E-5) GOTO 10
      MB=MB+1
      IP(MB)=K
10      CONTINUE

      IHJ=0
      DO 100 K=1,NE
      NM1=ME(K)
      IF(ABS(AE(1,NM1)-0.2).GT.1.0E-5) GO TO 100
      IHJ=IHJ+1
      IG(IHJ)=K
100      CONTINUE
      WRITE(16,*) ' IHJ,MB=',IHJ,MB
      WRITE(*,*) ' IHJ,MB=',IHJ,MB
110      format(10i6)
      RETURN
      END
c--------------------------------------------------------------
      SUBROUTINE  WILJG(NNM)
      DIMENSION SS(3),UU(3)
      DIMENSION ag(3,4001),na(3),ugtt(3)
      DIMENSION ILO(100),URD(990000,3)
      DIMENSION GMAX(500000),SMAX(500000,3),SMIN(500000,3)
      DIMENSION UVW1(500000,3),UVW2(500000,3),UMAX(500000,3),
     *UMIN(500000,3)
      DIMENSION STMB1(500000,3),STMB2(500000,3),
     *SSMB1(500000,3),SSMB2(500000,3)
      DIMENSION DU(500000,3),DU2(500000,3),DUMAX(500000,3),
     *DUMIN(500000,3)
      DIMENSION SSSM(500000,3),SSSN(500000,3)
      DIMENSION AGMA(500000,3),AGMI(500000,3)
      DIMENSION UVW(500000,3)
      DIMENSION C1(3),C2(3),C3(3),C4(3),L1(3),L2(3),L3(3),L4(3)
      DIMENSION AGX(500000,3),AGN(500000,3)
      DIMENSION INP(500000)
      DIMENSION ST1(500000),ST2(500000),SST1(500000),SST2(500000)
      COMMON /D11/KG,GK(9,40)
      COMMON /D12/DE(500000),OMGA
      COMMON /D14/GM(500000)
      COMMON /u012/u0(990000),u1(990000),u2(990000)
      COMMON /A1/SK(125930000),DM(990000)
      COMMON /A9/SC(125930000)
      COMMON /A6/MA(990000)
      COMMON /A3/NE,NP,NR
      COMMON /A4/N,NH,MX,JR(3,500000)
      COMMON /A5/R(990000)
      COMMON /D15/IG(50000),DU1(500000,3)
      COMMON /AB/seita,dt,ab0,ab1,ab2,ab3,ab4,ab5,ab6,ab7,ab8
      COMMON /STRM/SSS(500000,3)
      COMMON /ST/STR(500000,6),EPG(500000,6)
      COMMON /D13/GMV(500000)
      COMMON /RMD/RLD(16),W0(16)
      COMMON /NJX/NJJ,MBJ,HDMAX,HDAM(100),HHY(100),II0(100),HHY2(100)
      COMMON /EJOINT/IJKL(500000),VML(500000)
      COMMON /A8/ME(500000)
      COMMON /D16/DL(500000)
      COMMON /FA/ALFA,HDZ
      COMMON /MBJOINT/MB,IP(50000)/IFANG/IHJ
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /METRICAL/MER,NPOINT
      COMMON /RDE/RS(990000)
      COMMON /A15/STRZ(500000,6)
      COMMON /STM/STRMM(500000,3)
      COMMON /DAMPLY/NEROCK(30),DAMPC(20,50)
      COMMON /A13/ET(500000),UT(500000),GAMT(500000)
      COMMON /DJPLUS/SD(500000,5)
      COMMON /FSELE/NUMFS,IPLX(20000)
      COMMON /FACTOR11/ACCEFACTOR,STRAINFACTOR,EPS2
      FA=ATAN(1.0/ALFA)
      READ(7,*) KG,((GK(J,I),J=1,9),I=1,NNM)   
      write(8,710) KG
      write(8,720) ((GK(J,I),J=1,9),I=1,NNM)
      
710      FORMAT('*********** GK*********=',I5)
720      FORMAT(F9.2,F8.4,F6.2,F7.3,F8.3,2F10.1,2F8.3)
      READ(7,*) EPS1 !误差控制标准 
      WRITE(8,730) EPS1
730      FORMAT('     EPS1=',3F10.5)
      READ(7,*) MER,NPOINT 
      READ(7,*) (NEROCK(J),J=1,MER)
731   FORMAT(A10,30I3)
      WRITE(8,731) ' NEROCK=',(NEROCK(J),J=1,MER)
      READ(7,*) ((DAMPC(I,J),J=1,MER*2+1),I=1,NPOINT)
      WRITE(8,735) ((DAMPC(I,J),I=1,NPOINT),J=1,MER*2+1)
735      FORMAT(10X,11F10.6)
      read(9,*) dt,seita
      write(8,750) dt,seita
750      format(5x,'dt=',f10.6,'   seita=',f8.4)
      READ(9,*) NNDD,(ILO(I),I=1,NNDD)
      write(8,740) NNDD,(ILO(I),I=1,NNDD)
740      FORMAT(' NNDD=',I4,' ILO=',15I5)
      read(9,*) nt,(na(i),i=1,3)
      write(8,760) nt,(na(i),i=1,3)
760      format(3x,'nt=',i5,2x,'  na(1,2,3)=',3i4)

      do 10 i=1,nt
      do 10 j=1,3
10      ag(j,i)=0.0
      lt=3
      if(na(1).eq.1.AND.na(2).eq.0) lt=1
      if(na(2).eq.1.AND.na(3).eq.0) lt=2
      do i=1,nt
      read(9,*) (ag(j,i),j=1,lt) !unit:m/s/s
       do k=1,lt
        ag(k,i)=ag(k,i)*ACCEFACTOR !放大倍数
       enddo
      enddo
      close(7)
      close(9)

      tt=seita*dt
      ab0=6.0/tt/tt
      ab1=3.0/tt
      ab2=2.0*ab1
      ab3=tt/2.0
      ab4=ab0/seita
      ab5=-ab2/seita
      ab6=1.0-3.0/seita
      ab7=dt/2.0
      ab8=dt*dt/6.0
      write(16,790) seita,dt,ab0,ab1,ab2,ab3,ab4,ab5,ab6,ab7,ab8
790      format(' seita,dt,ab0,ab1,ab2,ab3,ab4,ab5,ab6,ab7,ab8='/
     #          6e13.5/5e13.5)

      DO 30 K=1,NE
30      DE(K)=0.0

      DO 40 I=1,N
      URD(I,1)=0.0
      URD(I,2)=0.0
      URD(I,3)=0.0
40      CONTINUE

      DO 50 K=1,NE
      GMAX(K)=0.0
      ST1(K)=0.0
      ST2(K)=0.0
      DO 50 J=1,3
      SMAX(K,J)=0.0
      SMIN(K,J)=0.0
50      CONTINUE

      DO 60 I=1,NP
      DO 60 J=1,3
      UMAX(I,J)=0.0
      UMIN(I,J)=0.0
      AGX(I,J)=0.0
      AGN(I,J)=0.0
60      CONTINUE

      DO 70 I=1,MB
      DO 70 J=1,3
      K=IP(I)
      SSMB1(K,J)=0.0
      SSMB2(K,J)=0.0
70      CONTINUE

      DO 80 I=1,IHJ
      DO 80 J=1,3
      K=IG(I)
      DUMAX(K,J)=0.0
      DUMIN(K,J)=0.0
80      CONTINUE

      G3=1.0
      ID=1
      DO 90 K=1,NE
      GM(K)=0.000002
90      CONTINUE
      CALL GDMODE
      WRITE(*,*) ' GDMODE END'
      OMGA=3.55
      DO 95 I=1,N
95      RS(I)=1.0
      DO 1000 JID=1,NNDD
      IMM=0
100      CALL KRGAM(NJJ)
      WRITE(*,*) 'KRGAM END'
      CALL DEIGV(1,IMM)
      OMGA=W0(1)
      WRITE(*,800) ID,G3,OMGA,IMM
      IF(IMM.EQ.1)WRITE(8,800) ID,G3,OMGA,IMM
800      format('     ID=',i4,'  G3,OMGA,IMM=',F10.4,F12.6,I5)
      CALL KRGAM(NJJ)
      do 110 i=1,n   
      j=ma(i)
      sk(j)=sk(j)+ab0*DM(i)
110      continue
      do 120 i=1,nh
      sk(i)=sk(i)+ab1*sc(i)
120      continue
      write(*,810)
810      FORMAT(' KRGAM      END')

      DO 130 K=1,NE
      GM(K)=0.0
      SST1(K)=0.0
      SST2(K)=0.0
      DO 130 J=1,3
      SSSM(K,J)=0.0
      SSSN(K,J)=0.0
130      CONTINUE
      do 140 i=1,n
      u0(i)=URD(I,1)
      u1(i)=URD(I,2)
      u2(i)=URD(I,3)
140      continue
      DO 150 I=1,NP
      DO 150 J=1,3
      UVW1(I,J)=0.0
      UVW2(I,J)=0.0
      AGMA(I,J)=0.0
      AGMI(I,J)=0.0
150      CONTINUE
      DO 160 I=1,IHJ
      DO 160 J=1,3
      K=IG(I)
      DU1(K,J)=0.0
      DU2(K,J)=0.0
160      CONTINUE
      DO 170 I=1,MB
      DO 170 J=1,3
      K=IP(I)
      STMB1(K,J)=0.0
      STMB2(K,J)=0.0
170      CONTINUE
      IF(ID.EQ.1) INT=1
      IF(ID.NE.1) INT=ILO(ID-1)+1
      NT=ILO(ID)      
      do 500 it=INT,nt
      do 180 i=1,3
      ugtt(i)=ag(i,it)+seita*(ag(i,it+1)-ag(i,it))
180      continue
      call RXY(ugtt,r)
      call SSORPCG
      do 190 i=1,n
      pp=u2(i)
      po=u1(i)
      u2(i)=ab4*(r(i)-u0(i))+ab5*u1(i)+ab6*u2(i)
      u1(i)=u1(i)+ab7*(u2(i)+pp)
      u0(i)=u0(i)+po*dt+ab8*(u2(i)+2.0*pp)
190      continue
      do 200 i=1,n
      r(i)=u0(i)
200      continue
      DO 210 I=1,NP
      DO 210 J=1,3
      IV=JR(J,I)
      IF(IV.EQ.0) UVW(I,J)=0.0
      IF(IV.NE.0) UVW(I,J)=R(IV)
210      CONTINUE
      DO 220 I=1,NP
      DO 220 J=1,3
      IF(UVW(I,J).GT.UVW1(I,J)) UVW1(I,J)=UVW(I,J)
      IF(UVW(I,J).LT.UVW2(I,J)) UVW2(I,J)=UVW(I,J)
220      CONTINUE
820      format(I5,3F8.3)
      DO 230 I=1,NP
      DO 230 J=1,3
      IV=JR(J,I)
      !!!!!!!!!!!!!ABSOLUTE ACCELERATION!!!!!!!!!!!!!
      IF(IV.EQ.0) UVW(I,J)=AG(J,IT+1)
      IF(IV.NE.0) UVW(I,J)=U2(IV)+AG(J,IT+1)
      !!!!!!!!!!!!!ABSOLUTE ACCELERATION!!!!!!!!!!!!!
230      CONTINUE
      DO 240 I=1,NP
      DO 240 J=1,3
      IF(UVW(I,J).GT.AGMA(I,J)) AGMA(I,J)=UVW(I,J)
      IF(UVW(I,J).LT.AGMI(I,J)) AGMI(I,J)=UVW(I,J)
240      CONTINUE
      write(*,830) IT,(UVW(6037,J),J=1,3)
      IF(IMM.EQ.1)THEN
      write(50,830) IT,(UVW(6037,J),J=1,3)
      write(51,830) IT,(UVW(5986,J),J=1,3)
      write(52,830) IT,(UVW(5921,J),J=1,3)
      ENDIF
830      format(I6,F10.2,F10.2,F10.2)
      CALL STRESS(NJJ)
      CALL MAIN(1) 
831   FORMAT(1X,5F10.3)
      IF(IMM.EQ.1)THEN
      WRITE(31,*)'THE TIEM STEP',IT
      DO I=1,NUMFS
      IL=IPLX(I)
      SD(IL,1)=STRZ(IL,1)+STR(iL,1)
      SD(IL,2)=STRZ(IL,2)+STR(iL,2)
      SD(IL,3)=STRZ(IL,4)+STR(iL,4)
      SD(IL,4)=STRMM(IL,1)+SSS(iL,1)
      SD(IL,5)=STRMM(IL,3)+SSS(iL,3)
      WRITE(31,831)(SD(IL,J),J=1,5)
      ENDDO
      ENDIF 
      IF(IHJ.GT.0)THEN
      DO 260 I=1,IHJ
      K=IG(I)
      III=IJKL(K)
      CALL  WHDUS(K,III,SS,UU)
      DO 280 J=1,3
      SSS(K,J)=SS(J)
      DU(K,J)=UU(J)*100.0 !displacement of joint unit is cm
      IF(DU(K,J).GT.DU1(K,J)) DU1(K,J)=DU(K,J)
      IF(DU(K,J).LT.DU2(K,J)) DU2(K,J)=DU(K,J)
280      CONTINUE
260      CONTINUE
      ENDIF
      DO 300 K=1,NE
      ST=0.5*(SSS(K,1)-SSS(K,3))
      IF(IMM.EQ.1)then
      ACCE=STR(K,4)
      IF(ABS(ACCE).LE.1E-4)ACCE=1.0
      IF(K.EQ.6929)write(53,*) IT,ST,ST*ACCE/(ABS(ACCE)),STR(K,4)
      IF(K.EQ.4059)write(54,*) IT,ST,ST*ACCE/(ABS(ACCE)),STR(K,4)
      IF(K.EQ.952)write(55,*) IT,ST,ST*ACCE/(ABS(ACCE)),STR(K,4)
      endif        
      IF(ST.GT.SST1(K)) SST1(K)=ST
      IF(ST.LT.SST2(K)) SST2(K)=ST
      DO 300 J=1,3
      IF(SSS(K,J).GT.SSSM(K,J)) SSSM(K,J)=SSS(K,J)
      IF(SSS(K,J).LT.SSSN(K,J)) SSSN(K,J)=SSS(K,J)
300      CONTINUE

      CALL MAIN(0)
      DO 310 K=1,NE
      IF(ABS(GMV(K)).GT.GM(K)) GM(K)=ABS(GMV(K))
310      CONTINUE
      IF(MB.EQ.0) GOTO 500
      CALL MBYL(1)      
      DO 320 I=1,MB
      K=IP(I)
      DO 340 J=1,3
      IF(STR(K,J).GT.STMB1(K,J)) STMB1(K,J)=STR(K,J)
      IF(STR(K,J).LT.STMB2(K,J)) STMB2(K,J)=STR(K,J)
340      CONTINUE
320      CONTINUE

500       CONTINUE
      CALL GDMODE
      KPH=0
      DO 360 K=1,NE
      IF(DL(K).GT.EPS1)  KPH=KPH+1
360      CONTINUE
      G3=FLOAT(KPH)/FLOAT(NE)
      IF(IMM.EQ.1)  GOTO 381
      IF(G3.GT.EPS2) GOTO 100
      IF(G3.LE.EPS2)IMM=IMM+1
      IF(G3.LE.EPS2.AND.IMM.EQ.1)GOTO 100
381   CONTINUE

      DO 370 K=1,NE
      IF(ABS(GM(K)).GT.GMAX(K)) GMAX(K)=ABS(GM(K))
      IF(ABS(SST1(K)).GT.ST1(K)) ST1(K)=ABS(SST1(K))
      IF(SST2(K).LT.ST2(K)) ST2(K)=SST2(K)
      DO 370 J=1,3
      IF(SSSM(K,J).GT.SMAX(K,J)) SMAX(K,J)=SSSM(K,J)
      IF(SSSN(K,J).LT.SMIN(K,J)) SMIN(K,J)=SSSN(K,J)
370      CONTINUE

      DO 380 I=1,IHJ
      DO 380 J=1,3
      K=IG(I)
      IF(DU1(K,J).GT.DUMAX(K,J)) DUMAX(K,J)=DU1(K,J)
      IF(DU2(K,J).LT.DUMIN(K,J)) DUMIN(K,J)=DU2(K,J)
      SMAX(K,J)=DUMAX(K,J)
      SMIN(K,J)=DUMIN(K,J) 
380      CONTINUE

      DO 390 I=1,NP
      DO 390 J=1,3
      IF(UVW1(I,J).GT.UMAX(I,J)) UMAX(I,J)=UVW1(I,J)
      IF(UVW2(I,J).LT.UMIN(I,J)) UMIN(I,J)=UVW2(I,J)
390      CONTINUE

      DO 400 I=1,MB
      DO 400 J=1,3
      K=IP(I)
      IF(STMB1(K,J).GT.SSMB1(K,J)) SSMB1(K,J)=STMB1(K,J)
      IF(STMB2(K,J).LT.SSMB2(K,J)) SSMB2(K,J)=STMB2(K,J)
      SMAX(K,J)=SSMB1(K,J)
      SMIN(K,J)=SSMB2(K,J)
400      CONTINUE

      DO 410 I=1,NP
      DO 410 J=1,3
      IF(AGMA(I,J).GT.AGX(I,J)) AGX(I,J)=AGMA(I,J)
      IF(AGMI(I,J).LT.AGN(I,J)) AGN(I,J)=AGMI(I,J)
410      CONTINUE

      DO 415 I=1,N
      URD(I,1)=U0(I)
      URD(I,2)=U1(I)
      URD(I,3)=U2(I)
415      CONTINUE

      ID=ID+1
      IF(ID.GT.NNDD) GO TO 420
1000      CONTINUE
420      CONTINUE
      DO 430 J=1,3
      C1(J)=0.0
      C2(J)=0.0
      C3(J)=0.0
      C4(J)=0.0
430      CONTINUE
      DO 440 I=1,NP
      DO 440 J=1,3
      IF(AGX(I,J).GT.C1(J)) L1(J)=I
      IF(AGX(I,J).GT.C1(J)) C1(J)=AGX(I,J)
      IF(AGN(I,J).LT.C2(J)) L2(J)=I
      IF(AGN(I,J).LT.C2(J)) C2(J)=AGN(I,J)
      IF(UMAX(I,J).GT.C3(J)) L3(J)=I
      IF(UMAX(I,J).GT.C3(J)) C3(J)=UMAX(I,J)
      IF(UMIN(I,J).LT.C4(J)) L4(J)=I
      IF(UMIN(I,J).LT.C4(J)) C4(J)=UMIN(I,J)
440      CONTINUE
      WRITE(8,*) ' NP, AGMAX,AGMIN,UMAX,UMIN (unit: m/s/s cm)'
      WRITE(8,840) (L1(I),C1(I),I=1,3)
      WRITE(8,840) (L2(I),C2(I),I=1,3)
      WRITE(8,840) (L3(I),C3(I)*100.0,I=1,3)
      WRITE(8,840) (L4(I),C4(I)*100.0,I=1,3)
840      FORMAT(' NP MAX.or.MIN=',3(I8,F13.3))
      WRITE(8,*) ' NP          {AG MAX}       {AG MIN} (unit:m/s/s)'
      IK=0
      DO 470 I=1,NP
      DO 450 JK=1,3
      JQ=JR(JK,I)
      IF(JQ.NE.0) GO TO 460
450      CONTINUE
      GO TO 470
460      IK=IK+1
      INP(IK)=I
470      CONTINUE
      WRITE(8,850) (I,(AGX(I,J),J=1,3),
     *                     (AGN(I,J),J=1,3),I=1,NP)
850      FORMAT(I5,6F7.2)
      WRITE(8,*) ' NP             {UMAX}         {UMIN} (unit:cm)'
      WRITE(8,860) (I,(UMAX(I,J)*100.0,J=1,3),
     *            (UMIN(I,J)*100.0,J=1,3),I=1,NP)
860      FORMAT(I5,6F15.2)
      WRITE(8,*) '     SSS MAX               SSS MIN (unit:t/m/m)'
      WRITE(8,870) (K,(SMAX(K,J),J=1,3),(SMIN(K,J),J=1,3),ME(K),K=1,NE)
      WRITE(8,*)  '  MBYL SSS(n,/,z) MAX (unit: t/m/m) MBYL SSS MIN'

      DO 510 K=1,MB
      L=IP(K)
      WRITE(8,871) K,L,(SSMB1(L,J),J=1,3),(SSMB2(L,J),J=1,3),ME(L)
510      CONTINUE
870      FORMAT(I5,6F15.2,I4)
871      FORMAT(2I5,6F15.2,I4)
      WRITE(8,*) ' K FANG  DU MAX(short,loang,n)  MIN(3) (unit:mm)'
      DO 490 I=1,IHJ
      K=IG(I)
      MN=ME(K)
      T=AE(1,MN)
      IF(ABS(T-0.1).LT.0.0001) GO TO 490
      DO 480 L=1,3
      DUMAX(K,L)=DUMAX(K,L)*10.0 !THE DISPLACMENT OF JOINT UNIT IS MM
      DUMIN(K,L)=DUMIN(K,L)*10.0
480      CONTINUE      
      WRITE(8,880) K,(DUMAX(K,J),J=1,3),(DUMIN(K,J)*(-1),J=1,3),MN
490      CONTINUE
880      FORMAT(I5,6E15.3,I4)
      DO K=1,NE
      WRITE(100,*)K,GMAX(k)
      WRITE(200,890)K,ST1(K),ST2(K),ME(K)
      ENDDO
890      FORMAT(I5,2F10.2,I4)
      DO I=1,NP
      WRITE(61,901)(AGX(I,J),J=1,3)
      WRITE(63,901)(UMAX(I,J)*(100.0),J=1,3)
      ENDDO
      DO K=1,NE
      WRITE(62,902)(SMAX(K,J)*(0.01),J=1,3),(SMin(K,J)*(-0.01),J=1,3)
      ENDDO
901   FORMAT(3F12.5)
902   FORMAT(6F12.5)
      RETURN
      END
c------------------------------------------------------------------
      SUBROUTINE RXY(UGTT,R)
      DIMENSION ugtt(3),R(990000),at(990000),bt(990000)
      COMMON /AB/seita,dt,ab0,ab1,ab2,ab3,ab4,ab5,ab6,ab7,ab8
      COMMON /A1/SK(125930000),DM(990000)
      COMMON /A9/SC(125930000)
      COMMON /A6/MA(990000)
      COMMON /A4/N,NH,MX,JR(3,500000)
      COMMON /u012/u0(990000),u1(990000),u2(990000)
      COMMON /A3/NE,NP,NR
      COMMON/STRU3/ICOL(125930000)
      do 40 i=1,n
      at(i)=ab0*u0(i)+ab2*u1(i)+2.0*u2(i)
      bt(i)=ab1*u0(i)+2.0*u1(i)+ab3*u2(i)
40      continue
        DO 125 I=1,N
        R(I)=0.0
125     CONTINUE       
C       
        R(1) = SC(1)*BT(1)
        DO 100 I=2,N
        I1=MA(I-1)
        I2=MA(I)
        I3=I2-I1
        DO 100 J=1,I3
        J1=I1+J
        J2=ICOL(J1)
        R(I)=R(I)+SC(J1)*BT(J2)
        IF(J2.NE.I) R(J2)=R(J2)+SC(J1)*BT(I)
100     CONTINUE
      do 50 i=1,n
      r(i)=r(i)+at(i)*DM(i)
50      continue
      do 75 i=1,np
      do 76 j=1,3
      iv=jr(j,i)
      if(iv.eq.0) go to 76
      r(iv)=r(iv)-DM(iv)*ugtt(j)
76      continue
75      continue
      END


      SUBROUTINE  GDFIRS
      COMMON /A3/NE,NP,NR
      COMMON /A2/IPE(500000,8)
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /STRM/SSS(500000,3)
      COMMON /D17/SG(500000)
      COMMON /EJOINT/IJKL(500000),VML(500000)
      DO 100 K=1,NE
        III=IJKL(K)
      V=VML(K)
      IF(ABS(V).LT.1.0E-5) GO TO 50
      A=0.0
      DO 30 I=1,3
      IF(SSS(K,I).GT.0.0) A=A+SSS(K,I)
30      CONTINUE
      A=A/3.0
      IF(A.LT.10.0) A=10.0
      SG(K)=A
      GO TO 90
50      SG(K)=SSS(K,3)
      IF(SG(K).LE.10.0) SG(K)=10.0
90      CONTINUE
100      CONTINUE

      DO 60 K=1,NE
      DO 60 J=1,3
      SSS(K,J)=0.0
60      CONTINUE

      RETURN
      END

      SUBROUTINE DAMPCURVE(K,GM,G,DAMPF)
      COMMON /A8/ME(500000)
      COMMON /D14/GMV(500000)
      COMMON /DAMPLY/NEROCK(30),DAMPC(20,50)
      COMMON /METRICAL/MER,NPOINT
      NM=ME(K)
      NJP=NPOINT
      DO 10 I=1,MER
      IF(NM.EQ.NEROCK(I)) L=I
      IF(NM.EQ.NEROCK(I)) GOTO 20
10      CONTINUE
20      CONTINUE
      J1=2*L
      J2=J1+1
      X=GMV(K)
      Y=0.0
      IF(X.LE.DAMPC(1,1)) G=GM*DAMPC(1,J1)
      IF(X.LE.DAMPC(1,1)) DAMPF=DAMPC(1,J2)
      IF(X.LE.DAMPC(1,1)) GOTO 50
      IF(X.GE.DAMPC(NJP,1)) G=GM*DAMPC(NJP,J1)
      IF(X.GE.DAMPC(NJP,1)) DAMPF=DAMPC(NJP,J2)
      IF(X.GE.DAMPC(NJP,1)) GOTO 50
      DO 30 I=1,NJP
      IF(X.GT.DAMPC(I,1).AND.X.LE.DAMPC(I+1,1)) M=I
      IF(X.GT.DAMPC(I,1).AND.X.LE.DAMPC(I+1,1)) GOTO 40
30      CONTINUE
      GOTO 50
40      CONTINUE
      X1=DAMPC(M,1)
      X2=DAMPC(M+1,1)
      X1=LOG10(X1)
      X2=LOG10(X2)
      Y1=DAMPC(M,J1)  
      Y2=DAMPC(M+1,J1)
      Y3=DAMPC(M,J2)  
      Y4=DAMPC(M+1,J2)
      AK=(Y2-Y1)/(X2-X1)
      X=LOG10(X)
      Y=Y1+AK*(X-X1)
      G=GM*Y
      AK=(Y4-Y3)/(X2-X1)
      DAMPF=Y3+AK*(X-X1)
50      CONTINUE
      write(16,100) k,GMV(K),y,gm,g,dampf,nm
100      format(' k,x,y,gam,g,dampf,NM=',i5,2f10.6,2f10.2,f10.4,i4)
      RETURN
      END

C----------------------------------------------
      SUBROUTINE GDMODE
C     CALCULATE THE DYNAMIC MODULI OF ROCKFILL ELEMENTS
      COMMON /D15/IG(50000),DU1(500000,3)
      COMMON /A13/ET(500000),UT(500000),GAMT(500000)
      COMMON /D12/DE(500000),OMGA
      COMMON /A8/ME(500000)
      COMMON /A3/NE,NP,NR
      COMMON /A2/IPE(500000,8)
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /D11/KG,GK(9,40)
      COMMON /D14/GMV(500000)
      COMMON /D16/DL(500000)
      COMMON /D17/SG(500000)
      COMMON /ECON/ETC(500000),UTC(500000)
      COMMON /EJOINT/IJKL(500000),VML(500000)
      IK=0
      DO 150 K=1,NE
      DL(K)=0.0
      GMV(K)=GMV(K)*0.65
      IM=ME(K)
      T=AE(1,IM)
      V=VML(K)
      IF(ABS(V).LT.1.0E-5)    GOTO 90
      IF(AE(2,IM).LT.15000.0) GOTO 30
      E=GK(1,IM)*ETC(K)*0.5
      U=GK(2,IM)
      DE(K)=GK(3,IM)
      GO TO 120
30    CONTINUE
      PE=SG(K)/10.0
      IF(PE.LE.1.0)PE=1.0
      CN=GK(2,IM)
      GM=GK(1,IM)*10.0*((PE)**CN)

      IF(KG.NE.1) GOTO 50

      GO TO 40
50      CONTINUE

      IF(KG.EQ.4) CALL DAMPCURVE(K,GM,G,DAMPF)
      IF(KG.EQ.4) goto 40


40    E=2.0*G*(1.0+UT(K))
      U=UT(K)
      DE(K)=DAMPF
!      IF(DE(K).LE.0.05)DE(K)=0.05
      GO TO 120
90    IK=IK+1
      IF(ABS(T-0.2).LT.1.0E-5) GO TO 100
      conff=1.0
      SKMAX=GK(1,IM)*(SG(K))**0.7*100.0
      TF=SG(K)*TAN(GK(2,IM))
      DELYX=ABS(DU1(K,1))
      DELYZ=ABS(DU1(K,2))
      E=SKMAX/(1.0+conff*DELYX*SKMAX*2.0/TF)
      U=SKMAX/(1.0+conff*DELYZ*SKMAX*2.0/TF)
      YXZ=SQRT(DELYX*DELYX+DELYZ*DELYZ)
      PSK=SKMAX/(1.0+conff*SKMAX*YXZ*2.0/TF)
      DE(K)=GK(3,IM)*(1.1-PSK/SKMAX)
      IF(E.GT.5500) E=5500.0
      IF(U.GT.5500) U=5500.0
      IF(E.LT.20.0) E=20.0
      IF(U.LT.20.0) U=20.0
      GO TO 120
100      CONTINUE
!      X=AE(2,IM)                                       !以下为缝单元
!      Y=AE(3,IM)
      E=ETC(K) !22.50*X
      U=UTC(K) !60.80*X+140.00*Y
      IF(E.GT.5500) E=5500.0
      IF(U.GT.5500) U=5500.0
!      IF(GAMT(K).GT.5500) GAMT(K)=5500.0
      IF(E.LT.20.0) E=20.0
      IF(U.LT.20.0) U=20.0
      IF(GAMT(K).LT.20) GAMT(K)=20.0
!      GAMT(K)=65.00*X+53.00*Y
!      IF(E.GT.ETC(K)) E=ETC(K)
!      IF(U.GT.UTC(K)) U=UTC(K)      
      DE(K)=GK(3,IM)
120   CONTINUE
      DL(K)=ABS((ET(K)-E)/E)
      write(16,170) K,IM,E,U,GAMT(K),DE(K),dl(k)
170      FORMAT('K,IM,E,U,GAMT(K),DE(K)=',2I5,3F12.2,2F18.4)
      ET(K)=E
      UT(K)=U
150      CONTINUE
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!    SUBRUTINRE INPUT   !!!!!!!!!!!!!!!!!!!!!!!!
!    FUNCTION      输入数据
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      SUBROUTINE INPUT(NM)
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /A2/IPE(500000,8)
      COMMON /A8/ME(500000)
      COMMON /A3/NE,NP,NR
      COMMON /CMN4/NEE,NET,NRW,NDE
      COMMON /NJX/NJJ,MBJ,HDMAX,HDAM(100),HHY(100),II0(100),HHY2(100)
      COMMON /A10/JC(500000)
      COMMON /A11/net1(100),nrw1(100),NRWD1(100)
      COMMON /FA/ALFA,HDZ
      COMMON /MODE/MX
      READ(15,*) NDE                    
      READ(15,*) NE,NP,NR
      READ(15,*) (IO,(IPE(I,J),J=1,8),ME(I),I=1,NE)
      READ(15,*) (IO,(COP(I,J),J=1,3),I=1,NP)      
      READ(15,*) (JC(I),I=1,NR)      
      READ(15,*) NJJ,MBJ,HDMAX,(HDAM(I),I=1,NJJ),(HHY(I),I=1,NJJ)
      READ(15,*) (HHY2(I),I=1,NJJ)
      READ(15,*) (II0(I),I=1,NJJ)
      READ(15,*) ALFA,HDZ
      READ(15,*) GAMWA
      READ(15,*) NM,MX,((AE(I,J),I=1,11),J=1,NM)
      READ(15,*) (NET1(I),I=1,NJJ),(NRW1(I),I=1,NJJ),(NRWD1(I),I=1,NJJ)      
      RETURN
      END

      SUBROUTINE WHDUS(K,III,SS,UU)
      DIMENSION XYZ(3,8),U(3)
      DIMENSION BAL(3,3),SS(3),UVW(8,3),UU(3)
      COMMON /CMN51/FUN(8),P(2,8),XJR(2,3)
      COMMON /A2/IPE(500000,8)
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /A5/R(990000)
      COMMON /A4/N,NH,MX,JR(3,500000)
      COMMON /A13/ET(500000),UT(500000),GAMT(500000)
      DO 5 I=1,3
      UU(I)=0.0
      SS(I)=0.0
5      CONTINUE
      DO 10 I=1,III
      IL=IPE(K,I)
      DO 10 J=1,3
      XYZ(J,I)=COP(IL,J)
      IV=JR(J,IL)
      IF(IV.EQ.0) UVW(I,J)=0.0
      IF(IV.NE.0) UVW(I,J)=R(IV)
10      CONTINUE
      IF(III.EQ.8) CALL WHD4FX(XYZ,BAL)
      IF(III.EQ.6) CALL CL3(XYZ,BAL)
      IM=III/2
      DO 80 J=1,3
      U(J)=0.0
      DO 80 I=1,IM
      IC=I+IM
      U(J)=U(J)+UVW(I,J)-UVW(IC,J)
80      CONTINUE
      DO 90 I=1,3
      U(I)=U(I)/FLOAT(IM)
90      CONTINUE
      DO 95 I=1,3
      UU(I)=0.0
      DO 95 J=1,3
      UU(I)=UU(I)+BAL(I,J)*U(J)
95      CONTINUE
      SS(1)=ET(K)*UU(1)
      SS(2)=UT(K)*UU(2)
      SS(3)=GAMT(K)*UU(3)
      RETURN
      END
C------------------------------------------------------
      SUBROUTINE  KRGAM(KK)
      DIMENSION XYZ(3,8)
      COMMON /A1/SK(125930000),SM(990000)
      COMMON /A9/SC(125930000)
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /A6/MA(990000)
      COMMON /A2/IPE(500000,8)
      COMMON /A3/NE,NP,NR
      COMMON /A4/N,NH,MX,JR(3,500000)
      COMMON /A5/R(990000)
      COMMON /CMN4/K,NET,NRW,NDE
      COMMON /GAUSS/RSTG(3),H(3)
      COMMON /CMN3/SKE(24,24),RF(24),CE(24)
      COMMON /NJX/NJJ,MBJ,HDMAX,HD(100),HY(100),II0(100),HHY2(100)
      COMMON /A13/ET(500000),UT(500000),GAMT(500000)
      COMMON /A8/ME(500000)
      COMMON /D12/DE(500000),OMGA
      COMMON /EJOINT/IJKL(500000),VML(500000)
      COMMON/STRU3/ICOL(125930000)
      IF(KK.GT.1) JK=II0(KK-1)
      IF(KK.EQ.1) JK=0
      DO 10 I=1,NH
      SK(I)=0.0
      SC(I)=0.0
10      CONTINUE
      DO 15 I=1,N
      SM(I)=0.0
15      CONTINUE
      KE=II0(KK)
      DO 20 K=1,KE
      NM=ME(K)
      AIX=AE(1,NM)
      DO 30 I=1,8
      IV=IPE(K,I)
      IF(IV.EQ.0) GO TO 30
      DO 28 J=1,3
      XYZ(J,I)=COP(IV,J)
28      CONTINUE
30      CONTINUE
      III=IJKL(K)
      V1=VML(K)
      E=ET(K)
      U=UT(K)
      GAM=GAMT(K)
      IF(III.EQ.4) CALL KRS(K,JK,GAM)
      IF(V1.LE.1.0E-5.AND.III.EQ.6) CALL SJXK(E,U,GAM,XYZ)
      IF(V1.LE.1.0E-5.AND.III.EQ.8.AND.ABS(AIX-0.1).LT.1.0E-5)
     *             CALL WHDMK(E,U,GAM,XYZ)
      IF(V1.LE.1.0E-5.AND.III.EQ.8.AND.ABS(AIX-0.2).LT.1.0E-5)
     *             CALL WHDFK(E,U,GAM,XYZ)
      IF(V1.GT.1.0E-5.AND.III.EQ.8) CALL STIF(K,JK,E,U,GAM,XYZ,III)
      IF(V1.GT.1.0E-5.AND.III.EQ.6) CALL STIF6(K,JK,E,U,GAM,XYZ,III)
        DO 45 II=1,8
        IV=IPE(K,II)
        IF(IV.EQ.0) GO TO 45

        IQ1=(II-1)*3
        DO 40 IJ=1,3
        IU=JR(IJ,IV)
        IP1=IQ1+IJ
        IF(IU.EQ.0) GOTO 40
        JV=MA(IU)
        DO 55 JI=1,8
        IW=IPE(K,JI)
        IF(IW.EQ.0) GO TO 55
        IQ2=(JI-1)*3
        DO 50 JJ=1,3
        IP2=IQ2+JJ
        IX=JR(JJ,IW)
        IF(IX.EQ.0)  GOTO 50
        IF(IX.GT.IU) GOTO 50              !只存放下三角
        if(iu.gt.1) il=MA(iu)-MA(iu-1)
        if(iu.eq.1) il=MA(1)
        do 25 k3=1,il
        if(iu.gt.1) k1=MA(iu-1)+k3
        if(iu.eq.1) k1=k3
        k2=ICOL(k1)
        if(k2.eq.ix) then
        sk(k1)=sk(k1)+ske(ip1,ip2)
        SC(K1)=SC(K1)+SKE(IP1,IP2)*DE(K)/OMGA
        goto 26
        endif
25      continue
26      continue
50        CONTINUE
55        CONTINUE
        SM(IU)=SM(IU)+ABS(CE(IP1))
         SC(JV)=SC(JV)+ABS(CE(IP1))*DE(K)*OMGA 
60      CONTINUE
          IF(NET.EQ.1) R(IU)=R(IU)+RF(IP1)  !Gravity
40        CONTINUE
45        CONTINUE
20      CONTINUE
      CALL FJZL
      RETURN
      END


      SUBROUTINE WHDMK(E,U,GAM,XYZ)
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
10      CONTINUE
      A(2,1)=2.0
      A(3,1)=2.0
      A(3,2)=1.0
      A(4,1)=1.0
      A(4,2)=2.0
      A(4,3)=2.0
      DO 20 I=1,3
      DO 20 J=I+1,4
       A(I,J)=A(J,I)
20      CONTINUE
C      IF(NEE.EQ.49) WRITE(16,150) B,AL,S,((A(I,J),J=1,4),I=1,4)
150      FORMAT(' B,L,S=',3F15.3/'  A(4,4)='/4F10.2/4F10.2/4F10.2/4F10.2)
      DO 30 I=1,4
      DO 30 J=1,4
      A(I,J)=A(I,J)*S/36.0
30      CONTINUE
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
C      IF(NEE.EQ.49) WRITE(16,180) ((A(I,J),J=1,8),I=1,8)
180      FORMAT(8F10.2)
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
90      CONTINUE
100      CONTINUE
110      CONTINUE
      RETURN
      END

      SUBROUTINE WHDFK(E,U,GAM,XYZ)
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
90      CONTINUE
100      CONTINUE
110      CONTINUE
      RETURN
      END

      SUBROUTINE KRS(K,JK,GAM)
      COMMON /CMN4/NEE,NET,NRW,NDE
      COMMON /CMN3/SKE(24,24),RF(24),CE(24)
      COMMON /A13/ET(500000),UT(500000),GAMT(500000)
      COMMON/B11/BI(4),CI(4),DI(4)
      DO 5 I=1,24
      RF(I)=0.0
      CE(I)=0.0
      DO 5 J=1,24
      SKE(I,J)=0.0
5      CONTINUE
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
      V1=UT(K)
      A1=V1/(1.0-V1)
      A2=(1.0-2.0*V1)/(1.0-V1)/2.0
      T1=(1.0-V1)*ET(K)/(1.0-2.0*V1)/V/36.0/(1.0+V1)
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
35      CONTINUE
      IF(NET.NE.1) GO TO 60
      IF(K.LE.JK) GO TO 60
      DO 40 I=1,4
      J=I*3-1
      RF(J)=0.25*W
40      CONTINUE
60      CONTINUE
      RETURN
      END
C-------------------------------------------------
      SUBROUTINE FOUR(K,V)
      DIMENSION X(4),Y(4),Z(4)
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /A2/IPE(500000,8)
      COMMON/B11/BI(4),CI(4),DI(4)
      DO 10 I=1,4
      KP=IPE(K,I)
      X(I)=COP(KP,1)
      Y(I)=COP(KP,2)
      Z(I)=COP(KP,3)
10      CONTINUE
      DO 20 J=1,4
      IF(J+1-4) 30, 30, 40
30    L=J+1
      GO TO 50
40    L=1
50    CONTINUE
      IF(J+2-4)60, 60, 70
60    M=J+2
      GO TO 80
70    M=J-2
80    CONTINUE
      IF(J+3-4)90, 90, 100
90    I=J+3
      GO TO 200
100   I=J-1
200   BI(J)=Y(L)*(Z(I)-Z(M))+Y(M)*(Z(L)-Z(I))+Y(I)*(Z(M)-Z(L))
      CI(J)=X(L)*(Z(M)-Z(I))+X(M)*(Z(I)-Z(L))+X(I)*(Z(L)-Z(M))
      DI(J)=X(L)*(Y(I)-Y(M))+X(M)*(Y(L)-Y(I))+X(I)*(Y(M)-Y(L))
20    CONTINUE
      V=VOLM(X,Y,Z)
      RETURN
      END
C-----------------------------------------------------------
      SUBROUTINE STIF(KA,JK,E,U,GAM,XYZ,III)
      DIMENSION XYZ(3,8),RJX(3,3),Q(3,8),BV(24),D(9)
      COMMON /CMN5/FUN(8),P(3,8),XJR(3,3)
      COMMON /CMN4/NEE,NET,NRW,NDE
      COMMON /GAUSS/RSTG(3),H(3)
      COMMON /CMN3/SKE(24,24),RF(24),CE(24)
      D1=E*(1.0-U)/((1.0+U)*(1.0-2.0*U))
      D2=E*U/((1.0+U)*(1.0-2.0*U))
      D3=E*0.5/(1.0+U)
      DO 5 I=1,24
      RF(I)=0.0
      CE(I)=0.0
      DO 5 J=1,24
      SKE(I,J)=0.0
5      CONTINUE
      DO 100 IT=1,3
      T=RSTG(IT)
      TH=H(IT)
      DO 90 IS=1,3
      S=RSTG(IS)
      SH=H(IS)
      DO 80 IR=1,3
      R=RSTG(IR)
      RH=H(IR)
      CALL RMSD(XYZ,Q,DET,R,S,T,RJX,III)
      DO 35 I=1,8
      LJ=(I-1)*3
      ALS=-FUN(I)*RH*SH*TH*DET*GAM
      DO 35 J=1,3
      IJ=LJ+J
      CE(IJ)=CE(IJ)+ALS/9.80665
35      CONTINUE
20      IF(NET.EQ.0) GO TO 40
      IF(KA.LE.JK) GO TO 40
      DO 30 I=1,8
      J=I*3-1
      RF(J)=RF(J)-FUN(I)*RH*SH*TH*DET*GAM
30      CONTINUE
40      K3=0
      DO 10 I=1,8
      K3=K3+3
      K2=K3-1
      K1=K2-1
      BV(K1)=Q(1,I)
      BV(K2)=Q(2,I)
      BV(K3)=Q(3,I)
10      CONTINUE
      DO 70 I=1,24
      DO 60 J=I,24
      SKE(I,J)=SKE(I,J)+BV(I)*BV(J)*DET*RH*SH*TH
60      CONTINUE
70      CONTINUE
80      CONTINUE
90      CONTINUE
100      CONTINUE
      DO 110 I=2,24
      M=I-1
      DO 110 J=1,M
      SKE(I,J)=SKE(J,I)
110      CONTINUE
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
115      CONTINUE
120      CONTINUE
      SKE(K1,L1)=D(1)*D1+(D(5)+D(9))*D3
      SKE(K2,L2)=D(5)*D1+(D(1)+D(9))*D3
      SKE(K3,L3)=D(9)*D1+(D(5)+D(1))*D3
      SKE(K1,L2)=D(2)*D2+D(4)*D3
      SKE(K2,L1)=D(4)*D2+D(2)*D3
      SKE(K2,L3)=D(6)*D2+D(8)*D3
      SKE(K3,L2)=D(8)*D2+D(6)*D3
      SKE(K3,L1)=D(7)*D2+D(3)*D3
      SKE(K1,L3)=D(3)*D2+D(7)*D3
125      CONTINUE
130      CONTINUE
      RETURN
      END
      SUBROUTINE SJXK(E,U,GAM,XYZ)
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
10    CONTINUE
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
90    CONTINUE
100   CONTINUE
110   CONTINUE
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE CL3(XYZ,BAL)
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
C----------------------------------------------------------------------
      SUBROUTINE  CH29(BAL,BAK,BAD)
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
35    CONTINUE
30    CONTINUE
      DO 40 I=1,3
      DO 40 J=1,3
      BAD(I,J)=0.0
      DO 42 II=1,3
      BAD(I,J)=BAD(I,J)+BAC(I,II)*BAL(II,J)
42    CONTINUE
40    CONTINUE
      RETURN
      END

      SUBROUTINE PFD(R,S,AAB,XYZ)
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
80      E(I)=0.0
      DO 90 I=1,3
      E(1)=E(1)+XJR(1,I)*XJR(1,I)
      E(2)=E(2)+XJR(2,I)*XJR(2,I)
      E(3)=E(3)+XJR(1,I)*XJR(2,I)
90      CONTINUE
      DET=E(1)*E(2)-E(3)*E(3)
      AAB=SQRT(DET)
      RETURN
      END
CC---------------------------------------------------------------------
      FUNCTION      VOLM(X,Y,Z)
C     IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(4),Y(4),Z(4)
      VOLM=((X(2)-X(1))*((Y(3)
     1-Y(1))*(Z(4)-Z(1))-(Y(4)-Y(1))*(Z(3)-Z(1)))
     1-(Y(2)-Y(1))*((X(3)-X(1))*(Z(4)-Z(1))-(X(4)-X(1))*(Z(3)-Z(1)))+
     1(Z(2)-Z(1))*((X(3)-X(1))*(Y(4)-Y(1))-(X(4)-X(1))*(Y(3)-Y(1))))
     1 /6.0
C     WRITE(16,6) (X(IO),Y(IO),Z(IO),IO=1,4)
6     FORMAT(2X,'X=',3E12.5)
      RETURN
      END
      SUBROUTINE V68(XYZ,III,V)
      DIMENSION XYZ(3,8),X(4),Y(4),Z(4),NOD(4,6),NOD1(4,3)
      DATA NOD/1,4,2,6,
     *             1,5,4,6,
     *             1,3,4,5,
     *             1,4,8,2,
     *             1,8,4,7,
     *             1,4,3,7/
      DATA NOD1/1,4,5,6,
     *              1,5,2,6,
     *              1,2,3,6/
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
100      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE WHD4FX(XYZ,BAL)
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
10      BAL(2,1)=BL
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
      SUBROUTINE  CH30(BAL,BAK,BAD)
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
35    CONTINUE
30    CONTINUE
      DO 40 I=1,3
      DO 40 J=1,3
      BAD(I,J)=0.0
      DO 42 II=1,3
      BAD(I,J)=BAD(I,J)+BAC(I,II)*BAL(II,J)
42    CONTINUE
40    CONTINUE
      RETURN
      END

      SUBROUTINE MBYL(ID)
      DIMENSION BAL(3,3),BAK(3,3),BAD(3,3),BLA(3,3),SZ(6)
      COMMON /A15/STRZ(500000,6)
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /A8/ME(500000)
      COMMON /A3/NE,NP,NR
      COMMON /ST/STR(500000,6),EPG(500000,6)
      COMMON /MBJOINT/MB,IP(50000)
      IF(ID.EQ.0) WRITE(16,*) ' K     STR(n,/,z)   (unit:t/m/m)   NM='
      IF(ID.EQ.0) WRITE(16,*)'法向  顺坡向  坝轴向'
      CALL MBFX(BAL)
      DO 30 I=1,3
      DO 30 J=1,3
      BLA(I,J)=BAL(J,I)
30      CONTINUE
      DO 40 I=1,3
      DO 40 J=1,3
      BAL(I,J)=BLA(I,J)
40      CONTINUE
      DO 60 IK=1,MB
      K=IP(IK)
      DO 10 J=1,6
      IF(ID.EQ.0) SZ(J)=STRZ(K,J)
      IF(ID.EQ.1) SZ(J)=STR(K,J)
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
      CALL CH30(BAL,BAK,BAD)
      STR(K,1)=BAD(1,1)
      STR(K,2)=BAD(2,2)
      STR(K,3)=BAD(3,3)
60      CONTINUE
      RETURN
      END

      SUBROUTINE MBFX(BAL)
      DIMENSION BAL(3,3)
      COMMON /FA/AN,HDZ

      DO 10 I=1,3
      DO 10 J=1,3
      BAL(I,J)=0.0
10      CONTINUE
      BAL(3,3)=1.0
      T=SQRT(AN*AN+1.0)
      BAL(1,1)=1.0/T
      BAL(1,2)=-AN/T
      BAL(2,1)=AN/T
      BAL(2,2)=1.0/T
      RETURN
      END

C-----------------------------------------------------------
      SUBROUTINE STIF6(KA,JK,E,U,GAM,XYZ,III)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
5      CONTINUE
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
60      CONTINUE
70      CONTINUE
80      CONTINUE
90      CONTINUE
100      CONTINUE
      DO 110 I=2,18
      M=I-1
      DO 110 J=1,M
      SKE(I,J)=SKE(J,I)
110      CONTINUE
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
115      CONTINUE
120      CONTINUE
      SKE(K1,L1)=D(1)*D1+(D(5)+D(9))*D3
      SKE(K2,L2)=D(5)*D1+(D(1)+D(9))*D3
      SKE(K3,L3)=D(9)*D1+(D(5)+D(1))*D3
      SKE(K1,L2)=D(2)*D2+D(4)*D3
      SKE(K2,L1)=D(4)*D2+D(2)*D3
      SKE(K2,L3)=D(6)*D2+D(8)*D3
      SKE(K3,L2)=D(8)*D2+D(6)*D3
      SKE(K3,L1)=D(7)*D2+D(3)*D3
      SKE(K1,L3)=D(3)*D2+D(7)*D3
125      CONTINUE
130      CONTINUE
      RETURN
      END


      
C-------------------------------------------------------
      SUBROUTINE STRESS(KK)
      DIMENSION XYZ(3,8),RJX(3,3),Q(3,8)
      DIMENSION UVW(3,8),X(9),C(6),MM(4),BP(6,3)
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /ST/STR(500000,6),EPG(500000,6)
      COMMON /A2/IPE(500000,8)
      COMMON /A8/ME(500000)
      COMMON /A3/NE,NP,NR
      COMMON /A5/R(990000)
      COMMON /CMN4/NEE,NET,NRW,NDE
      COMMON /A4/N,NH,MX,JR(3,500000)
      COMMON /NJX/NJJ,MBJ,HDMAX,HDAM(100),HHY(100),II0(100),HHY2(100)
      COMMON /A13/ET(500000),UT(500000),GAMT(500000)
      COMMON /B11/BI(4),CI(4),DI(4)
      COMMON /EJOINT/IJKL(500000),VML(500000)
      DATA  MM/1.0,-1.0,1.0,-1.0/
      DO 10 I=1,NE
      DO 10 J=1,6
      STR(I,J)=0.0
      EPG(I,J)=0.0
10      CONTINUE
      KE=II0(KK)
      DO 150 K=1,KE
      DO 20 I=1,8
      IV=IPE(K,I)
      DO 20 J=1,3
      XYZ(J,I)=COP(IV,J)
20      CONTINUE
      III=IJKL(K)
      V=VML(K)
      V1=ABS(V)
      IF(V1.LT.1.0E-5) GO TO 100
      E=ET(K)
      U=UT(K)
      D1=E*(1.0-U)/((1.0+U)*(1.0-2.0*U))
      D2=E*U/((1.0+U)*(1.0-2.0*U))
      D3=E*0.5/(1.0+U)
      DO 30 I=1,8
      IV=IPE(K,I)
      IF(IV.EQ.0) GO TO 30
      DO 32 J=1,3
      IW=JR(J,IV)
      IF(IW.EQ.0) UVW(J,I)=0.0
      IF(IW.NE.0) UVW(J,I)=R(IW)
32      CONTINUE
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
      IV=IPE(K,I)
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
40      CONTINUE
      C(1)=D1*X(1)+D2*(X(2)+X(3))
      C(2)=D1*X(2)+D2*(X(3)+X(1))
      C(3)=D1*X(3)+D2*(X(1)+X(2))
      C(4)=D3*(X(7)+X(4))
      C(5)=D3*(X(9)+X(6))
      C(6)=D3*(X(8)+X(5))
      DO 80 I=1,6
      STR(K,I)=-C(I)
80      CONTINUE
      DO 90 I=1,3
90      EPG(K,I)=-X(I)
      EPG(K,4)=-(X(7)+X(4))
      EPG(K,5)=-(X(9)+X(6))
      EPG(K,6)=-(X(8)+X(5))
      GO TO 100
200   CALL FOUR(K,V)
      DO 210 J=1,6
210   X(J)=0.0
      DO 220 I=1,3
      DO 220 J=1,6
220   BP(J,I)=0.0
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
230      CONTINUE
      C(1)=D1*X(1)+D2*(X(2)+X(3))
      C(2)=D1*X(2)+D2*(X(3)+X(1))
      C(3)=D1*X(3)+D2*(X(1)+X(2))
      C(4)=D3*X(4)
      C(5)=D3*X(5)
      C(6)=D3*X(6)
      DO 240 I=1,6
      STR(K,I)=-C(I)
240      EPG(K,I)=-X(I)
100      CONTINUE
150      CONTINUE
      RETURN
      END

C-------------------------------------------------------------
      


C-----------------------------------------------------
      SUBROUTINE MAIN(IM)
      DIMENSION SW(6),PS(3)
      COMMON /ST/STR(500000,6),EPG(500000,6)
      COMMON /A3/NE,NP,NR
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /STRM/SSS(500000,3)
      COMMON /A2/IPE(500000,8)
      COMMON /D13/GMV(500000)
      COMMON /EJOINT/IJKL(500000),VML(500000)
      COMMON /FACTOR11/ACCEFACTOR,STRAINFACTOR,EPS2
      DO 100 K=1,NE
      III=IJKL(K)
      V=VML(K)
      IF(ABS(V).LE.1.0E-5) GO TO 100
      DO 30 I=1,3
      PS(I)=0.0
30      CONTINUE
      IF(IM.EQ.1) GO TO 60
      DO 40 I=1,3
40      SW(I)=EPG(K,I)*1.0E6
      DO 50 I=4,6
50      SW(I)=EPG(K,I)*0.5*1.0E6
      GO TO 80
60    DO 70 I=1,6
         SW(I)=STR(K,I)
70      CONTINUE
80      B=-SW(1)-SW(2)-SW(3)
      C=SW(1)*SW(2)+SW(3)*SW(2)+SW(1)*SW(3)-SW(4)*SW(4)-SW(5)*SW(5)
     #-SW(6)*SW(6)
      D=SW(1)*SW(5)*SW(5)+SW(2)*SW(6)*SW(6)+SW(3)*SW(4)*SW(4)
     #-2.0*SW(4)*SW(5)*SW(6)-SW(1)*SW(2)*SW(3)
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
110      S=SQRT(-P/3.0)
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
      S=PS(2)
      PS(2)=PS(3)
      PS(3)=S
120      CONTINUE
      IF(IM.NE.1) GO TO 150
      DO 140 I=1,3
      SSS(K,I)=PS(I)
140      CONTINUE
      GO TO 100
150      GMV(K)=((PS(1)-PS(3))*1.0E-6)*STRAINFACTOR
C      IF(GMV(K).GT.0.1) WRITE(16,170) K,GMV(K)
100      CONTINUE
170      FORMAT(I6,E15.6)
      RETURN
      END

C------------------------------------------------
      
C-----------------------------------------------------------
      SUBROUTINE MJZDDH(KK,NO)
      COMMON /NJX/NJJ,MBJ,HDMAX,HDAM(100),HHY(100),II0(100),HHY2(100)
      COMMON /A3/NE,NP,NR
      COMMON /A4/N,NH,MX,JR(3,500000)
      COMMON /A2/IPE(500000,8)
      NEJK=II0(KK)
      NO=0
      DO 20 I=1,NEJK
      DO 10 J=1,8
      IV=IPE(I,J)
      IF(IV.EQ.0) GO TO 10
      IF(IV.GT.NO) NO=IV
10      CONTINUE
20      CONTINUE
      DO 50 I=1,NO
      DO 40 K=1,NEJK
      DO 30 J=1,8
      IV=IPE(K,J)
      IF(IV.EQ.0) GO TO 40
      IF(IV.EQ.I) GO TO 50
30      CONTINUE
40      CONTINUE
      JR(1,I)=0
      JR(2,I)=0
      JR(3,I)=0
50      CONTINUE
      DO 60 I=NO+1,NP
      JR(1,I)=0
      JR(2,I)=0
      JR(3,I)=0
60      CONTINUE
      RETURN
      END
C-----------------------------------------------------------
      SUBROUTINE MR
      COMMON /A10/JC(500000)
      COMMON /A3/NE,NP,NR
      COMMON /A4/N,NH,MX,JR(3,500000)
      COMMON /NJX/NJJ,MBJ,HDMAX,HDAM(100),HHY(100),II0(100),HHY2(100)
      DO 10 I=1,NP
      DO 10 J=1,3
      JR(J,I)=1
10      CONTINUE
      DO 20 I=1,NR
      K=JC(I)
      J=K/1000
      L=(K-J*1000)/100
      M=(K-J*1000-L*100)/10
      N=K-J*1000-L*100-M*10
      JR(1,J)=L
      JR(2,J)=M
      JR(3,J)=N
20      CONTINUE
      CALL MJZDDH(NJJ,NO)
      N=0
      DO 70 I=1,NO
      DO 60 J=1,3
      IF(JR(J,I)) 60,60,50
50      N=N+1
      JR(J,I)=N
60      CONTINUE
70      CONTINUE
      RETURN
      END
C--------------------------------------------------------

C-------------------      ------------------------------
        SUBROUTINE SORT(NA,N)
        DIMENSION NA(159)
        DO 5 I=1,N-1
        DO 5 J=I+1,N
        IF(NA(I).GT.NA(J)) THEN
          T=NA(I)
        NA(I)=NA(J)
        NA(J)=T
        ENDIF
5       ENDDO
C       
        RETURN
        END
C-------------------      ------------------------------
      SUBROUTINE   FPJD(R,S,T,DET,XYZ)
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
 10      CONTINUE
      DO 40 I=1,3
      DO 30 J=1,3
      DET=0.0
      DO 20 K=1,8
      DET=DET+P(I,K)*XYZ(J,K)
 20      CONTINUE
      XJR(I,J)=DET
 30      CONTINUE
 40      CONTINUE
      DET=XJR(1,1)*XJR(2,2)*XJR(3,3)
     *         +XJR(2,1)*XJR(3,2)*XJR(1,3)
     *         +XJR(1,2)*XJR(2,3)*XJR(3,1)
     *         -XJR(1,3)*XJR(2,2)*XJR(3,1)
     *         -XJR(3,2)*XJR(2,3)*XJR(1,1)
     *         -XJR(1,2)*XJR(2,1)*XJR(3,3)
      IF(DET.LT.1.0E-5) GO TO 100
Cc      WRITE(16,500) DET
      RETURN
 500      FORMAT(10X,'DET=***',F15.3)
 100      WRITE(16,600) NEE,R,S,T
      WRITE(16,500) DET
 600      FORMAT(1X,'FPJD ERROR***  ELEMENT(',I5,')',4X,'R=',F10.5,4X,'S=',
     *            F10.3,' T=',F10.5)
      STOP
      END
C C ---------------------- ------------------------------
      SUBROUTINE RMSD(XYZ,Q,DET,R,S,T,RJX,III)
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
 50      CONTINUE
      Q(J,I)=X
 60      CONTINUE
 70      CONTINUE
      RETURN
      END
C-------------------      ------------------------------
      SUBROUTINE FPJD6(R,S,T,DET,XYZ)
      DIMENSION  XYZ(3,8),ZETA(6)
      COMMON /CMN5/FUN(8),P(3,8),XJR(3,3)
      COMMON /CMN4/NEE,NET,NRW,NDE
      DATA ZETA/-1.0,-1.0,-1.0,1.0,1.0,1.0/
      DO 5 I=1,8
      FUN(I)=0.0
      DO 5 J=1,3
      P(J,I)=0.0
5      CONTINUE      
      DO 10 J=1,2
      K=(J-1)*3
      DO 10 L=1,3
      I=K+L
      ZETA0=ZETA(I)*T
      IF(L.EQ.1) FUN(I)=0.5*R*(1.0+ZETA0)
      IF(L.EQ.2) FUN(I)=0.5*S*(1.0+ZETA0)
      IF(L.EQ.3) FUN(I)=0.5*(1.0-R-S)*(1.0+ZETA0)
10      CONTINUE
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
15      CONTINUE
      DO 40 I=1,3
      DO 30 J=1,3
      DET=0.0
      DO 20 K=1,6
      DET=DET+P(I,K)*XYZ(J,K)
 20      CONTINUE
      XJR(I,J)=DET
 30      CONTINUE
 40      CONTINUE
      DET=XJR(1,1)*XJR(2,2)*XJR(3,3)
     *         +XJR(2,1)*XJR(3,2)*XJR(1,3)
     *         +XJR(1,2)*XJR(2,3)*XJR(3,1)
     *         -XJR(1,3)*XJR(2,2)*XJR(3,1)
     *         -XJR(3,2)*XJR(2,3)*XJR(1,1)
     *         -XJR(1,2)*XJR(2,1)*XJR(3,3)
      IF(DET.LT.1.0E-5) GO TO 100
Cc      WRITE(16,500) DET
      RETURN
 500      FORMAT(10X,'DET=***',F15.3)
 100      WRITE(16,600) NEE,R,S,T
      WRITE(16,500) DET
 600      FORMAT(1X,'FPJD6 ERROR***  ELEMENT(',I5,')',4X,'R=',F10.5,4X,'S=',
     *            F10.3,' T=',F10.5)
      STOP
      END

      SUBROUTINE DEIGV(MJ,IMM)
       DIMENSION AL(16),BT(16)
      DIMENSION T0(16)
      COMMON /A1/SK(125930000),DM(990000)
      COMMON /A6/MA(990000)
      COMMON /A3/NE,NP,NR
      COMMON /A4/N,NH,MX,JR(3,500000)
      COMMON /A5/R(990000)
      COMMON /RRR1/RR(990000,4),R5(990000)
      COMMON /RMD/RLD(16),W0(16)
      COMMON /RDE/RS(990000)
      DO 200 J=1,MJ
1      NUM=0
      BB=0.0
      DO 10 I=1,N
      RR(I,J)=RS(I)
10      CONTINUE
12      DO 30 I=1,N
      R(I)=DM(I)*RR(I,J)
30      CONTINUE
45      CALL SSORPCG
60      RLD(J)=R(1)
      DO 70 L=2,N
      IF(ABS(R(L)).GT.ABS(RLD(J))) RLD(J)=R(L)
70      CONTINUE
      DO 80 L=1,N
80      R(L)=R(L)/RLD(J)
      RE=0.0
      R3=0.0
      CALL VMULTI(R4)
      DO 140      I=1,N
      R3=R3+DM(I)*R(I)*R(I)
140      CONTINUE
      RLD(J)=R3/R4
      AA=RLD(J)
      RE=ABS(AA-BB)/AA
      IF(RE.GT.1.0E-9) BB=AA
      DO 160 I=1,N
160      RR(I,J)=R(I)
      NUM=NUM+1
      WRITE(*,*)RE
      IF(FLOAT(num/10).EQ.FLOAT(num)/10.0) WRITE(*,698) num
      IF(ABS(RE).LE.2.0E-5) GO TO 165
      IF(NUM.GT.50) GO TO 165
      GO TO 12
165      continue
      DO 210 I=1,N
210      R(I)=DM(I)*RR(I,J)
      CALL XTX(J,AJ)
      AL(J)=RLD(J)/AJ
      W0(J)=1.0/SQRT(AA)
      AA=W0(J)
      T0(J)=2.0*3.1415926/AA
      BB=T0(J)
      CC=1.0/BB
      WRITE(*,601) J,NUM,RE,AA,BB,CC
      IF(IMM.EQ.1)write(8,600) NUM,RE,AA,BB,CC
200      CONTINUE
         DO 123 I=1,N
123      RS(I)=RR(I,1)
601      FORMAT(1X,'MODE-',I3,' NUM,RE,W0,T0,F=',I3,E15.4,3F12.5)
550      FORMAT(3X,'MODE--',I3)
600      FORMAT(1X,'NUM=*RE=*W0=*T0,F=*',I6,E15.6,3F12.5)
698      format(5x,'num=',i5)
      RETURN
      END
      SUBROUTINE VMULTI(R4)
      COMMON /A4/N,NH,MX,JR(3,500000)
      COMMON /A5/R(990000)
      COMMON /A1/SK(125930000),DM(990000)
      COMMON /A6/MA(990000)
      COMMON/STRU3/ICOL(125930000)
      DIMENSION A(990000)


             R4=0.0     
        DO 125 I=1,N
        A(I)=0.0
125     CONTINUE       
C       
        A(1) = SK(1)*R(1)
        DO 100 I=2,N
        I1=MA(I-1)
        I2=MA(I)
        I3=I2-I1
        DO 100 J=1,I3
        J1=I1+J
        J2=ICOL(J1)
        A(I)=A(I)+SK(J1)*R(J2)
        IF(J2.NE.I) A(J2)=A(J2)+SK(J1)*R(I)
100     CONTINUE
       DO 130 I=1,N
       R4=R4+R(i)*A(I)
130    CONTINUE
C       
        RETURN
         END
C-------------------------------------------------------------

C--------------------SUB XTX----------------------
      SUBROUTINE XTX(JJ,AJ)
      COMMON /A4/N,NH,MX,JR(3,500000)
      COMMON /A5/R(990000)
      COMMON /RRR1/RR(990000,4),R5(990000)
      AJ=0.0
      DO 10 I=1,N
      AJ=AJ+RR(I,JJ)*R(I)
10      CONTINUE
      RETURN
      END

      SUBROUTINE MWI  ! 参见土石坝地震工程
      DIMENSION XYZ(3,8)
      DIMENSION KFACE(4,6),NJ(4),KFACE6(3,2)
      COMMON /ELE/IPEYW(500000,8)
      COMMON /RWM/RMW(990000)
      COMMON /A3/NE,NP,NR
      COMMON /FA/ALFA,HDZ
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /A8/ME(500000)
      COMMON /NJX/NJJ,MBJ,HDMAX,HDAM(100),HHY(100),II0(100),HHY2(100)
      COMMON /SLABADDMASS/IPSLAB(50000),NUH,MP(50000)
      DATA  KFACE/2,4,6,8,
     *            1,3,5,7,
     *            3,4,7,8,
     *            1,2,5,6,
     *            5,6,7,8,
     *            1,2,3,4/
      DATA KFACE6/4,5,6,1,2,3/
       
      Y0=HHY(NJJ)
      FA=ATAN(1.0/ALFA)
      WRITE(*,*) ' N,ALFA=',ALFA,FA

      DO 10 I=1,NP
      RMW(I)=0.0
10      CONTINUE

      DO 100 IW=1,NUH
      NEE=MP(IW)
      IM=ME(NEE)
      IF(AE(1,IM).LT.2.4.OR.AE(9,IM).NE.0)GOTO 100 !仅面板单元计算附加质量
      ND=IPSLAB(IW)       
      III=8
      IF(IPEYW(NEE,7).EQ.0) III=6
      DO 5 L=1,8
      DO 5 J=1,3
5      XYZ(J,L)=0.0
      DO 20 J=1,8
      IV=IPEYW(NEE,J)
      IF(IV.EQ.0) GO TO 20
      DO 25 L=1,3
      XYZ(L,J)=COP(IV,L)
25      CONTINUE
20      CONTINUE
      IF(III.EQ.8) CALL SA(ND,A,XYZ)
      IF(III.EQ.6) CALL SA3(ND,A,XYZ)
      II=III/2
      DO 30 L=1,II
      IF(III.EQ.8) J=KFACE(L,ND)
      IF(III.EQ.6) J=KFACE6(L,ND)
      IF(J.EQ.0) GO TO 30
      NJ(L)=IPEYW(NEE,J)
30      CONTINUE
      DO 50 L=1,III/2      
      J=NJ(L)
      IF(L.EQ.1) GO TO 56
      DO 55 M=1,L-1
      IF(NJ(M).EQ.J) GO TO 50
55      CONTINUE
56      CONTINUE
      YI=Y0-COP(J,2)
      IF(YI.GT.0.0) RMW(J)=RMW(J)+A/II
50      CONTINUE
60      CONTINUE
100      CONTINUE
      DO 70 I=1,NP
      IF(RMW(I).EQ.0.0) GOTO 70
      CALL HMAX(I,HB)
      H0=Y0-HB
      HI=Y0-COP(I,2)
      IF(HI.LE.0.0) GOTO 70
      RMW(I)=7.0/8.0*RMW(I)*SQRT(H0*HI)/9.80665*FA/1.570796327
70      CONTINUE
      DO 80 I=1,NP
      IF(ABS(RMW(I)).GE.1.0E-2) WRITE(16,110) I,RMW(I) 
80      CONTINUE
110      FORMAT(' DYNAMIC WATER PRESSURE ADD MASS=',I6,F12.3)
120      FORMAT(' NP=',I6,' Z  Y0,YI=',3F10.2)
      RETURN
      END

C     DETERMINE THE BOTTON LEVEL OF THE SECTION OF NODAL i 
      SUBROUTINE HMAX(IK,YMIN)
      DIMENSION KFACE(4,6),KFACE6(3,2)
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /ELE/IPEYW(500000,8)
      COMMON /A3/NE,NP,NR
      COMMON /A8/ME(500000)
      COMMON /EJOINT/IJKL(500000),VML(500000)
      COMMON /SLABADDMASS/IPSLAB(50000),NUH,MP(50000)
      DATA  KFACE/2,4,6,8,
     *            1,3,5,7,
     *            3,4,7,8,
     *            1,2,5,6,
     *            5,6,7,8,
     *            1,2,3,4/
      DATA KFACE6/4,5,6,1,2,3/
      YMIN=15000.0 !      elevation of the dam top
      ZI=COP(IK,3)
      DO 100 IW=1,NUH
      NEE=MP(IW)
      IM=ME(NEE)

      III=0
      DO 20 I=1,8
      IV=IPEYW(NEE,I)                                                                                                                                    
      IF(IV.EQ.0) GO TO 20
      III=III+1                                                                                                                                                      
20      CONTINUE

      IF(AE(1,IM).LT.2.4.OR.AE(9,IM).NE.0)GOTO 100
      ND=IPSLAB(IW)      
      II=III/2
      DO 30 L=1,II
      IF(III.EQ.8) J=KFACE(L,ND)
      IF(III.EQ.6) J=KFACE6(L,ND)
      IF(J.EQ.0) GO TO 30
      NJ=IPEYW(NEE,J)
      IF(ABS(ZI-COP(NJ,3)).GT.1.0E-01) GOTO 30
      IF(YMIN.GT.COP(NJ,2)) YMIN=COP(NJ,2)
30      CONTINUE
100      CONTINUE
      RETURN      
      END
C------------------------------------------------
      SUBROUTINE Sa(ND,A,XYZ)
      DIMENSION  XYZ(3,8),KCRD(6),FVAL(6),IPRM(3),RST(3)
      COMMON /GAUSS/RSTG(3),H(3)
      COMMON /CMN5/FUN(8),P(3,8),XJR(3,3)
      DATA KCRD/1,1,2,2,3,3/
      DATA IPRM/2,3,1/
      DATA FVAL/1.0,-1.0,1.0,-1.0,1.0,-1.0/
      A=0.0
      ML=KCRD(ND)
      MM=IPRM(ML)
      MN=IPRM(MM)
      RST(ML)=FVAL(ND)
      DO 20 LX=1,3
      RST(MM)=RSTG(LX)
      DO 20 LY=1,3
      RST(MN)=RSTG(LY)
      CALL FPJD(RST(1),RST(2),RST(3),DET,XYZ)
      E1=0.0
      E2=0.0
      E3=0.0
      DO 10 K=1,3
      E1=E1+XJR(MM,K)*XJR(MM,K)
      E2=E2+XJR(MN,K)*XJR(MN,K)
      E3=E3+XJR(MM,K)*XJR(MN,K)
10      CONTINUE
       E=SQRT(E1*E2-E3*E3)
      A=A+E*H(LX)*H(LY)
20      CONTINUE
      RETURN
      END 
C------------------------------------------------
      SUBROUTINE Sa3(ND,A,XYZ)
      DIMENSION  XYZ(3,8),FVAL(2)
      COMMON /GAUSS/RSTG(3),H(3)
      COMMON /CMN5/FUN(8),P(3,8),XJR(3,3)
      DATA  FVAL/1.0,-1.0/
      R=1/3
      S=R
      T=FVAL(ND)
      CALL FPJD6(R,S,T,DET,XYZ)
      E1=0.0
      E2=0.0
      E3=0.0
      DO 10 K=1,3
      E1=E1+XJR(1,K)*XJR(1,K)
      E2=E2+XJR(2,K)*XJR(2,K)
      E3=E3+XJR(1,K)*XJR(2,K)
10      CONTINUE
       E=SQRT(E1*E2-E3*E3)
      A=E*0.5
20      CONTINUE
      RETURN
      END 

      SUBROUTINE FJZL
      COMMON /RWM/RMW(990000)
      COMMON /A4/NJ,NH,MX,JR(3,500000)
      COMMON /A3/NE,NP,NR
      COMMON /A1/SK(125930000),DM(990000)
      DO 10 I=1,NP
      IF(RMW(I).EQ.0.0) GOTO 10
      DO 30 J=1,3
      IV=JR(J,I)
      IF(IV.EQ.0) GO TO 30
      DM(IV)=DM(IV)+RMW(I)
30      CONTINUE
10      CONTINUE
      RETURN
      END

      !!!!!!!!!!!!!!!!!!!!!!!!1
      SUBROUTINE SOLVEVOLUME
      DIMENSION XYZ(3,8)
      COMMON /EJOINT/IJKL(500000),VML(500000)
      COMMON /A7/COP(500000,3),AE(11,40)
      COMMON /A2/IPE(500000,8)
      COMMON /A3/NE,NP,NR
      COMMON /A8/ME(500000)

      DO 100 K=1,NE
      IM=ME(K)
      VML(K)=0.0
50      III=0
      DO 20 I=1,8
      IV=IPE(K,I)                                                                                                                                    
      IF(IV.EQ.0) GO TO 20
      III=III+1
      DO 30 J=1,3
      XYZ(J,I)=COP(IV,J)
30      CONTINUE                                                                                                                                                      
20      CONTINUE
40      IJKL(K)=III
      IF(AE(1,IM).LE.0.3.AND.AE(1,IM).GT.0.01)GOTO 101
      CALL V68(XYZ,III,V)
      V1=ABS(V)
      VML(K)=V1
!      WRITE(8,*)k,im,v1,AE(1,IM)
101      CONTINUE
100      CONTINUE
      RETURN
      END
      SUBROUTINE SSORPCG
C** IN THIS SUBROUTINE, THE STORAGE MODE OF THE MATRIX A IS AS FOLLOWS:
C** 1. THE LOW TRIANGLE MATRIX IS STORED;
C** 2. THE NONZERO ELEMENTS ARE STORED;
C** 3. THE ELEMENTS ARE STORED ACCORDING THE ROW;
      COMMON /A1/SK(125930000),DM(990000)
      COMMON /A6/MA(990000)
      COMMON /A4/N,NH,MX,JR(3,500000)
      COMMON /STRU3/ICOL(125930000)
      COMMON /A5/R(990000)

      DIMENSION C(990000),D(990000),E(990000),X(990000)
      DIMENSION G(990000),Y(990000),Z(990000) 
        
      OMIGA=1.2
      OMIGA1=(2-OMIGA)/OMIGA
      EPS=1E-2
      RLAMDA=1E3
      !DYNAMIC ANALYSIS
 
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
C** for the first diagonal, the coresponding nonzero non-diagonal element 
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
15       CONTINUE
       Y(I)=G(I)*OMIGA/SK(MA(I))
C** Z(I)=-VY
       Z(I)=-OMIGA1*(SK(MA(I)))*Y(I)
10      CONTINUE      
C** COMPUTE THE VECTOR D=W(-T)Z
      DO 20 I=N,1,-1
       IF(I.EQ.1) THEN
C** for the first diagonal, the coresponding nonzero non-diagonal element 
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
25     CONTINUE
20      CONTINUE
C** CHECK THE CONVERGENCE
60      DELTA=0.0
        
      DO I=1,N
       DELTA=DELTA+Y(I)**2*OMIGA1*(SK(MA(I)))
      ENDDO
      IF(DELTA.LE.EPS) THEN
       DO I=1,N
        R(I)=X(I)/RLAMDA
       ENDDO
!        WRITE(*,*)'END ITERATION STEP  ', KIter,'TORLERANCE  ',DELTA1
       RETURN
      ENDIF
C** SET THE TEMPERORY VARIBLE ZER
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
C** for the first diagonal, the coresponding nonzero non-diagonal element 
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
45       CONTINUE
       E(I)=C(I)*OMIGA/(SK(MA(I)))
       Y(I)=Y(I)+T*(D(I)+E(I))
       DELTA1=DELTA1+Y(I)**2*OMIGA1*(SK(MA(I)))
40      CONTINUE
        
      BETA0=DELTA1/DELTA
      DO I=1,N
       Z(I)=-OMIGA1*(SK(MA(I)))*Y(I)+BETA0*Z(I)
       D(I)=0.0
      ENDDO
C** D(K+1)=W(-T)Z(K+1)
      DO 50 I=N,1,-1
       IF(I.EQ.1) THEN
C** for the first diagonal, the coresponding nonzero non-diagonal element 
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
55     CONTINUE
50      CONTINUE
      KIter=KIter+1
      IF(MOD(KIter,50).EQ.0) THEN
!       WRITE(*,*)'ITERATION STEP IS ', KIter,'TORLERANCE IS ',DELTA1
      ENDIF
      GOTO 60
C       
      RETURN
      END
        
