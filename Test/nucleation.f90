      program nucleation; implicit none
      integer,parameter :: NKR = 33
      INTEGER NDROPMAX,IDROP,ICCN,INEXT,ISMALL,KR,NCRITI
      INTEGER ICEMAX,IMIN,IMAX,I,II,I0,I1
      REAL &
     &  SUP1,TT,RACTMAX,XKOE,R03,SUPCRITI,AKOE23,RCRITI,BKOE, &
     &  AKOE,CONCCCNIN,DEG01,ALN_IP
      REAL :: CCNCONC(NKR)
      REAL CCNCONC_BFNUCL

      REAL :: COL = 0.23105
      REAL :: RCCN(NKR),DROPRADII(NKR),FCCNR(NKR),FCCNR_old(NKR)
      REAL :: RACT(NKR),DROPCONC(NKR),DROPCONCN(NKR),temp(NKR)
      REAL DLN1,DLN2,FOLD_IP
      
      TT = 273.0!273.15+15.0
      SUP1 = 0.008 !for maritime in WDM6
      FCCNR(:) = 0.0
      do KR = 1,NKR
         read(14,"(2E13.5)") FCCNR(KR), RCCN(KR)
         !write(*,"(2E13.5)") FCCNR(KR), RCCN(KR)
         FCCNR_old(KR)=FCCNR(KR)
         read(18,"(2E13.5)") temp(KR), DROPRADII(KR)
      end do
      DEG01=1./3.

! calculation initial value of NDROPMAX - maximal number of drop bin
! which is activated

! initial value of NDROPMAX

      NDROPMAX=0

      DO KR=1,NKR
! initialization of bin radii of activated drops
         RACT(KR)=0.
! initialization of aerosol(CCN) bin concentrations
         CCNCONC(KR)=0.
! initialization of drop bin concentrations
         DROPCONCN(KR)=0.
      ENDDO
!       CCNCONC_BFNUCL - concentration of aerosol particles before
!                  nucleation
      CCNCONC_BFNUCL=0.
      DO I=1,NKR
         CCNCONC_BFNUCL=CCNCONC_BFNUCL+FCCNR(I)
      ENDDO

      CCNCONC_BFNUCL=CCNCONC_BFNUCL*COL

      IF(CCNCONC_BFNUCL.EQ.0.) THEN
         print*,"CCNCONC_BFNUCL.EQ.0."
         RETURN
      ELSE
         CALL BOUNDARY(IMIN,IMAX,FCCNR,NKR)
         CALL CRITICAL (AKOE,BKOE,TT,RCRITI,SUP1,DEG01)
         IF(RCRITI.GE.RCCN(IMAX))  RETURN
      END IF
 
      print*,"IMIN,  IMAX,  RCRITI"
      print*,IMIN,IMAX,RCRITI
      print*,"----------------------"
! calculation of CCNCONC(I) - aerosol(CCN) bin concentrations;
!                             I=IMIN,...,IMAX
! determination of NCRITI - number bin in which is located RCRITI
        IF (IMIN.EQ.1)THEN
         CALL CCNIMIN(IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC,COL, &
     &       FCCNR,NKR)
         CALL CCNLOOP(IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC,COL, &
     &       FCCNR,NKR)
        ELSE
         CALL CCNLOOP(IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC,COL, &
     &       FCCNR,NKR)
        END IF
      print*,"NCRITI = ",NCRITI
      print*,"          BIN,    CCNCONC,        FCCNR,       FCCNR_old"
      do KR = 1,NKR
         print*,KR,CCNCONC(KR),FCCNR(KR)*COL,FCCNR_old(KR)*COL
      end do
      print*,"----------------------"
! calculation CCNCONC_AFNUCL - ccn concentration after nucleation

!       CCNCONC_AFNUCL=0.

!       DO I=IMIN,IMAX
!          CCNCONC_AFNUCL=CCNCONC_AFNUCL+FCCNR(I)
!       ENDDO

!       CCNCONC_AFNUCL=CCNCONC_AFNUCL*COL

! calculation DEL_CCNCONC

!       DEL_CCNCONC=CCNCONC_BFNUCL-CCNCONC_AFNUCL
        CALL ACTIVATE(IMIN,IMAX,AKOE,BKOE,RCCN,RACT,RACTMAX,NKR)
      print*,"     RACT, RCCN"
      do KR = 1,NKR  
         print*,KR,RACT(KR),RCCN(KR)
      end do
      print*,"RACTMAX = ", RACTMAX
      print*,"----------------------"
      CALL DROPMAX(DROPRADII,RACTMAX,NDROPMAX,NKR)
      print*,"NDROPMAX bin = ", NDROPMAX
! put nucleated droplets into the drop bin according to radius
! change in drop concentration due to activation DROPCONCN(IDROP)
      ISMALL=NCRITI

      INEXT=ISMALL
!       ISMALL=1

!       INEXT=ISMALL
        DO IDROP=1,NDROPMAX !liquid bin
           DROPCONCN(IDROP)=0.
           DO I=ISMALL,IMAX
!the activated CCN bins be accumulated and transfer to certain liquid
!bin according to radius
              IF(RACT(I).LE.DROPRADII(IDROP)) THEN
                print*,IDROP,I,RACT(I)
                DROPCONCN(IDROP)=DROPCONCN(IDROP)+CCNCONC(I)
                INEXT=I+1
              ENDIF
           ENDDO
           print*,"--------"
           ISMALL=INEXT
        ENDDO
      print*,"          BIN,   DROPRADII,   DROPCONCN"
      do KR = 1,NKR
         print*,KR,DROPRADII(KR),DROPCONCN(KR)
      end do
      print*,"----------------------"
      print*,SUP1,TT



      stop
      end program nucleation
     !---------------------------------------------------------------
      SUBROUTINE BOUNDARY(IMIN,IMAX,FCCNR,NKR)
! IMIN - left CCN spectrum boundary
      IMPLICIT NONE
      INTEGER I,IMIN,IMAX,NKR
      REAL FCCNR(NKR)

      IMIN=0

      DO I=1,NKR
         IF(FCCNR(I).NE.0.) THEN
           IMIN=I
           GOTO 40
         ENDIF
      ENDDO

40    CONTINUE

! IMAX - right CCN spectrum boundary

      IMAX=0

      DO I=NKR,1,-1
         IF(FCCNR(I).NE.0.) THEN
           IMAX=I
           GOTO 41
         ENDIF
      ENDDO

41    CONTINUE
      RETURN
      END  SUBROUTINE BOUNDARY

      SUBROUTINE CRITICAL (AKOE,BKOE,TT,RCRITI,SUP1,DEG01)
! AKOE & BKOE - constants in Koehler equation
      IMPLICIT NONE
      REAL AKOE,BKOE,TT,RCRITI,SUP1,DEG01
      REAL RO_SOLUTE
      PARAMETER (RO_SOLUTE=2.16)
      INTEGER I
      REAL :: SUP1_TEST(99),RCRITI_TEST(99),RCRITI_TEST2(99),SUP1_I

      AKOE=3.3E-05/TT
      BKOE=2.*4.3/(22.9+35.5)
! new change 21.07.02                                         (begin)
      BKOE=BKOE*(4./3.)*3.141593*RO_SOLUTE
! new change 21.07.02                                           (end)

! table of critical aerosol radii

!       GOTO 992

! SUP1_TEST(I), %
       SUP1_TEST(1)=0.42!0.01
       SUP1_TEST(2)=0.13
       SUP1_TEST(3)=0.042
       SUP1_TEST(4)=0.013
       SUP1_TEST(5)=0.0042
       print*,"SUBROUTINE CRITICAL  r* (micro-m)          ","S*-1 (%)"
       DO I=1,5
          !SUP1_TEST(I+1)=SUP1_TEST(I)+0.01
          SUP1_I=SUP1_TEST(I)*0.01
          RCRITI_TEST(I)=(AKOE/3.)*(4./BKOE/SUP1_I/SUP1_I)**DEG01
          RCRITI_TEST2(I)=sqrt(3.0*BKOE/AKOE)
          print*,RCRITI_TEST(I)*1.d4,SUP1_TEST(I),RCRITI_TEST2(I)
       ENDDO

! RCRITI, cm - critical radius of "dry" aerosol

      RCRITI=(AKOE/3.)*(4./BKOE/SUP1/SUP1)**DEG01
      RETURN
      END  SUBROUTINE CRITICAL
     !------------------------------
      SUBROUTINE CCNIMIN(IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC,COL, &
     &       FCCNR,NKR)
! FOR    IMIN=1
      IMPLICIT NONE
      INTEGER IMIN,II,IMAX,NCRITI,NKR
      REAL RCRITI,COL
      REAL RCCN(NKR),FCCNR(NKR),CCNCONC(NKR)
      REAL RCCN_MIN
      REAL DLN1,DLN2,FOLD_IP
! rccn_min - minimum aerosol(ccn) radius
      RCCN_MIN=(1.23d-7)/10000.!RCCN(1)/10000.  !seed
      !RCCN_MIN=RCCN(1)/10000. 
! calculation of ccnconc(ii)=fccnr(ii)*col - aerosol(ccn) bin
!                                            concentrations,
!                                            ii=imin,...,imax
! determination of ncriti   - number bin in which is located rcriti
! calculation of ccnconc(ncriti)=fccnr(ncriti)*dln1/(dln1+dln2),
! where,
! dln1=Ln(rcriti)-Ln(rccn_min)
! dln2=Ln(rccn(1)-Ln(rcriti)
! calculation of new value of fccnr(ncriti)

!     IF(IMIN.EQ.1) THEN
        IF(RCRITI.LE.RCCN_MIN) THEN
          NCRITI=1
          DO II=NCRITI+1,IMAX
             CCNCONC(II)=COL*FCCNR(II)
             FCCNR(II)=0.
          ENDDO
          GOTO 42
        ENDIF
        IF(RCRITI.GT.RCCN_MIN.AND.RCRITI.LT.RCCN(IMIN)) THEN
          NCRITI=1
          DO II=NCRITI+1,IMAX
             CCNCONC(II)=COL*FCCNR(II)
             FCCNR(II)=0.
          ENDDO
          DLN1=ALOG(RCRITI)-ALOG(RCCN_MIN)
          DLN2=ALOG(RCCN(1))-ALOG(RCRITI)
          CCNCONC(NCRITI)=DLN2*FCCNR(NCRITI)
          FCCNR(NCRITI)=FCCNR(NCRITI)*DLN1/max((DLN1+DLN2),1.d-10)
          GOTO 42
! in case RCRITI.GT.RCCN_MIN.AND.RCRITI.LT.RCCN(IMIN)
        ENDIF
! in case IMIN.EQ.1
42      CONTINUE

        RETURN
        END SUBROUTINE CCNIMIN

        SUBROUTINE CCNLOOP(IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC,COL, &
     &       FCCNR,NKR)
        IMPLICIT NONE
        INTEGER I,IMIN,IMAX,NKR,II,NCRITI
        REAL COL
        REAL RCRITI,RCCN(NKR),CCNCONC(NKR),FCCNR(NKR)
        REAL DLN1,DLN2,FOLD_IP
        IF(IMIN.GT.1) THEN
          IF(RCRITI.LE.RCCN(IMIN-1)) THEN
            NCRITI=IMIN
            DO II=NCRITI,IMAX
               CCNCONC(II)=COL*FCCNR(II)
               FCCNR(II)=0.
            ENDDO
            GOTO 42
          ENDIF
          IF(RCRITI.LT.RCCN(IMIN).AND.RCRITI.GT.RCCN(IMIN-1)) &
     &    THEN
! this line eliminates bug you found (when IMIN=IMAX)
            NCRITI=IMIN
           !larger than IMIN+1 are all transferred
            DO II=NCRITI+1,IMAX
               CCNCONC(II)=COL*FCCNR(II)
               FCCNR(II)=0.
            ENDDO
            DLN1=ALOG(RCRITI)-ALOG(RCCN(IMIN-1))
            DLN2=COL-DLN1
            CCNCONC(NCRITI)=DLN2*FCCNR(NCRITI)
            FCCNR(NCRITI)=FCCNR(NCRITI)*DLN1/COL
            GOTO 42
! in case RCRITI.LT.RCCN(IMIN).AND.RCRITI.GT.RCCN(IMIN-1)
          ENDIF
! in case IMIN.GT.1
        ENDIF

! END of part of interest. so in case
!RCRITI.LT.RCCN(IMIN).AND.RCRITI.GT.RCCN(IMIN-1)
!we go to 42 and avoid the next loop
         DO I=IMIN,IMAX-1
           IF(RCRITI.EQ.RCCN(I)) THEN
             NCRITI=I+1
             DO II=I+1,IMAX
                CCNCONC(II)=COL*FCCNR(II)
                FCCNR(II)=0.
             ENDDO
             GOTO 42
           ENDIF
           IF(RCRITI.GT.RCCN(I).AND.RCRITI.LT.RCCN(I+1)) THEN
             NCRITI=I+1
           !larger than I+2 are all transferred
             IF(I.NE.IMAX-1) THEN
               DO II=NCRITI+1,IMAX
                  CCNCONC(II)=COL*FCCNR(II)
                  FCCNR(II)=0.
               ENDDO
             ENDIF
             DLN1=ALOG(RCRITI)-ALOG(RCCN(I))
             DLN2=COL-DLN1
             CCNCONC(NCRITI)=DLN2*FCCNR(NCRITI)
             FCCNR(NCRITI)=FCCNR(NCRITI)*DLN1/COL
             GOTO 42
! in case RCRITI.GT.RCCN(I).AND.RCRITI.LT.RCCN(I+1)
           END IF


         ENDDO
! cycle by I, I=IMIN,...,IMAX-1

42      CONTINUE
        RETURN
        END  SUBROUTINE CCNLOOP
       !--------------------------------------------------
        SUBROUTINE ACTIVATE(IMIN,IMAX,AKOE,BKOE,RCCN,RACT,RACTMAX,NKR)
        IMPLICIT NONE

        INTEGER IMIN,IMAX,NKR
        INTEGER I,I0,I1
        REAL RCCN(NKR)
        REAL  R03,SUPCRITI,RACT(NKR),XKOE
        REAL AKOE,BKOE,AKOE23,RACTMAX
! Spectrum of activated drops                                 (begin)
        DO I=IMIN,IMAX

! critical water supersaturations appropriating CCN radii

           XKOE=(4./27.)*(AKOE**3/BKOE)
           AKOE23=AKOE*2./3.
           R03=RCCN(I)**3
           SUPCRITI=SQRT(XKOE/R03)

! RACT(I) - radii of activated drops, I=IMIN,...,IMAX

           IF(RCCN(I).LE.(0.3E-5)) &
     &     RACT(I)=AKOE23/SUPCRITI
           IF(RCCN(I).GT.(0.3E-5).and.RCCN(I).LE.(2.E-4))&
     &     RACT(I)=5.*RCCN(I) !GCCN effect
           IF(5.*RCCN(I).GT.(10.E-4))&
     &     RACT(I)=RCCN(I)!10.E-4 !UGCCN effect 5 micro-m seed
        ENDDO
! cycle by I
! calculation of I0

        I0=IMIN

        DO I=IMIN,IMAX-1
           IF(RACT(I+1).LT.RACT(I)) THEN
             I0=I+1
             GOTO 45
           ENDIF
        ENDDO

 45     CONTINUE
! new changes 9.04.02                                         (begin)
        I1=I0-1
! new changes 9.04.02                                           (end)

        IF(I0.EQ.IMIN) GOTO 47

! new changes 9.04.02                                         (begin)

        IF(I0.EQ.IMAX) THEN
          RACT(IMAX)=RACT(IMAX-1)
          GOTO 47
        ENDIF

        IF(RACT(IMAX).LE.RACT(I0-1)) THEN
          DO I=I0,IMAX
             RACT(I)=RACT(I0-1)
          ENDDO
          GOTO 47
        ENDIF

! new changes 9.04.02                                           (end)

        DO I=I0+1,IMAX
           IF(RACT(I).GE.RACT(I0-1)) THEN
             I1=I
             GOTO 46
           ENDIF
        ENDDO
 46     CONTINUE

! spectrum of activated drops                                   (end)


! line interpolation RACT(I) for I=I0,...,I1

        DO I=I0,I1
           RACT(I)=RACT(I0-1)+(I-I0+1)*(RACT(I1)-RACT(I0-1)) &
     &                       /(I1-I0+1)
        ENDDO



  47    CONTINUE



        RACTMAX=0.

        DO I=IMIN,IMAX
           RACTMAX=AMAX1(RACTMAX,RACT(I))
        ENDDO
        RETURN

        END SUBROUTINE ACTIVATE

        SUBROUTINE DROPMAX(DROPRADII,RACTMAX,NDROPMAX,NKR)
        IMPLICIT NONE
        INTEGER IDROP,NKR,NDROPMAX
        REAL RACTMAX,DROPRADII(NKR)
! calculation of NDROPMAX - maximal number of drop bin which
! is activated

        NDROPMAX=1

        DO IDROP=1,NKR
           IF(RACTMAX.LE.DROPRADII(IDROP)) THEN
             NDROPMAX=IDROP
             GOTO 44
           ENDIF
        ENDDO
 44     CONTINUE
        RETURN
        END  SUBROUTINE DROPMAX

