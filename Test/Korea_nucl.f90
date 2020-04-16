      module constants; implicit none
      integer, parameter :: nkr  = 43
     !---------------------------------------------
      integer, parameter   :: nhydro   = 7
      integer, parameter   :: WATER    = 1
      integer, parameter   :: COLUMN   = 2
      integer, parameter   :: PLATE    = 3
      integer, parameter   :: DENDRITE = 4
      integer, parameter   :: SNOW     = 5
      integer, parameter   :: GRAUPEL  = 6
      integer, parameter   :: HAIL     = 7
     !---------------------------------------------
      real, dimension(nhydro,nkr)   :: r_hydro !radius
      integer, parameter   :: scal = 1
      real,    parameter   :: dlnr = log(2.0)/3.0/float(scal)
      real, parameter   :: Rd = 287. !J/(kg K)
      real, parameter   :: pi = 3.141593
      real, parameter   :: ro_solute = 2.16
      real, parameter   :: akoe = 3.3e-05
      real, parameter   :: bkoe = 2.0*4.3/(22.9+35.5) & !NaCl
                                       &*(4.0/3.0)*pi*ro_solute
      real, parameter :: rccn_min = 0.003e-4 ! cm
      real, parameter :: smax     = 0.20     ! 120% of RH
      real, parameter :: z0_ccn   = 2000.    ! 2.0 km
      real, parameter :: dz_ccn   = 2000.    ! 2.0 km
      real, parameter :: t00 = 273.15

      real, parameter :: Lv = 2.5008e6    ! latent heat for evaporation
      real, parameter :: Cp = 1004.7
      real, parameter :: Rv    = 461.51
      real, parameter :: eps   = Rd/Rv
      real, parameter :: tt_nucl_drop_min = -45.0
      real, parameter :: b_w = Lv/Rv
      real, parameter :: a_w = b_w/t00

      end module constants

      program Korea_nucl
      use constants; implicit none
      integer n
      real, dimension(nkr,nhydro) :: temp
      real, dimension(nkr) :: mass_hydro,fccn,ff
      real tt,qq,pp
      real tt1,qq1,pp1,rhoa1,rccni,fccni,ffwi


!--------------------------------------------------------------------
      open(unit=11,file="bulkradii43.asc",status="old" &
           &,action="read")
      read(11,900) temp
      close(11)
      r_hydro = transpose(temp)
      open(unit=11,file="masses43.asc",status="old" &
           &,action="read")
      read(11,900) mass_hydro 
      close(11)

! unit conversion
! ---------------------------------------------------------
! ffx : kg kg-1, ff : # m-3
! fccnx : # kg-1, fccn : # m-3

! caution: at the initialization, the unit of fccn is # m-3, not # kg-1.
      !do m=1, nhydro
        ! ff(m,:) = max(ffx(i,k,m,:)/mass_hydro*rhoa1,0.0)
      !end do
      !if(itimestep.gt.1) then
      !fccn = max(fccnx(i,k,:)*rhoa1,0.0)
! ---------------------------------------------------------
      call fccn_init(0., 0., 0, 0, fccn)
      !do n = 1,nkr
      !   write(*,*) r_hydro(WATER,n)*100. !m --> cm
      !end do
      tt = t00+15.0
      pp = 900.0*100.0 !N/m2
      rhoa1 = pp/tt/Rd
      print*,rhoa1
      fccn = fccn * rhoa1 !# kg-1 --> # m-3
      !do n = 1,nkr
      !   write(*,*) n, fccn(n) 
      !end do
      ff(:) = 0.0   
  
      !T_old = tt      !tt(i,k)
      rccni = r_hydro(WATER,1)
      qq = 0.013
      do n=nkr, 1, -1 ! large aerosol first

         fccni = fccn(n)

         if(fccni.gt.0.0) then
            tt1   = tt !tt(i,k) !K
            qq1   = qq !qq(i,k) !m-3 ???????????
            pp1   = pp !pp(i,k) !Pa
            rhoa1 = pp1/tt1/Rd
            ffwi  = ff(max(n-35,1))!ff(WATER,max(n-35,1))

            call water_nucl(tt1, pp1 ,qq1 ,rhoa1 ,fccni, rccni, ffwi)

            tt = tt1!tt(i,k) = tt1
            qq = qq1!qq(i,k) = qq1
            fccn(n) = fccni
            ff(max(n-35,1)) = ffwi !ff(WATER,max(n-35,1)) = ffwi
         end if
   
         rccni = rccni/exp(dlnr)
      end do
      do n = 1,nkr
         write(*,*) n, fccn(n), ff(n)
      end do


900   format(6e13.5)
!--------------------------------------------------------------------
      stop
      end program Korea_nucl
     !--------------------------------------------------------------- 
      subroutine fccn_init(alt, island, ivgtype, isurban, fccni)
      use constants
      implicit none

      real, intent(in) :: alt, island
      integer, intent(in) :: isurban, ivgtype
      real, dimension(nkr), intent(out) :: fccni

      real :: a, b, a1, n0ccn, kccn, r0, s_kr, factor
      integer :: kr

      a=akoe/(t00+15.0)
      b=bkoe
      a1=2.0*(a/3.0)**1.5/sqrt(b)

      !if(alt.le.z0_ccn) then
         factor = 1.0
      !else
      !  factor = exp(-(alt-z0_ccn)/dz_ccn)
      !end if
      n0ccn = 300.0e6!#kg-1 !n0ccn1 !ocean or urban or rural
      kccn  = 0.8!kccn1
      r0=r_hydro(WATER,1)*100. ! m -> cm

      do kr=nkr, 1, -1
         s_kr=a1/(r0**1.5)
         if(r0.le.rccn_min.or.s_kr.ge.smax) then
            fccni(kr)=0.0
         else
            fccni(kr)=(1.5*n0ccn*kccn*(100.*s_kr)**kccn)*dlnr*factor
         end if
         r0=r0/exp(dlnr)
      end do

      end subroutine fccn_init
     !--------------------------------------------------------------- 
      subroutine water_nucl(tt,pp,qq,rhoa,fccni,rccni,ffwi)
      use constants
      implicit none
      real esw !function
      real, intent(in) :: pp, rhoa, rccni
      real, intent(inout) :: tt, qq, fccni, ffwi
      real, dimension(nkr) :: mass_hydro 

      real :: supw, akoe1, bkoe1, rcrit, rccni_p, dmr, qq1, tt1, factor
      open(unit=11,file="masses43.asc",status="old" &
           &,action="read")
      read(11,900) mass_hydro 
      close(11)

      supw = (pp*qq/(eps+qq))/esw(tt)-1.0
      !print*,"supw",supw
      if(supw.le.0.0.or.tt.lt.t00+tt_nucl_drop_min) return

      akoe1=akoe/tt
      bkoe1=bkoe
      rcrit=(akoe1/3.0)*(4.0/bkoe1/(supw**2))**(1.0/3.0)*1.0e-2 ! cm -> m

      rccni_p = rccni/exp(dlnr)

      if(rcrit.lt.rccni_p) then
         factor = 1.0
      else if(rcrit.lt.rccni) then
         factor = (log(rccni)-log(rcrit))/(log(rccni)-log(rccni_p))
      else
         factor = 0.0
         return
      end if

      do
        dmr  = fccni*factor*mass_hydro(1)/rhoa
        qq1  = qq - dmr
        tt1  = tt + dmr*Lv/Cp
        supw = (pp*qq1/(eps+qq1))/esw(tt1)-1.0
        if(supw.ge.0.0) exit
           factor = factor * 0.95
      end do

      qq = qq1
      tt = tt1

      ffwi  = ffwi + fccni*factor
      fccni = fccni*(1.0-factor)
      
      900   format(6e13.5)
      end subroutine water_nucl
        real function esw(tt)
        use constants
        implicit none

        real, intent(in) :: tt

        esw = 611.*exp(a_w - b_w/tt)

        end function esw
 
