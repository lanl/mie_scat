      subroutine miegrn(wl, ag,sg)

c      use vcapcom
      implicit none
      include 'param.inc'
      real*8, intent(in) :: wl
      real*8, intent(out) :: ag(:),sg(:)
*-- input:
* wl: wavelength in micron 
*-- output:
* ag,sg: cross-sections, opacity in 1/cm/particule
* ag: absorption cross-section
* sg: scattering cross-section
      integer :: J

      real*8 qn,qk
      character*40 head
      integer nn,nk
      real*8 wln,xn,wlk,xk
      real*8, dimension(ndust,2) :: RI
      double precision qsca,qext,cssca,csabs

c     initialize vectors:
      ag(:) = 0.d0
      sg(:) = 0.d0
      do j=1,ndust
         call optcon(wl,j)
         RI(j,1)=qn
         RI(j,2)=qk
      enddo
      call mie_scat(wl,qsca,qext,RI,cssca,csabs)
      write(6,*)wl,qsca,qext
      sg(1)=cssca
      ag(1)=csabs

      return
      end
