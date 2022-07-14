
      subroutine mie_scat(wl,qsca,qext,RI,cssca,csabs)
        implicit none
        include 'param.inc'
  
        integer :: igfunc,nBins,N,iter,radIter,j
        double precision alpha_s,a_ts,a_cs, beta_s
        double precision :: gfunc_norm,gs_norm,rmean,ln_sigma
        common/grain_dist/gfunc_norm,rmean,ln_sigma,igfunc

        double precision x,y,qEtot,qStot,csTOT,caTOT,wRadT,wBinT
        complex m,b,i
        double precision, allocatable, dimension(:) :: Y_n, J_n
        complex, allocatable, dimension(:) :: a_n, b_n, AA_n

        real*8, dimension(ndust,2) :: RI
        real*8, dimension(4) :: abund
        double precision wl,nn,kk,qsca,qabs,qext,cssca,csabs
        double precision wln,wlk,xk,xn,pi,alpha,rad,radStep
        double precision wRad,wBin,binSize,radmin,radmax,qAdd

        nBins = 20
        i=cmplx(0,1)
        pi = 4.0 * atan(1.0)

        ! todo, add more options for abundances
        abund = (/5.0, 95.00,0.5,0.5 /)

        if(igfunc.eq.1) then ! ISM power law
           radmin = 6.25d-3 ! micron
           radmax = 2.5d-1 ! micron
           alpha = -3.5

         elseif(igfunc.eq.2) then ! log normal dist
           radmin = dexp(dlog(rmean) - 2.d0*ln_sigma)
           radmax = dexp(dlog(rmean) + 2.d0*ln_sigma)
           alpha = 1
 
         elseif(igfunc.eq.3) then ! Weingartner & Draine distribution
           radmin = 3.5d-4 ! micron
           alpha = -2.1
           radmax = 1.7d-3 + 7.0**0.333 * 1.d-1

         elseif(igfunc.eq.4) then ! Kim et al. distribution
           radmin = 2.25d-3
           radmax = 1.0
           alpha = 1.0

         elseif(igfunc.eq.5) then ! single grain
        
           radmin = 1 ! micron
           radmax = radmin
           nBins   = 1 ! num bins
           alpha  = 0.0 ! size dist. exp, n(a).propto.a^alpha

        endif
        
        write(*,*)'got igfunc',igfunc

        qEtot = 0.0
        qStot = 0.0
        caTOT = 0.0
        csTOT = 0.0
        wRadT = 0.0
        wBinT = 0.0
        radStep = (dlog(radmax)-dlog(radmin))/(nBins-1) 
        do j=1,2!ndust

           nn = RI(j,1)
           kk = RI(j,2)
           
           do radIter=1,nBins
              if(nBins.eq.1) then
                   wBin=1.
                   rad=radmin
              else
                   binSize =dexp(dlog(radmin)+radStep*(radIter))-
     &                 dexp(dlog(radmin)+radStep*(radIter-1))
                   rad = dexp(dlog(radmin)+radStep*(radIter-1))
                   wBin = abund(j)*rad**alpha*binSize 
                   !abund(j) used to be 100.
              endif
              wBinT = wBinT + wBin
              wRad = pi * (rad*1.d-6)**2*wBin
              wRadT = wRadT + wRad
              x = 2.* rad * pi/wl
              m = cmplx(nn,kk)
              b = cmplx(nn,-kk)

              y = sqrt(m*b)*x
              if(y<1) then 
                   N = 7.5*y+9.
              else if((y>1.).and.(y<100.)) then 
                   N = 1.25*y+15.5
              else if((y>100.).and.(y<50000.)) then 
                   N = 1.0625*y+28.5
              else 
                   N = 1.005*y+50.0
              endif

              allocate(Y_n(N+2),J_n(N+1),AA_n(N))
              allocate(a_n(N),b_n(N))

              Y_n(1)=-sqrt(2./(pi*x))*cos(x)
              Y_n(2)= sqrt(2./(pi*x))*(-cos(x)/x - sin(x))
              AA_n(N)=(N+1)/(m*x)
              J_n(1)=sqrt(2./(pi*x))*sin(x)
              do iter=2, N
                Y_n(iter+1)=(2.*(iter-1.)+1.)/x*Y_n(iter)-Y_n(iter-1)
                J_n(iter)= (2./(pi*x)+Y_n(iter)*J_n(iter-1))/Y_n(iter-1)
              enddo

              do iter=N-1,1,-1
                   AA_n(iter)=(iter+1)/(m*x)
     &              -(AA_n(iter+1)+(iter+1)/(m*x))**(-1)
              enddo

              qext=0.
              qsca=0.

              do iter=1,N
                 a_n(iter)=((AA_n(iter)/m+iter/x)
     &             *J_n(iter+1)-J_n(iter))/
     &             ((AA_n(iter)/m+iter/x)*J_n(iter+1)-J_n(iter)+ 
     &             i*((AA_n(iter)/m+iter/x)*Y_n(iter+1)-Y_n(iter)))
                 b_n(iter)=((AA_n(iter)*m+iter/x)
     &             *J_n(iter+1)-J_n(iter))/
     &             ((AA_n(iter)*m+iter/x)*J_n(iter+1)-J_n(iter)+
     &             i*((AA_n(iter)*m+iter/x)*Y_n(iter+1)-Y_n(iter)))
                 qAdd = (2.*iter+1.)*REAL(a_n(iter)+b_n(iter))
                 qext=qext+qAdd

                 qsca=qsca+(2.*iter+1.)*
     &           (a_n(iter)*conjg(a_n(iter))+b_n(iter)*conjg(b_n(iter)))
                 if(qAdd/qext.lt.1E-20) then
                     exit
                 endif
               enddo
               qEtot = qEtot + qext * 2./x**2. * wRad
               qStot = qStot + qsca*2./x**2. * wRad
               qabs = qEtot-qStot
               csTOT = pi*(rad)**2.*qsca*2./x**2. * wRad
               caTOT = caTOT + qabs*pi*(rad)**2.
               deallocate(Y_n,J_n,AA_n,a_n,b_n) 
           enddo
        enddo

        qext = qEtot/wRadT
        qsca = qStot/wRadT
        cssca = qStot/wBinT
        csabs = (qEtot-qStot)/wBinT 
      end subroutine mie_scat

