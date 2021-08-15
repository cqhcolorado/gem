subroutine pushe(icycle,irk,ncycle)
  use gem_com
  use gem_equil
  implicit none
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dgdtp,dpdzp,dadzp,aparp,vncp
  real :: dum,vxdum,dum1,dum2,eps,x,h_x,h_coll,hee,nue
  INTEGER :: m,i,j,k,l,n,ipover,ieover
  INTEGER :: np_old,np_new
  real :: vfac,kap,kapnp,kaptp,sz,vpar,ppar,vpdum,pidum,xnp
  real :: b,th,r,enerb,cost,sint,qr,laps,ter
  real :: xt,xs,yt,zt,xdot,ydot,zdot,pzdot,edot,pzd0,pzd1,vp0
  real :: wx0,wx1,wy0,wy1,wz0,wz1
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
  real :: grp,gxdgyp,psp,pzp,psip2p,bdcrvbp,curvbzp,dipdrp,bstar
  real :: myavptch,myaven
  integer :: myopz,myoen
  real :: x000,x001,x010,x011,x100,x101,x110,x111
  
  integer :: icycle,irk,ncycle,ierr
  real :: dte_cycle

  integer, save :: mme_tmp,ifirst=0
  real,dimension(:),allocatable :: x0e,y0e,z0e0,u0e0,w0e0,&
          mue0,index0,ipass0,mue00,xie0,z0e00, &
          eke0,pze0,u0e00,pif0e0,w5e0,w0e00
 
  save x0e,y0e,z0e0,u0e0,w0e0,mue0,index0,ipass0,mue00,xie0,z0e00,eke0,pze0,u0e00,pif0e0,w5e0,w0e00

  if(ifirst .ne. -99)then
     allocate(x0e(1:mmxe),y0e(1:mmxe),z0e0(1:mmxe),u0e0(1:mmxe),w0e0(1:mmxe),&
          mue0(1:mmxe),index0(1:mmxe),ipass0(1:mmxe),mue00(1:mmxe),xie0(1:mmxe),& 
          z0e00(1:mmxe),eke0(1:mmxe),pze0(1:mmxe),u0e00(1:mmxe),pif0e0(1:mmxe),w5e0(1:mmxe),w0e00(1:mmxe))
     ifirst=-99
     if(myid==0)then
       write(gemout,*)'save middle value'
       call flush(gemout)
     endif 
  endif
  if(irk==1)then
    dte_cycle=0.5*dte/real(ncycle)
    if(icycle==1)then
      !save data
!      !$acc update host(x2e,y2e,z2e,u2e,w2e,mue2,index,ipass,xie,z0e,eke,pze,u0e,pif0e,w5e,w0e)
      !$omp target update from(x2e,y2e,z2e,u2e,w2e,mue2,index,ipass,xie,z0e,eke,pze,u0e,pif0e,w5e,w0e)
      mme_tmp=mme
      x0e=x2e
      y0e=y2e
      z0e0=z2e
      u0e0=u2e
      w0e0=w2e
      mue0=mue2
      index0=index
      ipass0=ipass
      !mue00=mue3
      xie0=xie
      z0e00=z0e
      eke0=eke
      pze0=pze
      u0e00=u0e
      pif0e0=pif0e
      w5e0=w5e
      w0e00=w0e
      if(myid==0)then
        write(gemout,*)'save for RK-1'
        call flush(gemout)
      endif
    endif   
  else
    if(icycle==1)then
      !restore data
      mme=mme_tmp
      x2e=x0e
      y2e=y0e
      z2e=z0e0
      u2e=u0e0
      w2e=w0e
      mue2=mue0
      index=index0
      ipass=ipass0
      !mue3=mue00
      xie=xie0
      z0e=z0e00
      eke=eke0
      pze=pze0
      u0e=u0e00
      pif0e=pif0e0
      w5e=w5e0
      w0e=w0e00
!      !$acc update device(x2e,y2e,z2e,u2e,w2e,mue2,index,ipass,xie,z0e,eke,pze,u0e,pif0e,w5e,w0e)
     !$omp target update to(x2e,y2e,z2e,u2e,w2e,mue2,index,ipass,xie,z0e,eke,pze,u0e,pif0e,w5e,w0e)
      if(myid==0)then
        write(gemout,*)'restore for RK-0'
        call flush(gemout)
      endif
    endif
    dte_cycle=dte/real(ncycle) 
  endif
 
  myopz = 0
  myoen = 0
  myavptch = 0.
  myaven = 0.
  pidum = 1./(pi*2)**1.5*(vwidthe)**3
!  !$acc parallel loop gang vector
  !$omp target data map(tofrom:w2e(1:mme),u3e(1:mme),mue3(1:mme),z3e(1:mme),z2e(1:mme),x2e(1:mme),w3e(1:mme),y3e(1:mme),x3e(1:mme),y2e(1:mme),u2e(1:mme)) map(to:zeff(:),mue2(1:mme),jfn(:),ex(:,:,0:1),psip(:),capne(:),phincp(:),grdgt(:,:),grcgt(:,:),capte(:),f(:),dipdr(:),bdcrvb(:,:),dbdr(:,:),t0e(:),curvbz(:,:),gr(:,:),qhat(:,:),psip2(:),ez(:,:,0:1),psi(:),sf(:),bfld(:,:),nue0(:),dbdth(:,:),thfnz(:),dydr(:,:),ey(:,:,0:1),radius(:,:))
  !$omp target teams loop
  do m=1,mme
     r=x2e(m)-0.5*lx+lr0

     k = int(z2e(m)/delz)
     wz0 = ((k+1)*delz-z2e(m))/delz
     wz1 = 1-wz0
     th = wz0*thfnz(k)+wz1*thfnz(k+1)

     i = int((r-rin)/dr)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     k = int((th+pi)/dth)
     wz0 = (-pi+(k+1)*dth-th)/dth
     wz1 = 1.-wz0
     dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
          +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1)
     dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
          +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1)
     grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
          +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1)
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
     radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
          +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1)
     dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
          +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1)
     qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
          +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1)
     grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
          +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1)
     gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
          +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1)

     curvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
          +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1)
     bdcrvbp = wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
          +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1)
     grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
          +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1)

     fp = wx0*f(i)+wx1*f(i+1)
     jfnp = wz0*jfn(k)+wz1*jfn(k+1)
     psipp = wx0*psip(i)+wx1*psip(i+1)
     psp = wx0*psi(i)+wx1*psi(i+1)
     ter = wx0*t0e(i)+wx1*t0e(i+1)
     kaptp = wx0*capte(i)+wx1*capte(i+1)
     kapnp = wx0*capne(i)+wx1*capne(i+1)
     xnp = wx0*xn0e(i)+wx1*xn0e(i+1)
     vncp = wx0*phincp(i)+wx1*phincp(i+1)
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)
     b=1.-tor+tor*bfldp
     pzp = emass*u2e(m)/b+psp/br0
     xt = x2e(m)
     yt = y2e(m)

#ifdef JYCHENG
     include 'ppushlie.h'
#else
     i=int(xt/dx)
     j=int(yt/dy)
     k=0 !int(z3e(m)/dz)-gclr*kcnt

     wx0=real(i+1)-xt/dx
     wx1=1.-wx0
     wy0=real(j+1)-yt/dy
     wy1=1.-wy0
     wz0=real(gclr*kcnt+k+1)-z2e(m)/dz
     wz1=1.-wz0
     x000=wx0*wy0*wz0
     x001=wx0*wy0*wz1
     x010=wx0*wy1*wz0
     x011=wx0*wy1*wz1
     x100=wx1*wy0*wz0
     x101=wx1*wy0*wz1
     x110=wx1*wy1*wz0
     x111=wx1*wy1*wz1
     exp1 = x000*ex(i,j,k) + x100*ex(i+1,j,k)  &
           + x010*ex(i,j+1,k) + x110*ex(i+1,j+1,k) + &
           x001*ex(i,j,k+1) + x101*ex(i+1,j,k+1) +   &
           x011*ex(i,j+1,k+1) + x111*ex(i+1,j+1,k+1)

     eyp = x000*ey(i,j,k) + x100*ey(i+1,j,k)  &
           + x010*ey(i,j+1,k) + x110*ey(i+1,j+1,k) + &
           x001*ey(i,j,k+1) + x101*ey(i+1,j,k+1) +   &
           x011*ey(i,j+1,k+1) + x111*ey(i+1,j+1,k+1)

     ezp = x000*ez(i,j,k) + x100*ez(i+1,j,k)  &
           + x010*ez(i,j+1,k) + x110*ez(i+1,j+1,k) + &
           x001*ez(i,j,k+1) + x101*ez(i+1,j,k+1) +   &
           x011*ez(i,j+1,k+1) + x111*ez(i+1,j+1,k+1)

     !delbxp = x000*delbx(i,j,k) + x100*delbx(i+1,j,k)  &
     !      + x010*delbx(i,j+1,k) + x110*delbx(i+1,j+1,k) + &
     !      x001*delbx(i,j,k+1) + x101*delbx(i+1,j,k+1) +   &
     !      x011*delbx(i,j+1,k+1) + x111*delbx(i+1,j+1,k+1)

     !delbyp = x000*delby(i,j,k) + x100*delby(i+1,j,k)  &
     !      + x010*delby(i,j+1,k) + x110*delby(i+1,j+1,k) + &
     !      x001*delby(i,j,k+1) + x101*delby(i+1,j,k+1) +   &
     !      x011*delby(i,j+1,k+1) + x111*delby(i+1,j+1,k+1)

     !dgdtp = x000*dphidt(i,j,k) + x100*dphidt(i+1,j,k)  &
     !      + x010*dphidt(i,j+1,k) + x110*dphidt(i+1,j+1,k) + &
     !      x001*dphidt(i,j,k+1) + x101*dphidt(i+1,j,k+1) +   &
     !      x011*dphidt(i,j+1,k+1) + x111*dphidt(i+1,j+1,k+1)

     !dpdzp = x000*dpdz(i,j,k) + x100*dpdz(i+1,j,k)  &
     !      + x010*dpdz(i,j+1,k) + x110*dpdz(i+1,j+1,k) + &
     !      x001*dpdz(i,j,k+1) + x101*dpdz(i+1,j,k+1) +   &
     !      x011*dpdz(i,j+1,k+1) + x111*dpdz(i+1,j+1,k+1)

     !dadzp = x000*dadz(i,j,k) + x100*dadz(i+1,j,k)  &
     !      + x010*dadz(i,j+1,k) + x110*dadz(i+1,j+1,k) + &
     !      x001*dadz(i,j,k+1) + x101*dadz(i+1,j,k+1) +   &
     !      x011*dadz(i,j+1,k+1) + x111*dadz(i+1,j+1,k+1)

     !aparp = x000*apar(i,j,k) + x100*apar(i+1,j,k)  &
     !      + x010*apar(i,j+1,k) + x110*apar(i+1,j+1,k) + &
     !      x001*apar(i,j,k+1) + x101*apar(i+1,j,k+1) +   &
     !      x011*apar(i,j+1,k+1) + x111*apar(i+1,j+1,k+1)

     !phip = x000*phi(i,j,k) + x100*phi(i+1,j,k)  &
     !      + x010*phi(i,j+1,k) + x110*phi(i+1,j+1,k) + &
     !      x001*phi(i,j,k+1) + x101*phi(i+1,j,k+1) +   &
     !      x011*phi(i,j+1,k+1) + x111*phi(i+1,j+1,k+1)
#endif

     vfac = 0.5*(emass*u2e(m)**2 + 2.*mue2(m)*b)
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw
     kap = kapnp - (1.5-vfac/ter)*kaptp

     ppar = u2e(m)
     vpar = u2e(m)-qel/emass*aparp*ipara
     bstar = b*(1+emass*vpar/(qel*b)*bdcrvbp)
     enerb=(mue2(m)+emass*vpar*vpar/b)/qel*b/bstar*tor
     dum2 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
     vxdum = (eyp/b+vpar/b*delbxp*iflut)*dum2
     xdot = vxdum*nonline  &
          -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
     ydot = (-exp1/b+vpar/b*delbyp*iflut)*dum2*nonline &
          +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
          (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
          +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
          -emass*vpar**2/(qel*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*grcgtp*lr0/q0*qhatp

     zdot =  vpar*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
          +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*q0*br0*grcgtp/jfnp

     pzd0 =tor*(-mue2(m)/emass/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
          +mue2(m)*vpar/(qel*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
     pzd1 = qel/emass*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
          +qel/emass*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)
     pzdot = pzd0+pzd1

     vpdum = vpar
     edot = qel*(xdot*exp1+(ydot-vp0)*eyp+zdot*ezp) &
          +qel*pzdot*aparp &
          +qel*vpdum*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)   &
          -qel*vpar*delbxp*vp0

     x3e(m) = x2e(m) + dte_cycle*xdot
     y3e(m) = y2e(m) + dte_cycle*ydot
     z3e(m) = z2e(m) + dte_cycle*zdot
     u3e(m) = u2e(m) + dte_cycle*pzdot
     mue3(m) = mue2(m)

     eps = (b*mue2(m)+0.5*emass*u2e(m)*u2e(m))/ter
     x = sqrt(eps)
     h_x    = 4.0*x*x/(3.0*sqrt(pi))
     h_coll = h_x/sqrt(1.0+h_x**2)
     !         h_coll =
     !         exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*DERF(x)
     ! collision frequency for experimental profiles
     hee = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*erf(x)
     nue=rneu*(wx0*nue0(i)+wx1*nue0(i+1))*(wx0*zeff(i)+wx1*zeff(i+1)+hee)
     !         dum1 = 0.5*rneu/eps**1.5*(1+h_coll)
     dum1 = nue/(2*(eps+0.1))**1.5  !*(1+h_coll)
     !         if(x<0.3)dum1=0.0

     dum = 1-w2e(m)*nonline*0.
     if(eldu.eq.1)dum = (tge/ter)**1.5*exp(vfac*(1/tge-1./ter))
     vxdum = (eyp/b+vpdum/b*delbxp)*dum2
#ifndef DIRECT_F
     w3e(m)=w2e(m) + dte_cycle*(  &
          (vxdum*kap + edot/ter-dum1*ppar*aparp/ter)*xnp     &
          +isg*(-dgdtp-zdot*dpdzp+xdot*exp1+ydot*eyp)/ter*xnp-tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp)*dum
#endif
     !         if(x3e(m)>lx .or. x3e(m)<0.)w3e(m) = 0. 

     !         go to 333

!     ieover = 0
!     ipover = 0
!     if(abs(pzp-pze(m))>pzcrite)then
!       myopz = myopz+1
!        ipover = 1
!     end if
!     if(abs(vfac-eke(m))>encrit.and.abs(x2e(m)-lx/2)<(lx/2-1))then
!        myoen = myoen+1
!        ieover = 1
!        myaven = myaven+eke(m)
!        myavptch = myavptch+abs(vpar)/sqrt(2/emass*vfac)
!     end if
!     if(itube==1)goto 333
!     if(ieover==1.or.ipover==1)then
!        x3e(m) = xie(m)
!        z3e(m) = z0e(m)
!        r = x3e(m)-lx/2+lr0
!        k = int(z3e(m)/delz)
!        wz0 = ((k+1)*delz-z3e(m))/delz
!        wz1 = 1-wz0
!        th = wz0*thfnz(k)+wz1*thfnz(k+1)

!        i = int((r-rin)/dr)
!        wx0 = (rin+(i+1)*dr-r)/dr
!        wx1 = 1.-wx0
!        k = int((th+pi)/dth)
!        wz0 = (-pi+(k+1)*dth-th)/dth
!        wz1 = 1.-wz0
!        b = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
!             +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
!        u3e(m) = u0e(m)
!        u2e(m) = u3e(m)
!        w3e(m) = 0.
!        w2e(m) = 0.
!        x2e(m) = x3e(m)
!        z2e(m) = z3e(m)
!        mue3(m) = mue(m)
!        mue2(m) = mue(m)
!     end if
333  continue
     laps=anint((z3e(m)/lz)-.5)*(1-peritr)
     r=x3e(m)-0.5*lx+lr0
     i = int((r-rin)/dr)
     i = min(i,nr-1)
     i = max(i,0)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     y3e(m)=mod(y3e(m)-laps*2*pi*qr*lr0/q0*sign(1.0,q0)+8000.*ly,ly)
     if(x3e(m)>lx.and.iperidf==0)then
        x3e(m) = lx-1.e-8
        z3e(m)=lz-z3e(m)
        x2e(m) = x3e(m)
        z2e(m) = z3e(m)
        w2e(m) = 0.
        w3e(m) = 0.
     end if
     if(x3e(m)<0..and.iperidf==0)then
        x3e(m) = 1.e-8
        z3e(m)=lz-z3e(m)
        x2e(m) = x3e(m)
        z2e(m) = z3e(m)
        w2e(m) = 0.
        w3e(m) = 0.
     end if
     z3e(m)=mod(z3e(m)+8.*lz,lz)
     x3e(m)=mod(x3e(m)+800.*lx,lx)
     x3e(m) = min(x3e(m),lx-1.0e-8)
     y3e(m) = min(y3e(m),ly-1.0e-8)
     z3e(m) = min(z3e(m),lz-1.0e-8)

     u2e(m)=u3e(m)
     x2e(m)=x3e(m)
     y2e(m)=y3e(m)
     z2e(m)=z3e(m)
     w2e(m)=w3e(m)

  end do
!  !$acc end parallel
  !$omp end target teams loop
  !$omp end target data
  call MPI_ALLREDUCE(myopz,nopz,1,MPI_integer, &
       MPI_SUM, MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(myoen,noen,1,MPI_integer, &
       MPI_SUM, MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(myaven,aven,1,MPI_REAL8, &
       MPI_SUM, MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(myavptch,avptch,1,MPI_REAL8, &
       MPI_SUM, MPI_COMM_WORLD,ierr)
  aven = aven/(noen+0.1)
  avptch = avptch/(noen+0.1)

  np_old=mme
#ifndef OLD_PMOVE
  call test_init_pmove(z3e,np_old,lz,ierr)
  call test_pmove(x2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(x3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(y2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(y3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(z2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(z3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(u2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(u3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(w2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(w3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(mue2,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(mue3,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(index,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(ipass,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(mue,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(xie,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(z0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(eke,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(pze,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(u0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  

  call test_pmove(w0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
#ifdef DIRECT_F
  call test_pmove(w5e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(pif0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
#endif


#else
  call init_pmove(z3e,np_old,lz,ierr)
  call pmove(x2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(x3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(y2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(y3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(u2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(u3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(mue2,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(mue3,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(index,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(ipass,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(mue,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(xie,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(eke,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(pze,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(u0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit


  call pmove(w0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
#ifdef DIRECT_F
  call pmove(w5e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(pif0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
#endif
#endif
  call end_pmove(ierr)
  mme=np_new

  !      return
end subroutine pushe

#ifdef RK4
subroutine pushe_rk4_1step(icycle,irk,ncycle)
  use gem_com
  use gem_equil
  implicit none
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dgdtp,dpdzp,dadzp,aparp,vncp
  real :: dum,vxdum,dum1,dum2,eps,x,h_x,h_coll,hee,nue
  INTEGER :: m,i,j,k,l,n,ipover,ieover
  INTEGER :: np_old,np_new
  real :: vfac,kap,kapnp,kaptp,sz,vpar,ppar,vpdum,pidum,xnp
  real :: b,th,r,enerb,cost,sint,qr,laps,ter
  real :: xt,xs,yt,zt,xdot,ydot,zdot,pzdot,edot,pzd0,pzd1,vp0
  real :: wx0,wx1,wy0,wy1,wz0,wz1
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
  real :: grp,gxdgyp,psp,pzp,psip2p,bdcrvbp,curvbzp,dipdrp,bstar
  real :: myavptch,myaven
  integer :: myopz,myoen
  real :: x000,x001,x010,x011,x100,x101,x110,x111
  real :: x2_k0,y2_k0,z2_k0,u2_k0,mue2_k0
  real :: xdot_k1,ydot_k1,zdot_k1,pzdot_k1,dt_k1,x3_k1,y3_k1,z3_k1,u3_k1
  real :: xdot_k2,ydot_k2,zdot_k2,pzdot_k2,dt_k2,x3_k2,y3_k2,z3_k2,u3_k2
  real :: xdot_k3,ydot_k3,zdot_k3,pzdot_k3,dt_k3,x3_k3,y3_k3,z3_k3,u3_k3
  real :: xdot_k4,ydot_k4,zdot_k4,pzdot_k4,dt_k4

  integer :: icycle,irk,ncycle
  real :: dte_cycle

  integer, save :: mme_tmp,ifirst=0
  real,dimension(:),allocatable :: x0e,y0e,z0e0,u0e0,w0e0,&
          mue0,index0,ipass0,xie0,z0e00, &
          eke0,pze0,u0e00,pif0e0,w5e0,w0e00
 
  save x0e,y0e,z0e0,u0e0,w0e0,index0,ipass0,mue0,xie0,z0e00,eke0,pze0,u0e00,pif0e0,w5e0,w0e00

  
  if(ifirst .ne. -99)then
     allocate(x0e(1:mmxe),y0e(1:mmxe),z0e0(1:mmxe),u0e0(1:mmxe),w0e0(1:mmxe),&
          mue0(1:mmxe),index0(1:mmxe),ipass0(1:mmxe),xie0(1:mmxe),& 
          z0e00(1:mmxe),eke0(1:mmxe),pze0(1:mmxe),u0e00(1:mmxe),pif0e0(1:mmxe),w5e0(1:mmxe),w0e00(1:mmxe))
     ifirst=-99
     if(myid==0)then
       write(gemout,*)'save middle value'
       call flush(gemout)
     endif 
  endif

  if(irk==1)then
    dte_cycle=0.5*dte/real(ncycle)
    if(icycle==1)then
      !save data
      !$acc update host(x2e,y2e,z2e,u2e,w2e,mue2,index,ipass,xie,z0e,eke,pze,u0e,pif0e,w5e,w0e)
      mme_tmp=mme
      x0e=x2e
      y0e=y2e
      z0e0=z2e
      u0e0=u2e
      w0e0=w2e
      mue0=mue2
      index0=index
      ipass0=ipass
      !mue00=mue0
      xie0=xie
      z0e00=z0e
      eke0=eke
      pze0=pze
      u0e00=u0e
      pif0e0=pif0e
      w5e0=w5e
      w0e00=w0e
      if(myid==0)then
       write(gemout,*)'save data'
       call flush(gemout)
     endif
    endif   
  else
    if(icycle==1)then
      !restore data
      mme=mme_tmp
      x2e=x0e
      y2e=y0e
      z2e=z0e0
      u2e=u0e0
      w2e=w0e0
      mue2=mue0
      index=index0
      ipass=ipass0
      !mue0=mue00
      xie=xie0
      z0e=z0e00
      eke=eke0
      pze=pze0
      u0e=u0e00
      pif0e=pif0e0
      w5e=w5e0
      w0e=w0e00
      !$acc update device(x2e,y2e,z2e,u2e,w2e,mue2,index,ipass,xie,z0e,eke,pze,u0e,pif0e,w5e,w0e)
      if(myid==0)then
       write(gemout,*)'restore data'
       call flush(gemout)
     endif
    endif
    dte_cycle=dte/real(ncycle) 
  endif
 
  myopz = 0
  myoen = 0
  myavptch = 0.
  myaven = 0.
  pidum = 1./(pi*2)**1.5*(vwidthe)**3
  !$acc parallel loop gang vector
  do m=1,mme

!calculate k1

     x2_k0=x2e(m)
     y2_k0=y2e(m)
     z2_k0=z2e(m)
     u2_k0=u2e(m)
     mue2_k0=mue2(m)
   
     r=x2_k0-0.5*lx+lr0
     k = int(z2_k0/delz)
     wz0 = ((k+1)*delz-z2_k0)/delz
     wz1 = 1-wz0
     th = wz0*thfnz(k)+wz1*thfnz(k+1)

     i = int((r-rin)/dr)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     k = int((th+pi)/dth)
     wz0 = (-pi+(k+1)*dth-th)/dth
     wz1 = 1.-wz0
     dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
          +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1)
     dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
          +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1)
     grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
          +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1)
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
     radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
          +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1)
     dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
          +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1)
     qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
          +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1)
     grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
          +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1)
     gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
          +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1)

     curvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
          +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1)
     bdcrvbp = wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
          +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1)
     grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
          +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1)

     fp = wx0*f(i)+wx1*f(i+1)
     jfnp = wz0*jfn(k)+wz1*jfn(k+1)
     psipp = wx0*psip(i)+wx1*psip(i+1)
     psp = wx0*psi(i)+wx1*psi(i+1)
     ter = wx0*t0e(i)+wx1*t0e(i+1)
     kaptp = wx0*capte(i)+wx1*capte(i+1)
     kapnp = wx0*capne(i)+wx1*capne(i+1)
     xnp = wx0*xn0e(i)+wx1*xn0e(i+1)
     vncp = wx0*phincp(i)+wx1*phincp(i+1)
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)
     b=1.-tor+tor*bfldp
     pzp = emass*u2_k0/b+psp/br0
     xt = x2_k0
     yt = y2_k0
     zt = z2_k0
#ifdef JYCHENG
     include 'ppushlie.h'
#else
     i=int(xt/dx)
     j=int(yt/dy)
     k=int(zt/dz) !int(z3e(m)/dz)-gclr*kcnt

     wx0=real(i+1)-xt/dx
     wx1=1.-wx0
     wy0=real(j+1)-yt/dy
     wy1=1.-wy0
     wz0=real(k+1)-zt/dz
     wz1=1.-wz0
     x000=wx0*wy0*wz0
     x001=wx0*wy0*wz1
     x010=wx0*wy1*wz0
     x011=wx0*wy1*wz1
     x100=wx1*wy0*wz0
     x101=wx1*wy0*wz1
     x110=wx1*wy1*wz0
     x111=wx1*wy1*wz1
     exp1 = x000*ex_3d(i,j,k) + x100*ex_3d(i+1,j,k)  &
           + x010*ex_3d(i,j+1,k) + x110*ex_3d(i+1,j+1,k) + &
           x001*ex_3d(i,j,k+1) + x101*ex_3d(i+1,j,k+1) +   &
           x011*ex_3d(i,j+1,k+1) + x111*ex_3d(i+1,j+1,k+1)

     eyp = x000*ey_3d(i,j,k) + x100*ey_3d(i+1,j,k)  &
           + x010*ey_3d(i,j+1,k) + x110*ey_3d(i+1,j+1,k) + &
           x001*ey_3d(i,j,k+1) + x101*ey_3d(i+1,j,k+1) +   &
           x011*ey_3d(i,j+1,k+1) + x111*ey_3d(i+1,j+1,k+1)

     ezp = x000*ez_3d(i,j,k) + x100*ez_3d(i+1,j,k)  &
           + x010*ez_3d(i,j+1,k) + x110*ez_3d(i+1,j+1,k) + &
           x001*ez_3d(i,j,k+1) + x101*ez_3d(i+1,j,k+1) +   &
           x011*ez_3d(i,j+1,k+1) + x111*ez_3d(i+1,j+1,k+1)

#endif

     vfac = 0.5*(emass*u2_k0**2 + 2.*mue2_k0*b)
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw
     kap = kapnp - (1.5-vfac/ter)*kaptp

     ppar = u2_k0
     vpar = u2_k0-qel/emass*aparp*ipara
     bstar = b*(1+emass*vpar/(qel*b)*bdcrvbp)
     enerb=(mue2_k0+emass*vpar*vpar/b)/qel*b/bstar*tor
     dum2 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
     vxdum = (eyp/b+vpar/b*delbxp*iflut)*dum2
     xdot = vxdum*nonline  &
          -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
     ydot = (-exp1/b+vpar/b*delbyp*iflut)*dum2*nonline &
          +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
          (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
          +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
          -emass*vpar**2/(qel*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*grcgtp*lr0/q0*qhatp

     zdot =  vpar*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
          +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*q0*br0*grcgtp/jfnp

     pzd0 =tor*(-mue2_k0/emass/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
          +mue2_k0*vpar/(qel*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
     pzd1 = qel/emass*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
          +qel/emass*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)
     pzdot = pzd0+pzd1

     vpdum = vpar
     edot = qel*(xdot*exp1+(ydot-vp0)*eyp+zdot*ezp) &
          +qel*pzdot*aparp &
          +qel*vpdum*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)   &
          -qel*vpar*delbxp*vp0

     xdot_k1=xdot
     ydot_k1=ydot
     zdot_k1=zdot
     pzdot_k1=pzdot
     
     dt_k1=dte_cycle/2.

     x3_k1 = x2e(m) + dt_k1*xdot_k1
     y3_k1 = y2e(m) + dt_k1*ydot_k1
     z3_k1 = z2e(m) + dt_k1*zdot_k1
     u3_k1 = u2e(m) + dt_k1*pzdot_k1
     mue3(m) = mue2(m)

     laps=anint((z3_k1/lz)-.5)*(1-peritr)
     r=x3_k1-0.5*lx+lr0
     i = int((r-rin)/dr)
     i = min(i,nr-1)
     i = max(i,0)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     y3_k1=modulo(y3_k1-laps*2*pi*qr*lr0/q0*sign(1.0,q0),ly)
     x3_k1=modulo(x3_k1,lx)
     z3_k1=modulo(z3_k1,lz)
!calculate k2

     x2_k0=x3_k1
     y2_k0=y3_k1
     z2_k0=z3_k1
     u2_k0=u3_k1
     mue2_k0=mue2(m)
   
     r=x2_k0-0.5*lx+lr0
     k = int(z2_k0/delz)
     wz0 = ((k+1)*delz-z2_k0)/delz
     wz1 = 1-wz0
     th = wz0*thfnz(k)+wz1*thfnz(k+1)

     i = int((r-rin)/dr)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     k = int((th+pi)/dth)
     wz0 = (-pi+(k+1)*dth-th)/dth
     wz1 = 1.-wz0
     dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
          +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1)
     dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
          +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1)
     grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
          +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1)
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
     radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
          +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1)
     dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
          +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1)
     qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
          +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1)
     grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
          +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1)
     gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
          +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1)

     curvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
          +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1)
     bdcrvbp = wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
          +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1)
     grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
          +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1)

     fp = wx0*f(i)+wx1*f(i+1)
     jfnp = wz0*jfn(k)+wz1*jfn(k+1)
     psipp = wx0*psip(i)+wx1*psip(i+1)
     psp = wx0*psi(i)+wx1*psi(i+1)
     ter = wx0*t0e(i)+wx1*t0e(i+1)
     kaptp = wx0*capte(i)+wx1*capte(i+1)
     kapnp = wx0*capne(i)+wx1*capne(i+1)
     xnp = wx0*xn0e(i)+wx1*xn0e(i+1)
     vncp = wx0*phincp(i)+wx1*phincp(i+1)
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)
     b=1.-tor+tor*bfldp
     pzp = emass*u2_k0/b+psp/br0
     xt = x2_k0
     yt = y2_k0
     zt = z2_k0
#ifdef JYCHENG
     include 'ppushlie.h'
#else
     i=int(xt/dx)
     j=int(yt/dy)
     k=int(zt/dz) !int(z3e(m)/dz)-gclr*kcnt

     wx0=real(i+1)-xt/dx
     wx1=1.-wx0
     wy0=real(j+1)-yt/dy
     wy1=1.-wy0
     wz0=real(k+1)-zt/dz
     wz1=1.-wz0
     x000=wx0*wy0*wz0
     x001=wx0*wy0*wz1
     x010=wx0*wy1*wz0
     x011=wx0*wy1*wz1
     x100=wx1*wy0*wz0
     x101=wx1*wy0*wz1
     x110=wx1*wy1*wz0
     x111=wx1*wy1*wz1
     exp1 = x000*ex_3d(i,j,k) + x100*ex_3d(i+1,j,k)  &
           + x010*ex_3d(i,j+1,k) + x110*ex_3d(i+1,j+1,k) + &
           x001*ex_3d(i,j,k+1) + x101*ex_3d(i+1,j,k+1) +   &
           x011*ex_3d(i,j+1,k+1) + x111*ex_3d(i+1,j+1,k+1)

     eyp = x000*ey_3d(i,j,k) + x100*ey_3d(i+1,j,k)  &
           + x010*ey_3d(i,j+1,k) + x110*ey_3d(i+1,j+1,k) + &
           x001*ey_3d(i,j,k+1) + x101*ey_3d(i+1,j,k+1) +   &
           x011*ey_3d(i,j+1,k+1) + x111*ey_3d(i+1,j+1,k+1)

     ezp = x000*ez_3d(i,j,k) + x100*ez_3d(i+1,j,k)  &
           + x010*ez_3d(i,j+1,k) + x110*ez_3d(i+1,j+1,k) + &
           x001*ez_3d(i,j,k+1) + x101*ez_3d(i+1,j,k+1) +   &
           x011*ez_3d(i,j+1,k+1) + x111*ez_3d(i+1,j+1,k+1)

#endif

     vfac = 0.5*(emass*u2_k0**2 + 2.*mue2_k0*b)
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw
     kap = kapnp - (1.5-vfac/ter)*kaptp

     ppar = u2_k0
     vpar = u2_k0-qel/emass*aparp*ipara
     bstar = b*(1+emass*vpar/(qel*b)*bdcrvbp)
     enerb=(mue2_k0+emass*vpar*vpar/b)/qel*b/bstar*tor
     dum2 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
     vxdum = (eyp/b+vpar/b*delbxp*iflut)*dum2
     xdot = vxdum*nonline  &
          -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
     ydot = (-exp1/b+vpar/b*delbyp*iflut)*dum2*nonline &
          +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
          (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
          +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
          -emass*vpar**2/(qel*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*grcgtp*lr0/q0*qhatp

     zdot =  vpar*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
          +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*q0*br0*grcgtp/jfnp

     pzd0 =tor*(-mue2_k0/emass/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
          +mue2_k0*vpar/(qel*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
     pzd1 = qel/emass*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
          +qel/emass*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)
     pzdot = pzd0+pzd1

     vpdum = vpar
     edot = qel*(xdot*exp1+(ydot-vp0)*eyp+zdot*ezp) &
          +qel*pzdot*aparp &
          +qel*vpdum*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)   &
          -qel*vpar*delbxp*vp0

     xdot_k2=xdot
     ydot_k2=ydot
     zdot_k2=zdot
     pzdot_k2=pzdot
     
     dt_k2=dte_cycle/2.

     x3_k2 = x2e(m) + dt_k2*xdot_k2
     y3_k2 = y2e(m) + dt_k2*ydot_k2
     z3_k2 = z2e(m) + dt_k2*zdot_k2
     u3_k2 = u2e(m) + dt_k2*pzdot_k2
     mue3(m) = mue2(m)

     laps=anint((z3_k2/lz)-.5)*(1-peritr)
     r=x3_k2-0.5*lx+lr0
     i = int((r-rin)/dr)
     i = min(i,nr-1)
     i = max(i,0)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     y3_k2=modulo(y3_k2-laps*2*pi*qr*lr0/q0*sign(1.0,q0),ly)

     x3_k2=modulo(x3_k2,lx)
     z3_k2=modulo(z3_k2,lz)
!calculate k3

     x2_k0=x3_k2
     y2_k0=y3_k2
     z2_k0=z3_k2
     u2_k0=u3_k2
     mue2_k0=mue2(m)
   
     r=x2_k0-0.5*lx+lr0
     k = int(z2_k0/delz)
     wz0 = ((k+1)*delz-z2_k0)/delz
     wz1 = 1-wz0
     th = wz0*thfnz(k)+wz1*thfnz(k+1)

     i = int((r-rin)/dr)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     k = int((th+pi)/dth)
     wz0 = (-pi+(k+1)*dth-th)/dth
     wz1 = 1.-wz0
     dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
          +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1)
     dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
          +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1)
     grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
          +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1)
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
     radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
          +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1)
     dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
          +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1)
     qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
          +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1)
     grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
          +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1)
     gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
          +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1)

     curvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
          +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1)
     bdcrvbp = wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
          +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1)
     grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
          +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1)

     fp = wx0*f(i)+wx1*f(i+1)
     jfnp = wz0*jfn(k)+wz1*jfn(k+1)
     psipp = wx0*psip(i)+wx1*psip(i+1)
     psp = wx0*psi(i)+wx1*psi(i+1)
     ter = wx0*t0e(i)+wx1*t0e(i+1)
     kaptp = wx0*capte(i)+wx1*capte(i+1)
     kapnp = wx0*capne(i)+wx1*capne(i+1)
     xnp = wx0*xn0e(i)+wx1*xn0e(i+1)
     vncp = wx0*phincp(i)+wx1*phincp(i+1)
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)
     b=1.-tor+tor*bfldp
     pzp = emass*u2_k0/b+psp/br0
     xt = x2_k0
     yt = y2_k0
     zt = z2_k0
#ifdef JYCHENG
     include 'ppushlie.h'
#else
     i=int(xt/dx)
     j=int(yt/dy)
     k=int(zt/dz) !int(z3e(m)/dz)-gclr*kcnt

     wx0=real(i+1)-xt/dx
     wx1=1.-wx0
     wy0=real(j+1)-yt/dy
     wy1=1.-wy0
     wz0=real(k+1)-zt/dz
     wz1=1.-wz0
     x000=wx0*wy0*wz0
     x001=wx0*wy0*wz1
     x010=wx0*wy1*wz0
     x011=wx0*wy1*wz1
     x100=wx1*wy0*wz0
     x101=wx1*wy0*wz1
     x110=wx1*wy1*wz0
     x111=wx1*wy1*wz1
     exp1 = x000*ex_3d(i,j,k) + x100*ex_3d(i+1,j,k)  &
           + x010*ex_3d(i,j+1,k) + x110*ex_3d(i+1,j+1,k) + &
           x001*ex_3d(i,j,k+1) + x101*ex_3d(i+1,j,k+1) +   &
           x011*ex_3d(i,j+1,k+1) + x111*ex_3d(i+1,j+1,k+1)

     eyp = x000*ey_3d(i,j,k) + x100*ey_3d(i+1,j,k)  &
           + x010*ey_3d(i,j+1,k) + x110*ey_3d(i+1,j+1,k) + &
           x001*ey_3d(i,j,k+1) + x101*ey_3d(i+1,j,k+1) +   &
           x011*ey_3d(i,j+1,k+1) + x111*ey_3d(i+1,j+1,k+1)

     ezp = x000*ez_3d(i,j,k) + x100*ez_3d(i+1,j,k)  &
           + x010*ez_3d(i,j+1,k) + x110*ez_3d(i+1,j+1,k) + &
           x001*ez_3d(i,j,k+1) + x101*ez_3d(i+1,j,k+1) +   &
           x011*ez_3d(i,j+1,k+1) + x111*ez_3d(i+1,j+1,k+1)

#endif

     vfac = 0.5*(emass*u2_k0**2 + 2.*mue2_k0*b)
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw
     kap = kapnp - (1.5-vfac/ter)*kaptp

     ppar = u2_k0
     vpar = u2_k0-qel/emass*aparp*ipara
     bstar = b*(1+emass*vpar/(qel*b)*bdcrvbp)
     enerb=(mue2_k0+emass*vpar*vpar/b)/qel*b/bstar*tor
     dum2 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
     vxdum = (eyp/b+vpar/b*delbxp*iflut)*dum2
     xdot = vxdum*nonline  &
          -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
     ydot = (-exp1/b+vpar/b*delbyp*iflut)*dum2*nonline &
          +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
          (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
          +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
          -emass*vpar**2/(qel*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*grcgtp*lr0/q0*qhatp

     zdot =  vpar*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
          +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*q0*br0*grcgtp/jfnp

     pzd0 =tor*(-mue2_k0/emass/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
          +mue2_k0*vpar/(qel*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
     pzd1 = qel/emass*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
          +qel/emass*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)
     pzdot = pzd0+pzd1

     vpdum = vpar
     edot = qel*(xdot*exp1+(ydot-vp0)*eyp+zdot*ezp) &
          +qel*pzdot*aparp &
          +qel*vpdum*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)   &
          -qel*vpar*delbxp*vp0

     xdot_k3=xdot
     ydot_k3=ydot
     zdot_k3=zdot
     pzdot_k3=pzdot
     
     dt_k3=dte_cycle

     x3_k3 = x2e(m) + dt_k3*xdot_k3
     y3_k3 = y2e(m) + dt_k3*ydot_k3
     z3_k3 = z2e(m) + dt_k3*zdot_k3
     u3_k3 = u2e(m) + dt_k3*pzdot_k3
     mue3(m) = mue2(m)


     laps=anint((z3_k3/lz)-.5)*(1-peritr)
     r=x3_k3-0.5*lx+lr0
     i = int((r-rin)/dr)
     i = min(i,nr-1)
     i = max(i,0)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     y3_k3=modulo(y3_k3-laps*2*pi*qr*lr0/q0*sign(1.0,q0),ly)
     x3_k3=modulo(x3_k3,lx)
     z3_k3=modulo(z3_k3,lz)
!calculate k4

     x2_k0=x3_k3
     y2_k0=y3_k3
     z2_k0=z3_k3
     u2_k0=u3_k3
     mue2_k0=mue2(m)
   
     r=x2_k0-0.5*lx+lr0
     k = int(z2_k0/delz)
     wz0 = ((k+1)*delz-z2_k0)/delz
     wz1 = 1-wz0
     th = wz0*thfnz(k)+wz1*thfnz(modulo(k+1,ntheta))

     i = int((r-rin)/dr)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     k = int((th+pi)/dth)
     wz0 = (-pi+(k+1)*dth-th)/dth
     wz1 = 1.-wz0
     dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
          +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1)
     dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
          +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1)
     grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
          +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1)
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
     radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
          +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1)
     dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
          +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1)
     qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
          +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1)
     grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
          +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1)
     gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
          +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1)

     curvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
          +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1)
     bdcrvbp = wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
          +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1)
     grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
          +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1)

     fp = wx0*f(i)+wx1*f(i+1)
     jfnp = wz0*jfn(k)+wz1*jfn(k+1)
     psipp = wx0*psip(i)+wx1*psip(i+1)
     psp = wx0*psi(i)+wx1*psi(i+1)
     ter = wx0*t0e(i)+wx1*t0e(i+1)
     kaptp = wx0*capte(i)+wx1*capte(i+1)
     kapnp = wx0*capne(i)+wx1*capne(i+1)
     xnp = wx0*xn0e(i)+wx1*xn0e(i+1)
     vncp = wx0*phincp(i)+wx1*phincp(i+1)
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)
     b=1.-tor+tor*bfldp
     pzp = emass*u2_k0/b+psp/br0
     xt = x2_k0
     yt = y2_k0
     zt = z2_k0
#ifdef JYCHENG
     include 'ppushlie.h'
#else
     i=int(xt/dx)
     j=int(yt/dy)
     k=int(zt/dz) !int(z3e(m)/dz)-gclr*kcnt

     wx0=real(i+1)-xt/dx
     wx1=1.-wx0
     wy0=real(j+1)-yt/dy
     wy1=1.-wy0
     wz0=real(k+1)-zt/dz
     wz1=1.-wz0
     x000=wx0*wy0*wz0
     x001=wx0*wy0*wz1
     x010=wx0*wy1*wz0
     x011=wx0*wy1*wz1
     x100=wx1*wy0*wz0
     x101=wx1*wy0*wz1
     x110=wx1*wy1*wz0
     x111=wx1*wy1*wz1
     exp1 = x000*ex_3d(i,j,k) + x100*ex_3d(i+1,j,k)  &
           + x010*ex_3d(i,j+1,k) + x110*ex_3d(i+1,j+1,k) + &
           x001*ex_3d(i,j,k+1) + x101*ex_3d(i+1,j,k+1) +   &
           x011*ex_3d(i,j+1,k+1) + x111*ex_3d(i+1,j+1,k+1)

     eyp = x000*ey_3d(i,j,k) + x100*ey_3d(i+1,j,k)  &
           + x010*ey_3d(i,j+1,k) + x110*ey_3d(i+1,j+1,k) + &
           x001*ey_3d(i,j,k+1) + x101*ey_3d(i+1,j,k+1) +   &
           x011*ey_3d(i,j+1,k+1) + x111*ey_3d(i+1,j+1,k+1)

     ezp = x000*ez_3d(i,j,k) + x100*ez_3d(i+1,j,k)  &
           + x010*ez_3d(i,j+1,k) + x110*ez_3d(i+1,j+1,k) + &
           x001*ez_3d(i,j,k+1) + x101*ez_3d(i+1,j,k+1) +   &
           x011*ez_3d(i,j+1,k+1) + x111*ez_3d(i+1,j+1,k+1)

#endif

     vfac = 0.5*(emass*u2_k0**2 + 2.*mue2_k0*b)
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw
     kap = kapnp - (1.5-vfac/ter)*kaptp

     ppar = u2_k0
     vpar = u2_k0-qel/emass*aparp*ipara
     bstar = b*(1+emass*vpar/(qel*b)*bdcrvbp)
     enerb=(mue2_k0+emass*vpar*vpar/b)/qel*b/bstar*tor
     dum2 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
     vxdum = (eyp/b+vpar/b*delbxp*iflut)*dum2
     xdot = vxdum*nonline  &
          -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
     ydot = (-exp1/b+vpar/b*delbyp*iflut)*dum2*nonline &
          +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
          (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
          +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
          -emass*vpar**2/(qel*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*grcgtp*lr0/q0*qhatp

     zdot =  vpar*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
          +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*q0*br0*grcgtp/jfnp

     pzd0 =tor*(-mue2_k0/emass/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
          +mue2_k0*vpar/(qel*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
     pzd1 = qel/emass*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
          +qel/emass*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)
     pzdot = pzd0+pzd1

     vpdum = vpar
     edot = qel*(xdot*exp1+(ydot-vp0)*eyp+zdot*ezp) &
          +qel*pzdot*aparp &
          +qel*vpdum*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)   &
          -qel*vpar*delbxp*vp0

     xdot_k4=xdot
     ydot_k4=ydot
     zdot_k4=zdot
     pzdot_k4=pzdot
     
     dt_k4=dte_cycle/6.

!fianl rk-4
     x3e(m) = x2e(m) + dt_k4*(xdot_k1+xdot_k4+2*xdot_k2+2*xdot_k3)
     y3e(m)  = y2e(m) + dt_k4*(ydot_k1+ydot_k4+2*ydot_k2+2*ydot_k3)
     z3e(m) = z2e(m) + dt_k4*(zdot_k1+zdot_k4+2*zdot_k2+2*zdot_k3)
     u3e(m) = u2e(m) + dt_k4*(pzdot_k1+pzdot_k4+2*pzdot_k2+2*pzdot_k3)
     mue3(m) = mue2(m)






     !         if(x<0.3)dum1=0.0

     dum = 1-w2e(m)*nonline*1.
     if(eldu.eq.1)dum = (tge/ter)**1.5*exp(vfac*(1/tge-1./ter))
     vxdum = (eyp/b+vpdum/b*delbxp)*dum2
#ifndef DIRECT_F
     w3e(m)=w2e(m) + dte_cycle*(  &
          (vxdum*kap + edot/ter-dum1*ppar*aparp/ter)*xnp     &
          +isg*(-dgdtp-zdot*dpdzp+xdot*exp1+ydot*eyp)/ter*xnp-tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp)*dum
#endif

333  continue
     laps=anint((z3e(m)/lz)-.5)*(1-peritr)
     r=x3e(m)-0.5*lx+lr0
     i = int((r-rin)/dr)
     i = min(i,nr-1)
     i = max(i,0)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     y3e(m)=modulo(y3e(m)-laps*2*pi*qr*lr0/q0*sign(1.0,q0)+8000.*ly,ly)
     if(x3e(m)>lx.and.iperidf==0)then
        x3e(m) = lx-1.e-8
        z3e(m)=lz-z3e(m)
        x2e(m) = x3e(m)
        z2e(m) = z3e(m)
        w2e(m) = 0.
        w3e(m) = 0.
     end if
     if(x3e(m)<0..and.iperidf==0)then
        x3e(m) = 1.e-8
        z3e(m)=lz-z3e(m)
        x2e(m) = x3e(m)
        z2e(m) = z3e(m)
        w2e(m) = 0.
        w3e(m) = 0.
     end if
     z3e(m)=modulo(z3e(m)+8.*lz,lz)
     x3e(m)=modulo(x3e(m)+800.*lx,lx)
     x3e(m) = min(x3e(m),lx-1.0e-8)
     y3e(m) = min(y3e(m),ly-1.0e-8)
     z3e(m) = min(z3e(m),lz-1.0e-8)

     u2e(m)=u3e(m)
     x2e(m)=x3e(m)
     y2e(m)=y3e(m)
     z2e(m)=z3e(m)
     w2e(m)=w3e(m)

  end do
  !$acc end parallel
  call MPI_ALLREDUCE(myopz,nopz,1,MPI_integer, &
       MPI_SUM, MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(myoen,noen,1,MPI_integer, &
       MPI_SUM, MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(myaven,aven,1,MPI_REAL8, &
       MPI_SUM, MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(myavptch,avptch,1,MPI_REAL8, &
       MPI_SUM, MPI_COMM_WORLD,ierr)
  aven = aven/(noen+0.1)
  avptch = avptch/(noen+0.1)

  np_old=mme
#ifndef OLD_PMOVE
  call test_init_pmove(z3e,np_old,lz,ierr)
  call test_pmove(x2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(x3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(y2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(y3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(z2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(z3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(u2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(u3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(w2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(w3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(mue2,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(mue3,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(index,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(ipass,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(mue,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(xie,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(z0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(eke,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(pze,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(u0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  

  call test_pmove(w0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
#ifdef DIRECT_F
  call test_pmove(w5e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call test_pmove(pif0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
#endif


#else
  call init_pmove(z3e,np_old,lz,ierr)
  call pmove(x2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(x3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(y2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(y3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(u2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(u3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w2e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w3e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(mue2,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(mue3,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(index,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(ipass,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(mue,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(xie,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(eke,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(pze,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(u0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit


  call pmove(w0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
#ifdef DIRECT_F
  call pmove(w5e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(pif0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
#endif
#endif
  call end_pmove(ierr)
  mme=np_new

  !      return
end subroutine pushe_rk4_1step
#endif


subroutine update_electron_weight
  use gem_com
  use gem_equil
  use mpi
!#ifdef ADIOS2
!  use coupling_core_edge
!  #include "adios_macro.h"
!#endif
  implicit none

  integer :: m,k,i,j,ierr
  integer, save :: init=0
  real :: r,wz0,wz1,th,wx0,wx1,xnp,ter,bfldp,b,vfac,wy0,wy1,&
             x000,x001,x010,x011,x100,x101,x110,x111,phip,new_f0,w4,dw,&
             xt,yt
  real :: pif0e_tmp(mme)
  

  !$acc parallel loop gang vector
  do m=1,mme

     !find the local phi
     r=x2e(m)-0.5*lx+lr0

     k = int(z2e(m)/delz)
     wz0 = ((k+1)*delz-z2e(m))/delz
     wz1 = 1-wz0
     th = wz0*thfnz(k)+wz1*thfnz(k+1)

     i = int((r-rin)/dr)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     k = int((th+pi)/dth)
     wz0 = (-pi+(k+1)*dth-th)/dth
     wz1 = 1.-wz0
     k=modulo(k,ntheta)

     xnp=wx0*xn0e(i)+wx1*xn0e(i+1)
     ter=wx0*t0e(i)+wx1*t0e(i+1)
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
     b=1.-tor+tor*bfldp
     vfac = 0.5*(emass*u2e(m)**2 + 2.*mue2(m)*b)
  
     xt=x2e(m)
     yt=y2e(m)  
     i=int(xt/dx)
     j=int(yt/dy)
     k=0 !int(z3e(m)/dz)-gclr*kcnt

     wx0=real(i+1)-xt/dx
     wx1=1.-wx0
     wy0=real(j+1)-yt/dy
     wy1=1.-wy0
     wz0=real(gclr*kcnt+k+1)-z2e(m)/dz
     wz1=1.-wz0
     x000=wx0*wy0*wz0
     x001=wx0*wy0*wz1
     x010=wx0*wy1*wz0
     x011=wx0*wy1*wz1
     x100=wx1*wy0*wz0
     x101=wx1*wy0*wz1
     x110=wx1*wy1*wz0
     x111=wx1*wy1*wz1

     phip = x000*phi(i,j,k) + x100*phi(i+1,j,k)  &
           + x010*phi(i,j+1,k) + x110*phi(i+1,j+1,k) + &
           x001*phi(i,j,k+1) + x101*phi(i+1,j,k+1) +   &
           x011*phi(i,j+1,k+1) + x111*phi(i+1,j+1,k+1)
     if(init==0 .and. iget==0)then
        pif0e(m)=xnp*sqrt(1./ter**3.)*exp(-vfac/ter) 
     endif
     new_f0=xnp*sqrt(1./ter**3.)*exp(-vfac/ter)

     !pif0e_tmp(m)=new_f0
     w4=1.0-new_f0/pif0e(m)*exp(phip/ter)
     !w4=1.0-exp(phip/ter)
     dw=w4-w5e(m)
     w3e(m)=w3e(m)+dw*w0e(m)
     w5e(m)=w4
     
  enddo
  !$acc end parallel

  init=1

  if(myid==0)then
    write(gemout,*)'xn0e=',xn0e
    write(gemout,*)'t0e=',t0e
    call flush(gemout)
    init=2 
  endif
  call mpi_barrier(mpi_comm_world,ierr)
end subroutine update_electron_weight

