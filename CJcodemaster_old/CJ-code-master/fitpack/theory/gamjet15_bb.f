      subroutine d0gamjet(s,pt,nregion,ndrv,als,ans,scale)
      implicit double precision (a-h,o-z)
      dimension arg1(18), arg2(18), ygl(2), ygh(2), yjl(2), yjh(2),
     2yg(18), yj(18)
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
c      akfac1(x)=1.7052-.0025083*x+4.4876e-6*x**2-1.2827e-9*x**3
c      akfac2(x)=1.8628-.0044022*x+1.1747e-5*x**2-1.1336e-8*x**3
c      akfac3(x)=2.3962-.0130730*x+6.5069e-5*x**2+4.8698e-8*x**3
c      akfac4(x)=3.5016-.0330140*x+1.9910e-4*x**2-3.7492e-7*x**3
c  Updated using CJ12_min   5/12/14 -- jfo
      akfac1(x)=1.7175-.0036757*x+1.2366e-5*x**2-1.475e-8*x**3
      akfac2(x)=1.8551-.0045853*x+1.2705e-5*x**2-1.177e-8*x**3
      akfac3(x)=2.1176-.017025*x+1.3386e-4*x**2-3.0046e-7*x**3
      akfac4(x)=3.3173-.030816*x+1.7245e-4*x**2-2.7553e-7*x**3
c
      akfac1half(x) = 1.7424-0.0042506*x+6.118e-06*x**2+3.5531e-09*x**3
      akfac2half(x) = 1.9811-0.00788*x+2.489e-05*x**2-2.7256e-08*x**3
      akfac3half(x) = 2.5134-0.020312*x+0.00012151*x**2-1.1743e-07*x**3
      akfac4half(x) = 3.615-0.033696*x+0.00016777*x**2-2.3399e-07*x**3
      akfac1two(x) = 1.7284-0.0017206*x+3.6976e-06*x**2-1.7463e-09*x**3
      akfac2two(x) = 1.8538-0.0026899*x+3.481e-06*x**2+5.4257e-09*x**3
      akfac3two(x) = 2.4908-0.013999*x+0.00010096*x**2-1.0863e-07*x**3
      akfac4two(x) = 2.8335-0.010584*x-1.115e-05*x**2+2.471e-07*x**3
      if(nregion.eq.1)then
         ygl(1)=0.
         ygh(1)=1.
         yjl(1)=0.
         yjh(1)=.8
         ygl(2)=-1.
         ygh(2)=0.
         yjl(2)=-.8
         yjh(2)=0.
      else if(nregion.eq.2)then
         ygl(1)=0.
         ygh(1)=1.
         yjl(1)=-.8
         yjh(1)=0.
         ygl(2)=-1.
         ygh(2)=0.
         yjl(2)=0.
         yjh(2)=.8
      else if(nregion.eq.3)then
         ygl(1)=0.
         ygh(1)=1.
         yjl(1)=1.5
         yjh(1)=2.5
         ygl(2)=-1.
         ygh(2)=0.
         yjl(2)=-2.5
         yjh(2)=-1.5
      else if(nregion.eq.4)then
         ygl(1)=0.
         ygh(1)=1.
         yjl(1)=-2.5
         yjh(1)=-1.5
         ygl(2)=-1.
         ygh(2)=0.
         yjl(2)=1.5
         yjh(2)=2.5
      endif
      ans=0.
      do 100 j=1,2
      call gq11(ygl(j),ygh(j),1,yg,arg1,temp)
      do 200 l=1,6
      call gq11(yjl(j),yjh(j),1,yj,arg2,temp)
      do 210 k=1,6
      call dpcross(s,pt,yg(l),yj(k),ndrv,als,arg2(k),scale)
 210  continue
      call gq11(yjl(j),yjh(j),1,yj,arg2,temp)
      call gq11(yjl(j),yjh(j),0,yj,arg2,arg1(l))
 200  continue
      call gq11(ygl(j),ygh(j),1,yg,arg1,temp)
      call gq11(ygl(j),ygh(j),0,yg,arg1,temp)
      ans=ans+temp/(ygh(j)-ygl(j))/(yjh(j)-yjl(j))
  100  continue
c
c  Must average over the two results since we are averaging over the 
c  rapidity intervals. This makes it equivalent to the PSS Monte Carlo 
c  result.
c   
      ans=ans/2.
      if(IORD.eq.0)return
      if(nregion.eq.1)then
         if(scale.eq.1.d0)then
            ans=ans*akfac1(pt)
         else if(scale.eq.0.5d0)then
            ans=ans*akfac1half(pt)
         else if(scale.eq.2.d0)then
            ans=ans*akfac1two(pt)
         else
            print*,'error in scale choice)'
         endif
      else if(nregion.eq.2)then
         if(scale.eq.1.d0)then
            ans=ans*akfac2(pt)
         else if(scale.eq.0.5d0)then
            ans=ans*akfac2half(pt)
         else if(scale.eq.2.d0)then
            ans=ans*akfac2two(pt)
         else
            print*,'error in scale choice)'
         endif
      else if(nregion.eq.3)then
         if(scale.eq.1.d0)then
            ans=ans*akfac3(pt)
         else if(scale.eq.0.5d0)then
            ans=ans*akfac3half(pt)
         else if(scale.eq.2.d0)then
            ans=ans*akfac3two(pt)
         else
            print*,'error in scale choice)'
         endif
      else if(nregion.eq.4)then
         if(scale.eq.1.d0)then
            ans=ans*akfac4(pt)
         else if(scale.eq.0.5d0)then
            ans=ans*akfac4half(pt)
         else if(scale.eq.2.d0)then
            ans=ans*akfac4two(pt)
         else
            print*,'error in scale choice)'
         endif
      endif
      return
      end


      subroutine dpcross(s,pt,yg,yj,ndrv,als,ans,scale)
      implicit double precision (a-h,o-z)
      dimension cod(3), sigd(3), zc(18), arg(18), cob(13), sigb(13)
      common/curpar/xc(100)
c      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD 
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
      common/qs/q2,q2f,q2r
c      common/alphas/iloop,flavor,ilog,alam5
      common/isolate/rgam,ehadcut,iso,zcl,ifscale
      data pi/3.14159/
      data itime/0/
      qs=pt**2
      q2=qs*scale
      q2f=q2
      q2r=q2
      ifscale=0
      rs=sqrt(s)
      alam5=xc(1)
      alpha=1./137.
c      as=alpha_s(2,q2r,alam5,neff)
      as=alpha_s(IORD+1,q2r,alam5,neff)
      xa=pt/rs*(exp(yg)+exp(yj))
      if(xa.gt.1.d0)xa=0.9999d0
      xb=pt/rs*(exp(-yg)+exp(-yj))
      if(xb.gt.1.d0)xb=0.9999d0
      sh=xa*xb*s
      th=-pt*rs*xa*exp(-yg)
      uh=-pt*rs*xb*exp(yg)
      call coefd(xa,xb,ndrv,als,cod)
      call sigmad(sh,th,uh,sigd)
      facd=386.e6*pi*alpha*as/sh**2*2.*pt*xa*xb
      ansd=0.
      do 100 j=1,3
      ansd=ansd+facd*cod(j)*sigd(j)
      if(ansd.eq.0.d0)then
         print*,'direct',facd,cod(j),sigd(j),xa,xb
      endif
 100  continue
      sh0=sh
      th0=th
      uh0=uh
      xa0=xa
      xb0=xb
      n=1
      jmax=6*n
      zcmin=1./1.07
      call gq11(zcmin,1.d0,n,zc,arg,ansb)
      do 200 j=1,jmax
      xa=xa0/zc(j)
      xb=xb0/zc(j)
      if(xa.gt.1.d0)xa=0.9999d0
      if(xb.gt.1.d0)xb=0.9999d0
      sh=sh0/zc(j)**2
      th=th0/zc(j)**2
      uh=uh0/zc(j)**2
      call coefb(xa,xb,zc(j),ndrv,als,cob)
      call sigmab(sh,th,uh,sigb)
      facb=386.e6*pi*as**2/sh**2*2.*pt*xa*xb/zc(j)**2
      arg(j)=0.
      do 201 l=1,13
      arg(j)=arg(j)+facb*cob(l)*sigb(l)
c      if(arg(j).eq.0.d0)then
c         print*,'brem',facb,cob(l),sigb(l),xa,xb,zc(j)
c      endif
 201  continue
 200  continue
      call gq11(zcmin,1.d0,0,zc,arg,ansb)
      ans=ansd+ansb
      return
      end
      subroutine coefd(xa,xb,ndrv,als,cod)
      implicit double precision (a-h,o-z)
      dimension cod(3), qa(-5:5), qb(-5:5), eqs(-5:5)
      data eqs/1.d0,4.d0,1.d0,1.d0,4.d0,0.d0,4.d0,1.d0,1.d0,4.d0,1.d0/
      idum=1
      dum=0.d0
      call fsupdf(ndrv,xa,als,u,d,ub,db,sb,cb,bb,glue)
      qa(-5)=bb/xa
      qa(-4)=cb/xa
      qa(-3)=sb/xa
      qa(-2)=db/xa
      qa(-1)=ub/xa
      qa(0)=glue/xa
      qa(1)=u/xa
      qa(2)=d/xa
      qa(3)=sb/xa
      qa(4)=cb/xa
      qa(5)=bb/xa
c
c  Convert to pbar for D0
c 
      call fsupdf(ndrv,xb,als,u,d,ub,db,sb,cb,bb,glue)
      qb(-5)=bb/xb
      qb(-4)=cb/xb
      qb(-3)=sb/xb
      qb(-2)=d/xb
      qb(-1)=u/xb
      qb(0)=glue/xb
      qb(1)=ub/xb
      qb(2)=db/xb
      qb(3)=sb/xb
      qb(4)=cb/xb
      qb(5)=bb/xb
      do 100 j=1,3
      cod(j)=0.d0
 100  continue
      do 200 l=1,5
      cod(1)=cod(1)+eqs(l)*qa(0)*(qb(l)+qb(-l))/9.
      cod(2)=cod(2)+eqs(l)*qb(0)*(qa(l)+qa(-l))/9.
      cod(3)=cod(3)+eqs(l)*(qa(l)*qb(-l)+qa(-l)*qb(l))/9.
 200  continue
      if(cod(1).eq.0.d0.or.cod(2).eq.0.d0.or.cod(3).eq.0.d0)then
         print*,qa
         print*,qb
      endif
      return
      end
      subroutine sigmad(s,t,u,sigd)
      implicit double precision (a-h,o-z)
      dimension sigd(3)
      sigd(1)=-(s/u+u/s)/3.
      sigd(2)=-(s/t+t/s)/3.
      sigd(3)=(u/t+t/u)*8./9.
      return
      end
      subroutine coefb(xa,xb,z,ndrv,als,cob)
      implicit double precision (a-h,o-z)
      dimension cob(13), qa(-5:5), qb(-5:5), dc(-6:6)
      common/qs/q2,q2f,q2r
c      COMMON/ALPHAS/ILOOP,FLAVOR,ILOG,ALAM5
      COMMON/ISOLATE/rgam,EhadCUT,iso,ZCL,IFSCALE
      call fsupdf(ndrv,xa,als,u,d,ub,db,sb,cb,bb,glue)
      qa(-5)=bb/xa
      qa(-4)=cb/xa
      qa(-3)=sb/xa
      qa(-2)=db/xa
      qa(-1)=ub/xa
      qa(0)=glue/xa
      qa(1)=u/xa
      qa(2)=d/xa
      qa(3)=sb/xa
      qa(4)=cb/xa
      qa(5)=bb/xa
c
c  Convert to pbar for D0
c
      call fsupdf(ndrv,xb,als,u,d,ub,db,sb,cb,bb,glue)
      qb(-5)=bb/xb
      qb(-4)=cb/xb
      qb(-3)=sb/xb
      qb(-2)=d/xb
      qb(-1)=u/xb
      qb(0)=glue/xb
      qb(1)=ub/xb
      qb(2)=db/xb
      qb(3)=sb/xb
      qb(4)=cb/xb
      qb(5)=bb/xb
      call dkgam(5,z,dc)
      do 100 j=1,13
      cob(j)=0.d0
 100  continue
      do 200 k=1,5
      do 200 l=1,5
      if(k.eq.l) goto 200
      cob(1)=cob(1)+qa(k)*qb(l)*dc(k)+qa(-k)*qb(-l)*dc(-k)
      cob(2)=cob(2)+qa(k)*qb(l)*dc(l)+qa(-k)*qb(-l)*dc(-l)
      cob(3)=cob(3)+qa(k)*qb(-l)*dc(k)+qa(-k)*qb(l)*dc(-k)
      cob(4)=cob(4)+qa(k)*qb(-l)*dc(-l)+qa(-k)*qb(l)*dc(l)
 200  continue
      do 210 l=1,5
      cob(5)=cob(5)+qa(l)*qb(l)*dc(l)+qa(-l)*qb(-l)*dc(-l)
      cob(7)=cob(7)+qa(l)*qb(-l)*dc(l)+qa(-l)*qb(l)*dc(-l)
      cob(8)=cob(8)+qa(l)*qb(-l)*dc(-l)+qa(-l)*qb(l)*dc(l)
      cob(9)=cob(9)+dc(0)*(qa(l)*qb(-l)+qa(-l)*qb(l))
      cob(10)=cob(10)+qa(0)*qb(0)*(dc(l)+dc(-l))
      cob(11)=cob(11)+qb(0)*(qa(l)*dc(l)+qa(-l)*dc(-l))
     2               +qa(0)*dc(0)*(qb(l)+qb(-l))
      cob(12)=cob(12)+qb(0)*dc(0)*(qa(l)+qa(-l))
     2               +qa(0)*(qb(l)*dc(l)+qb(-l)*dc(-l)) 
 210  continue
      cob(13)=qa(0)*qb(0)*dc(0)
      return
      end
      subroutine sigmab(s,t,u,sigb)
      implicit double precision (a-h,o-z)
      dimension sigb(13)
      s2=s*s
      t2=t*t
      u2=u*u
      sigb(1)=4./9.*(s2+u2)/t2
      sigb(2)=4./9.*(s2+t2)/u2
      sigb(3)=4./9.*(s2+u2)/t2
      sigb(4)=4./9.*(s2+t2)/u2
      sigb(5)=4./9.*((s2+u2)/t2+(s2+t2)/u2)-8./27.*s2/t/u
      sigb(6)=4./9.*(t2+u2)/s2
      sigb(7)=4./9.*((s2+u2)/t2+(u2+t2)/s2)-8./27.*u2/s/t
      sigb(8)=4./9.*((s2+t2)/u2+(u2+t2)/s2)-8./27.*t2/s/u
      sigb(9)=32./27.*(t/u+u/t)-8./3.*(t2+u2)/s2
      sigb(10)=1./6.*(t/u+u/t)-3./8.*(t2+u2)/s2
      sigb(11)=-4./9.*(s/u+u/s)+(s2+u2)/t2
      sigb(12)=-4./9.*(s/t+t/s)+(s2+t2)/u2
      sigb(13)=9./2.*(3.-t*u/s2-s*u/t2-s*t/u2)
      return
      end
