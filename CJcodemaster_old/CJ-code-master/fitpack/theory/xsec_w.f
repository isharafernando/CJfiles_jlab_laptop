      subroutine xsec_w(ibeam,iord,ich,rs,ans)
      implicit double precision (a-h,o-z)
      dimension u1(60),arg0(60),u2(60),arg1(60),ans1(60)
      s=rs**2
      GF=1.16637e-05
      xmw=80.4
      pi=3.14159
c      s2thw=0.2231
      s2thw=0.232
c      cu2=0.25+(0.5-4./3.*s2thw)**2
c      cd2=0.25+(0.5-2./3.*s2thw)**2
      BR=10.8e-02
      sig0=BR*pi/3.*sqrt(2.)*GF*xmw**2/s*389.4e03
      tau=xmw**2/s
c      call setctq6(200)
      as=0.
      if(iord.eq.2)as=alpha_s(2,xmw**2,.226d0,neff)/2./pi
      n=6
      npt=6*n
      xl=log(tau)
      call gq11(xl,0.d0,n,u1,arg0,ans)
      do j=1,npt
         x1=exp(u1(j))
         x20=tau/x1
         call pdffac(ibeam,ich,x1,x20,xmw,qq,qg)
         arg0(j)=sig0*qq*(1.
     2                 +as*8.*.5*4./3.*(log(1.-tau/x1))**2
     2                 +as*4./3.*(2./3.*pi**2-8.))
         u2l=log(tau/x1)
         call gq11(u2l,0.d0,n,u2,arg1,ans1(j))
         do k=1,npt
            zz=exp(u2(k))
            x2=tau/x1/zz
            call pdffac(ibeam,ich,x1,x2,xmw,qq1,qg1)
            arg1(k)=sig0*as*4./3.*(qq1*4.*(1.+zz**2)-
     2      zz*qq*8.)*(log(1.-zz))/(1.-zz)
            arg1(k)=arg1(k)+sig0*as*qq1*dq(zz)
            arg1(k)=arg1(k)+sig0*as*qg1*dg(zz)
         enddo
         call gq11(u2l,0.d0,0,u2,arg1,ans1(j))
         arg0(j)=arg0(j)+ans1(j)
      enddo
      call gq11(xl,0.d0,n,u1,arg0,ans)
      call gq11(xl,0.d0,0,u1,arg0,ans)
      return
      end
      subroutine pdffac(ibeam,ich,x1,x2,q,qq,qg)
      implicit double precision (a-h,o-z)
      character prefix*50
      dimension q1(-5:5),q2(-5:5)
      data itime/0/
      COMMON/CURPAR/CURPAR(100)  
      COMMON/Q2STUF/ Q02,Q2MAX
      vud2=.97419**2
      vus2=.2257**2
      vub2=.00359**2
      vcd2=.2256**2
      vcs2=.9733**2
      vcb2=.0415**2      
      als=log(log(q**2/curpar(1)**2)/log(q02/curpar(1)**2))
      call fsupdf(0,x1,als,u,d,ub,db,sb,cb,bb,glue)
      q1(-5)=bb/x1
      q1(-4)=cb/x1
      q1(-3)=sb/x1
      q1(-2)=db/x1
      q1(-1)=ub/x1
      q1(0)=glue/x1
      q1(1)=u/x1
      q1(2)=d/x1
      q1(3)=q1(-3)
      q1(4)=q1(-4)
      q1(5)=q1(-5)
      call fsupdf(0,x2,als,u,d,ub,db,sb,cb,bb,glue)
      q2(-5)=bb/x2
      q2(-4)=cb/x2
      q2(-3)=sb/x2
      q2(-2)=db/x2
      q2(-1)=ub/x2
      q2(0)=glue/x2
      q2(1)=u/x2
      q2(2)=d/x2
      q2(3)=q2(-3)
      q2(4)=q2(-4)
      q2(5)=q2(-5)
      if(ibeam.eq.2)then
         dum=q2(2)
         q2(2)=q2(-2)
         q2(-2)=dum
         dum=q2(1)
         q2(1)=q2(-1)
         q2(-1)=dum
      endif
      if(ich.eq.1)then
         qq=vud2*(q1(1)*q2(-2)+q1(-2)*q2(1))
     2     +vus2*(q1(1)*q2(-3)+q1(-3)*q2(1))
     3     +vcd2*(q1(4)*q2(-2)+q1(-2)*q2(4))
     4     +vcs2*(q1(4)*q2(-3)+q1(-3)*q2(4))
         qg=q1(0)*(vud2*(q2(1)+q2(-2))+vus2*(q2(1)+q2(-3))
     2            +vcd2*(q2(4)+q2(-2))+vcs2*(q2(4)+q2(-3)))
     3     +q2(0)*(vud2*(q1(1)+q1(-2))+vus2*(q1(1)+q1(-3))
     4            +vcd2*(q1(4)+q1(-2))+vcs2*(q1(4)+q1(-3)))
      else if(ich.eq.-1)then
         qq=vud2*(q1(-1)*q2(2)+q1(2)*q2(-1))
     2     +vus2*(q1(-1)*q2(3)+q1(3)*q2(-1))
     3     +vcd2*(q1(-4)*q2(2)+q1(2)*q2(-4))
     4     +vcs2*(q1(-4)*q2(3)+q1(3)*q2(-4))
         qg=q1(0)*(vud2*(q2(-1)+q2(2))+vus2*(q2(-1)+q2(3))
     2            +vcd2*(q2(-4)+q2(2))+vcs2*(q2(-4)+q2(3)))
     3     +q2(0)*(vud2*(q1(-1)+q1(2))+vus2*(q1(-1)+q1(3))
     4            +vcd2*(q1(-4)+q1(2))+vcs2*(q1(-4)+q1(3)))
      endif
      return
      end
      function dg(z)
      implicit double precision (a-h,o-z)
      dg=0.5*((z**2+(1.-z)**2)*log((1.-z)**2/z)+1./2.+3.*z-3.5*z**2)
      return
      end
      function dq(z)
      implicit double precision (a-h,o-z)
      dq=-2.*(1.+z**2)/(1.-z)*log(z)
      dq=dq*4./3.
      return
      end  
