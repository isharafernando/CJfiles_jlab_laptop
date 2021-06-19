      SUBROUTINE DKGAM(IBREM,Z,DC)
      implicit double precision (a-h,o-z)
      DIMENSION DC(-6:6)
      COMMON/QS/Q2,q2f,q2r
      COMMON/ALPHAS/ILOOP,FLAVOR,ILOG,ALAM5
      COMMON/ISOLATE/rgam,EhadCUT,iso,ZCL,IFSCALE
      common/invar/xa,xb,zc,v,vol
c
c  Modified 1/24/02 to use BFG Sets 1 and 2 as options
c
      DV=0.
      DS=0.
      DG=0.
      GO TO (6,7,8,11,12),IBREM
    6 CONTINUE
      DO 10 LF=-6,6
   10 DC(LF)=0.
      RETURN
    7 CONTINUE
      DV=(1.+(1.-Z)**2)/Z
      GO TO 9
    8 CONTINUE
      DV=(2.21-1.28*Z+1.29*Z**2)/(1.-1.63*dLOG(1.-Z))
     2/Z**.951
      DS=.002*(1.-Z)**2/Z**2.54
      DG=.194/8.*(1.-Z)**1.03/Z**1.97
    9 CONTINUE
      QM=IFSCALE*RGAM**2*Q2f+(1-IFSCALE)*Q2f
      ALQ=dLOG(QM/ALAM5**2)
      FAC=ALQ/(137.*2.*3.14159)
      D1=FAC*(4./9.*DV+DS)
      D2=FAC*(DV/9.+DS)
      DG=DG*FAC
      DC(1)=D1
      DC(2)=D2
      DC(3)=D2
      DC(4)=D1
      DC(5)=0.
      DC(6)=0.
      DC(-1)=D1
      DC(-2)=D2
      DC(-3)=D2
      DC(-4)=D1
      DC(-5)=0.
      DC(-6)=0.
      DC(0)=DG
      RETURN
 11   continue
      QM=IFSCALE*RGAM**2*Q2f+(1-IFSCALE)*Q2f
      call fonfra(z,1,qm,g,u,d,s,c,b)
      goto 13
 12   continue
      QM=IFSCALE*RGAM**2*Q2f+(1-IFSCALE)*Q2f
      call fonfra(z,2,qm,g,u,d,s,c,b)
 13   continue
      dc(1)=u/z
      dc(2)=d/z
      dc(3)=s/z
      dc(4)=c/z
      dc(5)=b/z
      dc(6)=0.d0
      dc(0)=g/z
      dc(-1)=dc(1)
      dc(-2)=dc(2)
      dc(-3)=dc(3)
      dc(-4)=dc(4)
      dc(-5)=dc(5)
      dc(-6)=dc(6)
      return
      end
