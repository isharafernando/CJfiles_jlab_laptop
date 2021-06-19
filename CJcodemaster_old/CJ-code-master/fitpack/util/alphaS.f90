!-----------------------------------------------------------
! Alpha_S calculation as in PEGASUS
! by: N.Sato ca. 2016
!     (some mods to adapt to CJ code by A.Accardi Feb 2018)
!-----------------------------------------------------------

module params
implicit none
real*8,parameter::pi=3.14159265359d0
real*8,parameter::mz2=91.187D0**2
real*8 :: mc2=1.69D0
real*8 :: mb2=20.25
real*8,parameter::mt2=172.9D0**2
real*8 :: alphaS_old=-1d0  ! initial value to force self-initialization
real*8 :: alphaSMZ=0.118044d0 ! default PDG value
real*8,parameter::alphaS0=0.371912360748d0
real*8,dimension(7,0:3)::beta
integer::iorder=1  ! order of evolution 0:LO 1:NLO 2:NNLO
real*8::Q20,a0,ac,ab,at,az
end module params

function get_Nf(Q2)
  use params
  implicit none
  integer::Nf,get_Nf
  real*8::Q2
  Nf=3
  if (Q2>=mc2) Nf=Nf+1
  if (Q2>=mb2) Nf=Nf+1
  if (Q2>=mt2) Nf=Nf+1
  get_Nf=Nf
end function get_Nf

function get_betaf(a,Nf)
  use params
  implicit none
  integer::Nf
  real*8::betaf,get_betaf,a
  betaf = -beta(Nf,0)
  if (iorder>=1)  betaf=betaf-a*beta(Nf,1)
  if (iorder>=2)  betaf=betaf-a*beta(Nf,2)
  get_betaf=betaf*a**2
end function get_betaf

function evolve_a(Q2ini,a,Q2,Nf)
! Runge-Kutta implemented in pegasus  
  use params
  implicit none
  real*8::evolve_a,Q2ini,a,Q2,LR,a_
  integer::Nf,k
  real*8::XK0,XK1,XK2,XK3
  real*8::get_betaf
  LR = dlog(Q2/Q2ini)/20D0
  a_=a
  do k=1,20
     XK0 = LR * get_betaf(a_,Nf)
     XK1 = LR * get_betaf(a_ + 0.5 * XK0,Nf)
     XK2 = LR * get_betaf(a_ + 0.5 * XK1,Nf)
     XK3 = LR * get_betaf(a_ + XK2,Nf)
     a_=a_ + (XK0 + 2.* XK1 + 2.* XK2 + XK3) * 0.166666666666666
  end do
  evolve_a=a_
end function evolve_a

function get_a(Q2)
  use params
  implicit none
  real*8::get_a,Q2ini,aini,Q2,a
  real*8::evolve_a
  integer:: Nf
  if (mb2<=Q2) then 
     Q2ini = mb2
     aini  = ab
     Nf    = 5
  else if (mc2<=Q2 .and. Q2<mb2) then 
     Q2ini = mc2
     aini  = ac
     Nf    = 4
  else if (Q2<mc2) then 
     Q2ini = Q20
     aini  = a0
     Nf    = 3
  end if
  a=evolve_a(Q2ini,aini,Q2,Nf)
  get_a=a
end function get_a


subroutine set_alphaSMZ(asmz)
  use params
  implicit none
  real*8::asmz
  alphaSMZ = asmz 
end subroutine set_alphaSMZ

subroutine set_alphaS_pars(iord,xmc,xmb)
  use params
  implicit none
  integer iord
  real*8::xmc,xmb
  iorder = iord
  mc2 = xmc**2
  mb2 = xmb**2
end subroutine set_alphaS_pars


function get_alphaS(Q2)
  use params
  implicit none
  real*8::get_alphaS,get_a,Q2

  if (alphas_old.ne.alphaSMZ) then
     call setup_aS()
  end if
  get_alphaS = get_a(Q2)*4d0*pi
end function get_alphaS


function get_lambdaQCD()
  use params
  ! returns the 5 flavors LambdaQCD using Eq.(3.22) of
  ! Deur, Brodski, de Teramond, Prog.Part.Nucl.Phys.90(2016)1 
  implicit none
  real*8::get_lambdaQCD,get_alphaS
  real*8::l2,as ! Lambda_QCD^2, alpha_s(MZ)
  integer::Nf=5 ! 5 flavors
  real*8::b0,b1 ! beta_0 and beta_1 (5 flavors)
  real*8::c0,c1

  b0 = beta(Nf,0)
  as = get_alphaS(mz2)
  c0 = 4*pi/(b0*as)
  if (iorder==0) then
     l2 = mz2 * dexp(-c0)
  else if (iorder==1) then
     b1 = beta(Nf,1)
     c1 = b1/(4*pi*b0)
     l2 = mz2 * dexp(-c0) * (c0 + c1)**c1
  else
     print*, 'ERROR(get_lambdaQCD): not immplemented for NNLO and beyond'
  end if
  
  get_lambdaQCD = l2**0.5
  
end function get_lambdaQCD


subroutine setup_aS()
  use params
  implicit none
  integer::Nf
  real*8::evolve_a
  do Nf=3,7
     beta(Nf,0)=11D0-2D0/3D0*Nf
     beta(Nf,1)=102D0-38D0/3D0*Nf 
     beta(Nf,2)=2857D0/2.0-5033D0/18D0*Nf+325D0/54D0*Nf**2 
  enddo
  Q20=1D0
  az=alphaSMZ/(4D0*pi)
  ab=evolve_a(mz2,az,mb2,5)
  ac=evolve_a(mb2,ab,mc2,4)
  a0=evolve_a(mc2,ac,Q20,3)
  
  !Q20=mc2
  !ac=alphaS0
  !ab=evolve_a(mc2,ac,mb2,4)

end subroutine setup_aS




















