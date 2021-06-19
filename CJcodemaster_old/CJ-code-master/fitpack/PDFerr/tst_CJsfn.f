      Program tst_CJpdf
      implicit real*8 (a-h,o-z)
      dimension q(4), sfn(1:9)
      data q/1.3, 10.,31.6228,100./
      print*,'ISET = '
      read*,iset
      call setCJ(iset)
      open(unit=1,file='tst_CJsfn.out', status='unknown')
      do j=1,4
         write(1,2) q(j)**2
         write(1,5)
 5       format(' Structure functions')
 2       format(' Q**2=',f10.2,' GeV**2')
         write(1,4)
 4       format(5x,' x',8x,'F2p',9x,'F2n',9x,'F2d'
     &        ,8x,'F2p(LT)',5x,'F2n(LT)',5x,'F2d(LT)'
     &        ,5x,'F2p(0)',6x,'F2n(0)',6x,'F2d(0)'
     &        )
         do l=1,19
            x=l*.05
            do k=1,9
               sfn(k)=CJsfn(k,x,q(j))
            enddo
            write(1,3) x, sfn
 3          format(f10.3, 11e12.3)
         enddo
      enddo
      print*, 'Output PDFs in tst_CJsfn.out'
      call exit
      end
