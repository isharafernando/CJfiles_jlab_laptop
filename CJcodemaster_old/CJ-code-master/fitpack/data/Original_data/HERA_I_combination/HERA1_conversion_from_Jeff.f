      program datfile
c      character*124 title
      dimension dum(124)
      open(unit=15,file='tmp.dat', status='unknown')
      open(unit=16,file='HERA_cc_el_cor',status='unknown')
c      n=0
c      read(15,*,end=999)title
 200  continue
      read(15,*,end=999) (dum(i),i=1,124)
      dum(7)=dum(7)*dum(5)/100.
      dum(10)=dum(10)*dum(5)/100.
      cor=sqrt(dum(10)**2-dum(7)**2)
      unc=dum(8)*dum(5)/100.
      print*,dum(1)
      write(16,998)dum(3),dum(2),dum(4),dum(5),dum(7),
     2cor,unc,dum(11),dum(12),dum(13),(dum(j),j=15,124)
 998  format(3e10.3,117e12.4)
      goto 200
 999  continue
c      print*,n
      call exit
      end
