! A helper to read in the tables
program rw

  integer,parameter:: nx=91,nq=111
    open(1,file="f2_tm1ht1th0os1.dat")
      do iq=1,nq
        do ix=1,nx
           read(1,*) x,q2,fp,fn,fd
!           write(2,'(f4.2,2x,f4.1,3(1pe12.4))') x,q2,fp,fn,fd
           write(2,'(f4.2,2x,f4.1,3(1pe12.4))') x,q2,fp,fn,fd
        end do
      end do
    close(1)

end program
