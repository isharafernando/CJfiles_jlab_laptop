************************************************************************
*
* I.2 integer-string conversion
      subroutine inttochar(integ,car,min)
*     programmer: Alberto Accardi
*     date: 10/05/00
*
*  A. COMMENTARY
*
*     Converts an integer INTEG to a string CAR of length 9.
*     The number is written in the last 9-min positions of CAR, i.e
*     it s stored in       
*                     car(min:9)
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables
      
      integer integ, min, i, digit, zero, rest

      character car*9

      logical started, neg

*    *** constants
      integer max
      parameter(max=9)
      
*
*  C. ACTION
*

      zero = ichar('0')
      min=max

      neg = .false.
      if (integ.lt.0) then 
         neg = .true.
         integ=-integ
      end if

      started = .false.
      rest = 0
      do i = max, 1, -1
         digit = int(real(integ)/real(10**(I-1))-rest*10)
         rest = rest*10+digit
         car(max-i+1:max-i+1) = char(digit+zero)
         if ((.not.started).and.(digit.gt.0)) then
            started = .true.
            min = max-i+1
         end if
      end do
      if (neg) then
         min = min -1
         car(min:min) = '-'
      end if

      return
      end


