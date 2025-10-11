function cdbonnpot(lp,l,kp,k,jj,S,Tz)
use cdbonnmodule, only : cdbonn
implicit none
real*8 :: cdbonnpot
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call cdbonn

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        cdbonnpot=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        cdbonnpot=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        cdbonnpot=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        cdbonnpot=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        cdbonnpot=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        cdbonnpot=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        cdbonnpot=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        cdbonnpot=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function loemn450(lp,l,kp,k,jj,S,Tz)
use lo450module, only : lo450
implicit none
real*8 :: loemn450
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call lo450

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        loemn450=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        loemn450=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        loemn450=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        loemn450=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        loemn450=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        loemn450=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        loemn450=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        loemn450=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function loemn500(lp,l,kp,k,jj,S,Tz)
use lo500module, only : lo500
implicit none
real*8 :: loemn500
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call lo500

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        loemn500=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        loemn500=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        loemn500=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        loemn500=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        loemn500=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        loemn500=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        loemn500=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        loemn500=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function loemn550(lp,l,kp,k,jj,S,Tz)
use lo550module, only : lo550
implicit none
real*8 :: loemn550
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call lo550

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        loemn550=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        loemn550=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        loemn550=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        loemn550=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        loemn550=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        loemn550=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        loemn550=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        loemn550=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function nloemn450(lp,l,kp,k,jj,S,Tz)
use nlo450module, only : nlo450
implicit none
real*8 :: nloemn450
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call nlo450

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        nloemn450=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        nloemn450=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        nloemn450=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        nloemn450=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        nloemn450=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        nloemn450=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        nloemn450=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        nloemn450=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function nloemn500(lp,l,kp,k,jj,S,Tz)
use nlo500module, only : nlo500
implicit none
real*8 :: nloemn500
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call nlo500

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        nloemn500=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        nloemn500=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        nloemn500=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        nloemn500=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        nloemn500=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        nloemn500=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        nloemn500=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        nloemn500=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function nloemn550(lp,l,kp,k,jj,S,Tz)
use nlo550module, only : nlo550
implicit none
real*8 :: nloemn550
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call nlo550

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        nloemn550=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        nloemn550=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        nloemn550=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        nloemn550=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        nloemn550=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        nloemn550=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        nloemn550=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        nloemn550=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function n2loemn450(lp,l,kp,k,jj,S,Tz)
use n2lo450module, only : n2lo450
implicit none
real*8 :: n2loemn450
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call n2lo450

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        n2loemn450=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        n2loemn450=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        n2loemn450=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        n2loemn450=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        n2loemn450=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        n2loemn450=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        n2loemn450=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        n2loemn450=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function n2loemn500(lp,l,kp,k,jj,S,Tz)
use n2lo500module, only : n2lo500
implicit none
real*8 :: n2loemn500
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call n2lo500

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        n2loemn500=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        n2loemn500=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        n2loemn500=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        n2loemn500=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        n2loemn500=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        n2loemn500=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        n2loemn500=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        n2loemn500=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function n2loemn550(lp,l,kp,k,jj,S,Tz)
use n2lo550module, only : n2lo550
implicit none
real*8 :: n2loemn550
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call n2lo550

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        n2loemn550=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        n2loemn550=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        n2loemn550=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        n2loemn550=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        n2loemn550=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        n2loemn550=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        n2loemn550=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        n2loemn550=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function n3loemn450(lp,l,kp,k,jj,S,Tz)
use n3lo450newmodule, only : n3lo450new
implicit none
real*8 :: n3loemn450
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call n3lo450new

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        n3loemn450=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        n3loemn450=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        n3loemn450=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        n3loemn450=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        n3loemn450=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        n3loemn450=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        n3loemn450=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        n3loemn450=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function n3loemn500(lp,l,kp,k,jj,S,Tz)
use n3lo500newmodule, only : n3lo500new
implicit none
real*8 :: n3loemn500
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call n3lo500new

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        n3loemn500=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        n3loemn500=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        n3loemn500=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        n3loemn500=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        n3loemn500=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        n3loemn500=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        n3loemn500=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        n3loemn500=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function n3loemn550(lp,l,kp,k,jj,S,Tz)
use n3lo550newmodule, only : n3lo550new
implicit none
real*8 :: n3loemn550
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call n3lo550new

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        n3loemn550=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        n3loemn550=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        n3loemn550=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        n3loemn550=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        n3loemn550=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        n3loemn550=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        n3loemn550=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        n3loemn550=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function n4loemn450(lp,l,kp,k,jj,S,Tz)
use n4lo450module, only : n4lo450
implicit none
real*8 :: n4loemn450
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call n4lo450

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        n4loemn450=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        n4loemn450=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        n4loemn450=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        n4loemn450=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        n4loemn450=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        n4loemn450=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        n4loemn450=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        n4loemn450=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function n4loemn500(lp,l,kp,k,jj,S,Tz)
use n4lo500module, only : n4lo500
implicit none
real*8 :: n4loemn500
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call n4lo500

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        n4loemn500=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        n4loemn500=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        n4loemn500=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        n4loemn500=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        n4loemn500=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        n4loemn500=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        n4loemn500=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        n4loemn500=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function n4loemn550(lp,l,kp,k,jj,S,Tz)
use n4lo550module, only : n4lo550
implicit none
real*8 :: n4loemn550
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call n4lo550

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        n4loemn550=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        n4loemn550=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        n4loemn550=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        n4loemn550=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        n4loemn550=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        n4loemn550=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        n4loemn550=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        n4loemn550=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function nnloopt(lp,l,kp,k,jj,S,Tz)
use nnlo_optmodule, only : nnlo_opt
implicit none
real*8 :: nnloopt
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call nnlo_opt

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        nnloopt=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        nnloopt=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        nnloopt=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        nnloopt=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        nnloopt=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        nnloopt=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        nnloopt=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        nnloopt=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function nnlosat(lp,l,kp,k,jj,S,Tz)
use nnlo_satmodule, only : nnlo_sat
implicit none
real*8 :: nnlosat
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call nnlo_sat

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        nnlosat=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        nnlosat=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        nnlosat=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        nnlosat=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        nnlosat=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        nnlosat=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        nnlosat=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        nnlosat=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function n3loem(lp,l,kp,k,jj,S,Tz)
use n3lomodule, only : n3lo
implicit none
real*8 :: n3loem
real*8,intent(in) :: kp,k
integer,intent(in) :: lp,l,jj,S,Tz

integer :: j,inn
real*8 :: v,xmev,ymev
logical :: heform,sing,trip,coup
common /cpot/   v(6),xmev,ymev
common /cstate/ j,heform,sing,trip,coup
common /cnn/ inn

heform=.false.
sing=.true.
trip=.true.
coup=.true.

j=jj
inn=Tz+2
xmev=kp
ymev=k

call n3lo

if(j .eq. 0)then
    if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
        n3loem=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
        n3loem=v(3)
    else
        write(*,*) "error1"
        return
    end if
else if(j .gt. 0)then
    if( (S .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
        n3loem=v(1)
    else if(( (S .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
        n3loem=v(2)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
        n3loem=v(3)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
        n3loem=v(4)
    else if(( (S .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
        n3loem=v(5)
    else if(( (S .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
        n3loem=v(6)
    else
        write(*,*) "error2"
        return
    end if
else
    write(*,*) "error3"
    return
end if

return
end function


function idahopot(lpot,lp,l,s,j,Tz,r)
    use idahopwmodule, only : idahopw
    implicit none
    real*8 :: idahopot
    integer,intent(in) :: lpot,lp,l,j,s,Tz
    real*8,intent(in) :: r
    integer :: t,t1z,t2z
    real*8 :: vpw(2,2)

    if(Tz .eq. -1)then
        t1z=1
        t2z=1
    else if(Tz .eq. 1)then
        t1z=-1
        t2z=-1
    else if(Tz .eq. 0)then
        t1z=-1
        t2z=1
    else
        write(*,*) "error"
        return
    end if

    t = mod(l+s+1,2)

    call idahopw(lpot,l,s,j,t,t1z,t2z,r,vpw)

    if(j .eq. 0)then
        if( (s .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
            idahopot=vpw(1,1)
        else if(( (s .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
            idahopot=vpw(1,1)
        else
            write(*,*) "error"
            return
        end if
    else if(j .gt. 0)then
        if( (s .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
            idahopot=vpw(1,1)
        else if(( (s .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
            idahopot=vpw(1,1)
        else if(( (s .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
            idahopot=vpw(2,2)
        else if(( (s .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
            idahopot=vpw(1,1)
        else if(( (s .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
            idahopot=vpw(2,1)
        else if(( (s .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
            idahopot=vpw(1,2)
        else
            write(*,*) "error"
            return
        end if
    else
        write(*,*) "error"
        return
    end if

    return
end function


function av18pot(lpot,lp,l,s,j,Tz,r)
    use av18pwmodule, only : av18pw
    implicit none
    real*8 :: av18pot
    integer,intent(in) :: lpot,lp,l,j,s,Tz
    real*8,intent(in) :: r
    integer :: t,t1z,t2z
    real*8 :: vpw(2,2)

    if(Tz .eq. -1)then
        t1z=1
        t2z=1
    else if(Tz .eq. 1)then
        t1z=-1
        t2z=-1
    else if(Tz .eq. 0)then
        t1z=-1
        t2z=1
    else
        write(*,*) "error"
        return
    end if

    t = mod(l+s+1,2)

    call av18pw(lpot,l,s,j,t,t1z,t2z,r,vpw)

    if(j .eq. 0)then
        if( (s .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
            av18pot=vpw(1,1)
        else if(( (s .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
            av18pot=vpw(1,1)
        else
            write(*,*) "error"
            return
        end if
    else if(j .gt. 0)then
        if( (s .eq. 0) .and. (lp .eq. j) .and. (l .eq. j) ) then
            av18pot=vpw(1,1)
        else if(( (s .eq. 1) .and. (lp .eq. j) .and. (l .eq. j) ))then
            av18pot=vpw(1,1)
        else if(( (s .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j+1) ))then
            call av18pw(lpot,l-2,s,j,t,t1z,t2z,r,vpw)
            av18pot=vpw(2,2)
        else if(( (s .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j-1) ))then
            av18pot=vpw(1,1)
        else if(( (s .eq. 1) .and. (lp .eq. j+1) .and. (l .eq. j-1) ))then
            av18pot=vpw(2,1)
        else if(( (s .eq. 1) .and. (lp .eq. j-1) .and. (l .eq. j+1) ))then
            call av18pw(lpot,l-2,s,j,t,t1z,t2z,r,vpw)
            av18pot=vpw(1,2)
        else
            write(*,*) "error"
            return
        end if
    else
        write(*,*) "error"
        return
    end if

    return
end function

function scs(cutnum,OSTAT,lp,l,kp,k,jj,S,Tz)
    use ICHIRAL,only: CHIRALMOMPWD
    implicit none
    real*8 :: scs
    real*8,intent(in) :: kp,k
    integer,intent(in) :: lp,l,jj,S,Tz,OSTAT,cutnum


    real*8 PPUN1,PPUN2,POTEN(6)
    CHARACTER(len=2) ::FORCE
    if (Tz .eq. -1)then
        FORCE='pp'
    else if (Tz .eq. 0 )then
        FORCE='np'
    else if (Tz .eq. 1 )then
        FORCE='nn'
    else
        write(*,*) 'Wrong Tz!'
    end if
    PPUN1=kp*0.001d0
    PPUN2=k*0.001d0

    CALL CHIRALMOMPWD(OSTAT,FORCE,PPUN1,PPUN2,jj,CUTNUM,POTEN)

    if(jj .eq. 0)then
        if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
            scs=POTEN(1)*1.0D-6
        else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
            scs=POTEN(6)*1.0D-6
        else
            write(*,*) "error1"
            return
        end if
    else if(jj .gt. 0)then
        if( (S .eq. 0) .and. (lp .eq. jj) .and. (l .eq. jj) ) then
            scs=POTEN(1)*1.0D-6
        else if(( (S .eq. 1) .and. (lp .eq. jj) .and. (l .eq. jj) ))then
            scs=POTEN(2)*1.0D-6
        else if(( (S .eq. 1) .and. (lp .eq. jj+1) .and. (l .eq. jj+1) ))then
            scs=POTEN(6)*1.0D-6
        else if(( (S .eq. 1) .and. (lp .eq. jj-1) .and. (l .eq. jj-1) ))then
            scs=POTEN(3)*1.0D-6
        else if(( (S .eq. 1) .and. (lp .eq. jj+1) .and. (l .eq. jj-1) ))then
            scs=-POTEN(5)*1.0D-6
        else if(( (S .eq. 1) .and. (lp .eq. jj-1) .and. (l .eq. jj+1) ))then
            scs=-POTEN(4)*1.0D-6
        else
            write(*,*) "error2"
            return
        end if
    else
        write(*,*) "error3"
        return
    end if

    return
end function

function sms(cutnum,OSTAT,lp,l,kp,k,jj,S,Tz)
    use SMSCHIRAL,only: CHIRALMOMPWD
    implicit none
    real*8 :: sms
    real*8,intent(in) :: kp,k
    integer,intent(in) :: lp,l,jj,S,Tz,OSTAT,cutnum


    real*8 PPUN1,PPUN2,POTEN(6)
    CHARACTER(len=2) ::FORCE
    if (Tz .eq. -1)then
        FORCE='pp'
    else if (Tz .eq. 0 )then
        FORCE='np'
    else if (Tz .eq. 1 )then
        FORCE='nn'
    else
        write(*,*) 'Wrong Tz!'
    end if
    PPUN1=kp*0.001d0
    PPUN2=k*0.001d0
    ! CALL CHIRALMOMPWD(OSTAT,FORCE,ka,kb,jch_BHF,CUTNUM,POTENT)
    CALL CHIRALMOMPWD(OSTAT,FORCE,PPUN1,PPUN2,jj,CUTNUM,POTEN)

    if(jj .eq. 0)then
        if( (S .eq. 0) .and. (lp .eq. 0) .and. (l .eq. 0) ) then
            sms=POTEN(1)*1.0D-6
        else if(( (S .eq. 1) .and. (lp .eq. 1) .and. (l .eq. 1) ))then
            sms=POTEN(6)*1.0D-6
        else
            write(*,*) "wrong channel with j=0"
            return
        end if
    else if(jj .gt. 0)then
        if( (S .eq. 0) .and. (lp .eq. jj) .and. (l .eq. jj) ) then
            sms=POTEN(1)*1.0D-6
        else if(( (S .eq. 1) .and. (lp .eq. jj) .and. (l .eq. jj) ))then
            sms=POTEN(2)*1.0D-6
        else if(( (S .eq. 1) .and. (lp .eq. jj+1) .and. (l .eq. jj+1) ))then
            sms=POTEN(6)*1.0D-6
        else if(( (S .eq. 1) .and. (lp .eq. jj-1) .and. (l .eq. jj-1) ))then
            sms=POTEN(3)*1.0D-6
        else if(( (S .eq. 1) .and. (lp .eq. jj+1) .and. (l .eq. jj-1) ))then
            sms=-POTEN(5)*1.0D-6
        else if(( (S .eq. 1) .and. (lp .eq. jj-1) .and. (l .eq. jj+1) ))then
            sms=-POTEN(4)*1.0D-6
        else
            write(*,*) "wrong channel with j>0"
            return
        end if
    else
        write(*,*) "wrong j"
        return
    end if

    return
end function
