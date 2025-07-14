      module idahopwmodule
      contains

c*********************************************************************
c
c        This code package computes the
c
c        Idaho local position-space two-nucleon potentials from
c        leading to fourth order of chiral effective field theory.
c
c        this package is self-contained and includes
c        all subroutines needed.
c        all codes are consistently in double precision.
c        when working on an UNIX/LINUX system, it is recommended
c        to compile this code with the  "static"  option
c        ( " -s "  on some compilers).
c
c        the program writes out a few comments onto unit=6.
c        if you want to change this, search for `kwrite'
c        and change it from its current setting  kwrite=6.
c*********************************************************************
c
c        October 1, 2022
c
c*********************************************************************
c
c        authors:    S.K. Saha, D.R. Entem, R. Machleidt, and Y. Nosyk
c                     department of physics
c                     university of idaho
c                     moscow, idaho 83844-0903
c                     u. s. a.
c                     e-mail:  machleid@uidaho.edu
c
c        Reference:
c        S. K. Saha, D. R. Entem, R. Machleidt, and Y. Nosyk,
c        "Local position-space two-nucleon potentials from leading to
c        fourth order of chiral effective field theory,"
c        Phys. Rev. C 107, 034002 (2023).
c
c*********************************************************************
c*********************************************************************
c
c        This package includes two major subroutines, namely,
c        subroutine idahopw(lpot,l,s,j,t,t1z,t2z,r,vpw)
c        and
c        subroutine idahoop(lpot,r,vnn)
c
c        idahopw provides the potentials in partial-wave decomposition.
c        idahoop provides the potentials in operator format.
c
c        Besides this, this package includes many auxiliar subroutines
c        (mainly mathematical functions) that the user can ignore.
c
c        The potentials calculated by these codes represent the strong
c        nuclear interaction between two nucleons. Electromagnetic (em)
c        interactions are not included. For the interaction between
c        protons, the user may add the Coulomb potential.
c        However, no em interactions must be added to the np potentials.
c        For explanations, see Sec. III.B of
c        Saha et al., PRC 107, 034002 (2023).
c
c
c *id* idahopw ********************************************************
c        This subroutine calculates the Idaho potentials in
c        partial-wave decomposition.
c
c        arguments:
c        lpot: parameter for potential choice
c            lpot=1   N3LO potential with cutoffs (Rpi,Rct)=(1.0,0.70)fm
c            lpot=2   N3LO potential with cutoffs (Rpi,Rct)=(1.1,0.72)fm
c            lpot=3   N3LO potential with cutoffs (Rpi,Rct)=(1.2,0.75)fm
c            lpot=4   NNLO potential with cutoffs (Rpi,Rct)=(1.0,0.70)fm
c            lpot=5   NNLO potential with cutoffs (Rpi,Rct)=(1.1,0.72)fm
c            lpot=6   NNLO potential with cutoffs (Rpi,Rct)=(1.2,0.75)fm
c            lpot=7   NLO potential with cutoffs (Rpi,Rct)=(1.0,0.70)fm
c            lpot=8   NLO potential with cutoffs (Rpi,Rct)=(1.1,0.72)fm
c            lpot=9   NLO potential with cutoffs (Rpi,Rct)=(1.2,0.75)fm
c            lpot=10  LO potential with cutoffs (Rpi,Rct)=(1.0,0.70)fm
c            lpot=11  LO potential with cutoffs (Rpi,Rct)=(1.1,0.72)fm
c            lpot=12  LO potential with cutoffs (Rpi,Rct)=(1.2,0.75)fm
c
c            lpot=101  v8' version of N3LO potential with
c                                       cutoffs (Rpi,Rct)=(1.0,0.70)fm
c            lpot=102  v8' version of N3LO potential with
c                                       cutoffs (Rpi,Rct)=(1.1,0.72)fm
c            lpot=103  v8' version of N3LO potential with
c                                       cutoffs (Rpi,Rct)=(1.2,0.75)fm
c            (for the definition of v8', see Pudliner et al.,
c             PRC 56, 1720 (1997), Eq. (2.19) therein.
c             the v8' potentials are charge-independent.
c             Note that NNLO and lower-order potentials are
c             v8 potentials from the outset.)
c
c        further arguments:
c        l:    total orbital ang. momentum of the 2 nucleons (0,1,2,...)
c        s:    total spin of the two nucleons (0 or 1)
c        j:    total angular momentum of the two nucleons (0,1,2,...)
c        t:    total isospin of the two nucleons (0 or 1)
c        t1z:  isospin of nucleon 1 (1 for p, -1 for n)
c        t2z:  isospin of nucleon 2 (1 for p, -1 for n)
c        r:    relative distance between the two nucleons in fm
c        vpw:  returned potential in MeV (2x2 array)
c              content of vpw( , ):
c              using, for given j, the symbolic notation v(s;l,l')
c              for the two-nucleon matrix elements, then
c              single channel cases (s=0,1;l=j):
c              vpw(1,1) = v(0;j,j) or v(1;j,j)
c              vpw(2,1) = 0.d0
c              vpw(1,2) = 0.d0
c              vpw(2,2) = 0.d0
c              coupled channel case (s=1;l.ne.j):
c              vpw(1,1) = v(1;j-1,j-1)
c              vpw(2,1) = v(1;j+1,j-1)
c              vpw(1,2) = v(1;j-1,j+1)=vpw(2,1)
c              vpw(2,2) = v(1;j+1,j+1)
c
c
      subroutine idahopw(lpot,l,s,j,t,t1z,t2z,r,vpw)
c
c
      implicit real*8 (a-h,o-z)
c
      common /cwrite/kwrite
      common/copera/opmat(5,7),jj
c
      integer l,s,j,t,t1z,t2z
      dimension vpw(2,2)
c
      dimension vnn(19)
      dimension vveven(5),vvodd(5),vj(5)
      dimension vvcd15(5),vvcd16(5),vvcd17(5)
      dimension vvca18(5),vvca19(5)
      dimension vvcd(5),vvca(5)
      data jold/-1/
      logical index/.true./
      save
c
c
10001 format (//' idahopw: Idaho chiral local position-space',
     1 ' NN potential'/' ',55(1h-))
10002 format (' version used: partial waves, lpot =',i3,', tauz1,2 =',
     1 2i3/' ',55(1h-))
      if (index) then
      kwrite=6
      write (kwrite,10001)
      write (kwrite,10002) lpot,t1z,t2z
      index=.false.
      end if
c
c
      if (j.ne.jold) then
      jold=j
      jj=j
      call opera
      end if
c
c
      do 10 k=1,5
      vveven(k)=0.d0
      vvodd(k)=0.d0
      vvcd15(k)=0.d0
      vvcd16(k)=0.d0
      vvcd17(k)=0.d0
      vvca18(k)=0.d0
      vvca19(k)=0.d0
      vvcd(k)=0.d0
      vvca(k)=0.d0
      vj(k)=0.d0
   10 continue
c
c
      call idahoop(lpot,r,vnn)
c
c
      do 30 iop=1,7
      ieven=2*iop
      iodd=2*iop-1
      do 20 k=1,5
      vveven(k)=vveven(k)+opmat(k,iop)*vnn(ieven)
      vvodd(k)=vvodd(k)+opmat(k,iop)*vnn(iodd)
   20 continue
   30 continue
c
c
      ttt=3.d0*t1z*t2z
      tpt=t1z+t2z
c
c
      do 35 k=1,5
      vvcd15(k)=opmat(k,1)*vnn(15)
      vvcd16(k)=opmat(k,2)*vnn(16)
      vvcd17(k)=opmat(k,3)*vnn(17)
      vvcd(k)=vvcd15(k)+vvcd16(k)+vvcd17(k)
      vvca18(k)=opmat(k,1)*vnn(18)
      vvca19(k)=opmat(k,2)*vnn(19)
      vvca(k)=(vvca18(k)+vvca19(k))*tpt
   35 continue
c
c
c
c
      do 40 k=1,5
         if (mod(jj,2).eq.0) then
c           j is even
            if (k.eq.2) then
            vveven(k)=-3.d0*vveven(k)
            vvcd(k)=vvcd(k)*(ttt+3.d0)
            else
            vvcd(k)=vvcd(k)*(ttt-1.d0)
            end if
         else
c           j is odd
            if (k.ne.2) then
            vveven(k)=-3.d0*vveven(k)
            vvcd(k)=vvcd(k)*(ttt+3.d0)
            else
            vvcd(k)=vvcd(k)*(ttt-1.d0)
            end if
         end if
   40 continue
c
c
      do 50 k=1,5
      vj(k)=vveven(k)+vvodd(k)+vvcd(k)+vvca(k)
   50 continue
c
c
      do 60 i=1,2
      do 60 ii=1,2
      vpw(i,ii)=0.d0
   60 continue
c
      if (s.eq.0) then
      vpw(1,1)=vj(1)
      else
         if (j.eq.0) then
         vpw(1,1)=vj(4)
         else
            if (j.eq.l) then
            vpw(1,1)=vj(2)
            else
            vpw(1,1)=vj(3)
            vpw(2,2)=vj(4)
            vpw(2,1)=vj(5)
            vpw(1,2)=vj(5)
            end if
        end if
      end if
c
c
      return
      end
c
c
c *id* idahoop ********************************************************
c        This subroutine calculates the Idaho potentials in
c        operator format.
c
c        arguments for idahoop:
c        lpot: parameter for potential choice
c            lpot=1   N3LO potential with cutoffs (Rpi,Rct)=(1.0,0.70)fm
c            lpot=2   N3LO potential with cutoffs (Rpi,Rct)=(1.1,0.72)fm
c            lpot=3   N3LO potential with cutoffs (Rpi,Rct)=(1.2,0.75)fm
c            lpot=4   NNLO potential with cutoffs (Rpi,Rct)=(1.0,0.70)fm
c            lpot=5   NNLO potential with cutoffs (Rpi,Rct)=(1.1,0.72)fm
c            lpot=6   NNLO potential with cutoffs (Rpi,Rct)=(1.2,0.75)fm
c            lpot=7   NLO potential with cutoffs (Rpi,Rct)=(1.0,0.70)fm
c            lpot=8   NLO potential with cutoffs (Rpi,Rct)=(1.1,0.72)fm
c            lpot=9   NLO potential with cutoffs (Rpi,Rct)=(1.2,0.75)fm
c            lpot=10  LO potential with cutoffs (Rpi,Rct)=(1.0,0.70)fm
c            lpot=11  LO potential with cutoffs (Rpi,Rct)=(1.1,0.72)fm
c            lpot=12  LO potential with cutoffs (Rpi,Rct)=(1.2,0.75)fm
c
c            lpot=101  v8' version of N3LO potential with
c                                       cutoffs (Rpi,Rct)=(1.0,0.70)fm
c            lpot=102  v8' version of N3LO potential with
c                                       cutoffs (Rpi,Rct)=(1.1,0.72)fm
c            lpot=103  v8' version of N3LO potential with
c                                       cutoffs (Rpi,Rct)=(1.2,0.75)fm
c            (for the definition of v8', see Pudliner et al.,
c             PRC 56, 1720 (1997), Eq. (2.19) therein.
c             the v8' potentials are charge-independent.
c             Note that NNLO and lower-order potentials are
c             v8 potentials from the outset.)
c
c        r:    relative distance between the two nucleons in fm
c
c        vnn: 19 component array containing the potentials (in MeV)
c             attached to the operators with the following
c             order of operators i in vnn(i):
c             i:    1=1                              2=t1.t2
c                   3=s1.s2                          4=(s1.s2)(t1.t2)
c                   5=S12 [=3(s1.r)(s2.r)-s1.s2]     6=S12(t1.t2)
c                   7=L.S                            8=L.S(t1.t2)
c                   9=L**2                          10=L**2(t1.t2)
c                  11=L**2(s1.s2)                  12=L**2(s1.s2)(t1.t2)
c                  13=(L.S)**2                      14=(L.S)**2(t1.t2)
c                  15=T12 [=3*t1z*t2z-t1.t2]        16=(s1.s2)T12
c                  17=S12*T12                       18=t1z+t2z
c                  19=(s1.s2)(t1z+t2z)
c
c                  where s1=sigma_1, t1=tau_1, t1z=tau_1(z), etc.
c
c
      subroutine idahoop(lpot,r,vnn)
c
c
      implicit real*8(a-h,o-z)
c
      common /cwrite/kwrite
c
      common /xcchr/c(20,300),wn,tlamb,ga,fpi,
     1            cb1,cb2,cb3,cb4,cd12,cd3,cd5,cd145,
     2            ic(20,300),mgg(300),ime,
     3            indt(300)
      logical indt
c
      dimension vnn(19)
      dimension vvv(300),vop(30),vst(30),vvop(30)
      dimension vpi(30),vcd(30)
c
      logical indpar/.false./
c
      data hbarc/197.32698d0/
      data pi/3.141592653589793d0/
      data gammae/0.5772156649015329d0/
c
      data xx7/-1.d0/
      data enn8/-1.d0/,xx8/-1.d0/
      data rrcut/-1.d0/,xxl/-1.d0/
      save
c
c
c
c
c         call subroutine chrpar once and only once
c
c
      if (indpar) go to 50
      indpar=.true.
c
      kwrite=6
c
c         headline
c
      write (kwrite,10001)
10001 format (//' idahoop: Idaho chiral local position-space',
     1 ' NN potential'/' ',55(1h-))
      write (kwrite,10002) lpot
10002 format (' version used: operators, lpot =',i3/' ',34(1h-)//)
c
c
      if (lpot.lt.1.or.lpot.gt.103) go to 9001
      if (lpot.gt.12.and.lpot.lt.101) go to 9001
c
c
      if (lpot.gt.100) then
      llpot=lpot-100
      else
      llpot=lpot
      end if
c
c
      call chrpar(llpot)
c
c
      pih=pi*0.5d0
      pi2=pi*pi
      sqrpi=dsqrt(pi)
      ga2 = ga * ga
      ga4 = ga2 * ga2
c
c
   50 xl=r
      xr = xl/hbarc
      xr3=xr*xr*xr
      xr4=xr3*xr
      xr5=xr4*xr
      xr6=xr3*xr3
      xr7=xr6*xr
c
c
      do 65 mes=1,300
      vvv(mes)=0.d0
   65 continue
c
c
c         contributions
c         -------------
c         -------------
c
c
c
c
      do 8500 mes=1,ime
c
      mg=mgg(mes)
c
c
      if (mg.lt.1.or.mg.gt.24) then
      write (kwrite,19000) mg
19000 format(/////' error in idahoop: contribution  ',i4,'   does not
     1 exist in this program.'/' execution terminated.'////)
      stop
      end if
c
c
c
c        create functions and factors
c        ----------------------------
c
c
      c4 = c(4,mes)
      xrc = xr * c4
      xrcd=xrc*2.d0
      xrcq=xrc*4.d0
      xrci=1.d0/xrc
      xrc2=xrc*xrc
      xrc3=xrc2*xrc
      xrc4=xrc2*xrc2
      xrc5=xrc3*xrc2
      xrc6=xrc3*xrc3
      xrc7=xrc4*xrc3
      xrcp1=xrc+1.d0
      xrcp1s=xrcp1*xrcp1
c
c
      if(xrc.gt.200.d0) then
         yukx = 0.d0
         else
         yukx = dexp(-xrc)
         end if
      if(xrcd.gt.200.d0) then
         yukx2 = 0.d0
         else
         yukx2 = dexp(-xrcd)
         end if
c
c
         yukx4 = dexp(xrcq)
c
c
      fac = c(1,mes)
      facpq = c(1,mes)
      mi = 1
      mm = 5
      go to 99
c
c
   95 mi = mi + 3
      mm = mm + 3
c
c
   99 itype = ic(mi,mes)
      if(itype.eq.0) go to 8000
      icase = ic(mi+1,mes)
c
c
 7999 go to(100,200,300,9002,9002,9002,700,800,900,9002,
     1 9002,9002,1300,1400,1500,1600,1700,1800,1900,2000,
     2 9002,9002,9002,9002,2500,2600,9002,2800,2900,3000,
     3 3100,3200,9002,9002,9002,3600,3700,9002,3900,4000,
     4 4100,4200,4300,9002,9002,9002,9002,9002,9002,9002,
     5 5100,5200,5300,5400),itype
c
c
c        for 1PE
c        -------
  100 fac = fac*c(mm,mes)*yukx/xr
      facpq=0.d0
      go to 95
c
c        for 1PE
c        -------
  200 fac = fac*c(mm,mes)*yukx/xr3
     1 *(3.d0+3.d0*xrc+xrc2)
      facpq=0.d0
      go to 95
c
c
c
  300 continue
      nexp = ic(mi+2,mes)
      go to (320,9002),nexp
c
  320 if (c(mm+1,mes).ne.rrcut.or.xl.ne.xxl) then
      rrcut=c(mm+1,mes)
      xxl=xl
      rrcut2=rrcut*rrcut
      rrcut3=rrcut2*rrcut
      rrcut5=rrcut3*rrcut2
      rrcut7=rrcut5*rrcut2
      rrcut9=rrcut7*rrcut2
      rrcut11=rrcut9*rrcut2
      xxl2=xl*xl
      xxl3=xl*xl*xl
      xxl4=xl*xl*xl*xl
      xxlcut2=xxl2/rrcut2
      f01=hbarc/(pi*sqrpi*rrcut2*rrcut)*dexp(-xxlcut2)
      fqq1=f01*(3.d0-2.d0*xxlcut2)*2.d0/rrcut2
c
      if(xxlcut2.gt.200.d0) then
         yukxrs = 0.d0
         else
         yukxrs = dexp(-xxlcut2)
         end if
c
      crs0=yukxrs/rrcut3
      crs1=-2.d0*xxl*yukxrs/rrcut5
      crs2=-2.d0*yukxrs/rrcut5+4.d0*xxl2*yukxrs/rrcut7
      crs3=12.d0*xxl*yukxrs/rrcut7-8.d0*xxl3*yukxrs/rrcut9
      crs4=12.d0*yukxrs/rrcut7-48.d0*xxl2*yukxrs/rrcut9+16.d0
     1      *xxl4*yukxrs/rrcut11
      end if
c
c
      go to (321,322,9002,9002,325,326,327,9002,9002,9002,
     1        331,332,333,334,335,336,337,338),icase
c
  321 fac = fac*f01
      facpq=0.d0
      go to 95
c
  322 fac = fac*fqq1
      facpq=0.d0
      go to 95
c
  325 fac = fac*f01*2.d0/rrcut2
      facpq=0.d0
      go to 95
c
  326 fac = fac*fqq1/3.d0
      facpq=0.d0
      go to 95
c
  327 fac = -fac*f01*4.d0*xxlcut2/(3.d0*rrcut2)
      facpq=0.d0
      go to 95

  331 fac = fac*c(mm,mes)*crs0
      facpq=0.d0
      go to 95
c
  332 fac = fac*c(mm,mes)*(-crs2-(2.d0/xxl)*crs1)
      facpq=0.d0
      go to 95
c
  333 fac = fac*c(mm,mes)*(crs4+(4.d0/xxl)*crs3)
      facpq=0.d0
      go to 95
c
  334 fac = fac*c(mm,mes)*(-crs2+(1.d0/xxl)*crs1)
      facpq=0.d0
      go to 95
c
  335 fac = fac*c(mm,mes)*(crs4+(1.d0/xxl)*crs3
     1      -(6.d0/xxl2)*crs2+(6.d0/xxl3)*crs1)
      facpq=0.d0
      go to 95
c
  336 fac =-fac*c(mm,mes)*(1.d0/xxl)*crs1
      facpq=0.d0
      go to 95
c
  337 fac = fac*c(mm,mes)*((1.d0/xxl)*crs3
     1      +(2.d0/xxl2)*crs2-(2.d0/xxl3)*crs1)
      facpq=0.d0
      go to 95
c
  338 fac = fac*c(mm,mes)*(1.d0/xxl2)*(-crs2+(1.d0/xxl)
     1      *crs1)
      facpq=0.d0
      go to 95
c
c
  700 enn=c(mm,mes)
      rcut=c(mm+1,mes)
      xx=(xl/rcut)**(2.d0*enn)
      if (xx.ne.xx7) then
      xx7=xx
         if(xx.ge.0.1d0) then
            cut7=1.d0-dexp(-xx)
            else
            cut7=oneme(xx)
         endif
      endif
      fac=fac*cut7
      facpq=0.d0
      go to 95
c
  800 enn=c(mm,mes)
      rcut=c(mm+1,mes)
      xx=(xl/rcut)**2.d0
      if (enn.ne.enn8.or.xx.ne.xx8) then
      enn8=enn
      xx8=xx
         if(xx.ge.0.1d0) then
            cut8=(1.d0-dexp(-xx))**enn
            else
            cut8=(oneme(xx))**enn
         endif
      endif
      fac=fac*cut8
      facpq=0.d0
      go to 95
c
c        Piarulli's pi ff
  900 enn=c(mm,mes)
      rcut=c(mm+1,mes)
      xx1=xl/rcut
      xxn=xx1**(2.d0*enn)
      cut9=1.d0-1.d0/(xxn*dexp(2.d0*(xx1-1.d0))+1.d0)
      fac=fac*cut9
      facpq=0.d0
      go to 95
c
c
 1300 fac=fac*c(mm,mes)/xr4
     1   *((1.d0+2.d0*ga2*(5.d0+2.d0*xrc2)
     2   -ga4*(23.d0+12.d0*xrc2))*dbesk1(xrcd)
     3   +xrc*(1.d0+10.d0*ga2-ga4*(23.d0+4.d0*xrc2))
     4   *dbesk0(xrcd))
      facpq=0.d0
      go to 95
c
 1400 go to (1410,1420),icase
 1410 fac=fac*c(mm,mes)/xr4
     1   *(3.d0*xrc*dbesk0(xrcd)
     2   +(3.d0+2.d0*xrc2)*dbesk1(xrcd))
      facpq=0.d0
      go to 95
c
 1420 fac=fac*c(mm,mes)/xr4
     1   *(-12.d0*xrc*dbesk0(xrcd)
     2   -(15.d0+4.d0*xrc2)*dbesk1(xrcd))
      facpq=0.d0
      go to 95
c
 1500 fac=fac*c(mm,mes)*yukx2/xr6
     1   *((2.d0*cb1+3.d0*ga2/(16.d0*wn))*xrc2*xrcp1s
     2   +ga2*xrc5/(32.d0*wn)
     3   +(cb3+3.d0*ga2/(16.d0*wn))
     4   *(6.d0+12.d0*xrc+10.d0*xrc2+4.d0*xrc3+xrc4))
      facpq=0.d0
      go to 95
c
 1600 fac=fac*c(mm,mes)*yukx2/xr5
     1   *(2.d0*(3.d0*ga2-2.d0)
     2   *(6.d0*xrci+12.d0+10.d0*xrc+4.d0*xrc2+xrc3)
     3   +ga2*xrc*(2.d0+4.d0*xrc+2.d0*xrc2+3.d0*xrc3))
      facpq=0.d0
      go to 95
c
 1700 go to (1710,1720),icase
 1710 fac=fac*c(mm,mes)*yukx2/xr5
     1   *(6.d0*xrci+12.d0+11.d0*xrc+6.d0*xrc2+2.d0*xrc3)
      facpq=0.d0
      go to 95
c
 1720 fac=fac*c(mm,mes)*yukx2/xr5
     2   *(12.d0*xrci+24.d0+20.d0*xrc+9.d0*xrc2+2.d0*xrc3)
      facpq=0.d0
      go to 95
c
 1800 go to (1810,1820),icase
 1810 fac=fac*c(mm,mes)*yukx2/xr6
     1   *((cb4+1.d0/(4.d0*wn))
     2   *xrcp1*(3.d0+3.d0*xrc+2.d0*xrc2)
     3   -ga2/(16.d0*wn)
     4   *(18.d0+36.d0*xrc+31.d0*xrc2+14.d0*xrc3+2.d0*xrc4))
      facpq=0.d0
      go to 95
c
 1820 fac=fac*c(mm,mes)*yukx2/xr6
     1   *(-(cb4+1.d0/(4.d0*wn))*xrcp1*(3.d0+3.d0*xrc+xrc2)
     2   +ga2/(32.d0*wn)
     3   *(36.d0+72.d0*xrc+52.d0*xrc2+17.d0*xrc3+2.d0*xrc4))
      facpq=0.d0
      go to 95
c
 1900 fac=fac*c(mm,mes)*yukx2/xr6
     1   *xrcp1
     2   *(2.d0+2.d0*xrc+xrc2)
      facpq=0.d0
      go to 95
c
 2000 fac=fac*c(mm,mes)*yukx2/xr6
     1   *xrcp1s
      facpq=0.d0
      go to 95
c
 2500 fac=fac*c(mm,mes)*yukx2/xr6
     1   *((2.d0*cb1)*xrc2*xrcp1s
     2   +(cb3)
     3   *(6.d0+12.d0*xrc+10.d0*xrc2+4.d0*xrc3+xrc4))
      facpq=0.d0
      go to 95
c
 2600 fac=fac*c(mm,mes)*yukx2/xr6
     1   *((3.d0*ga2/(16.d0*wn))*xrc2*xrcp1s
     2   +ga2*xrc5/(32.d0*wn)
     3   +(3.d0*ga2/(16.d0*wn))
     4   *(6.d0+12.d0*xrc+10.d0*xrc2+4.d0*xrc3+xrc4))
      facpq=0.d0
      go to 95
c
 2800 go to (2810,2820),icase
 2810 fac=fac*c(mm,mes)*yukx2/xr6
     1   *((cb4)
     2   *xrcp1*(3.d0+3.d0*xrc+2.d0*xrc2))
      facpq=0.d0
      go to 95
c
 2820 fac=fac*c(mm,mes)*yukx2/xr6
     1   *(-(cb4)*xrcp1*(3.d0+3.d0*xrc+xrc2))
      facpq=0.d0
      go to 95
c
 2900 go to (2910,2920),icase
 2910 fac=fac*c(mm,mes)*yukx2/xr6
     1   *((1.d0/(4.d0*wn))
     2   *xrcp1*(3.d0+3.d0*xrc+2.d0*xrc2)
     3   -ga2/(16.d0*wn)
     4   *(18.d0+36.d0*xrc+31.d0*xrc2+14.d0*xrc3+2.d0*xrc4))
      facpq=0.d0
      go to 95
c
 2920 fac=fac*c(mm,mes)*yukx2/xr6
     1   *(-(1.d0/(4.d0*wn))*xrcp1*(3.d0+3.d0*xrc+xrc2)
     2   +ga2/(32.d0*wn)
     3   *(36.d0+72.d0*xrc+52.d0*xrc2+17.d0*xrc3+2.d0*xrc4))
      facpq=0.d0
      go to 95
c
 3000 fac=fac*c(mm,mes)*(1.d0/xrc5)*((3.d0*cb2*cb2+20.d0*cb2*cb3
     1    +60.d0*cb3*cb3+4.d0*(2.d0*cb1+cb3)**2*xrc2)*xrc*dbesk1(xrcd)
     2	  +2.d0*(3.d0*cb2*cb2+20.d0*cb2*cb3+60.d0*cb3*cb3
     3	  +2.d0*(2.d0*cb1+cb3)*(cb2+6.d0*cb3)*xrc2)*dbesk2(xrcd))
      facpq=0.d0
      go to 95
c
 3100 go to (3110,3120),icase
 3110 fac=fac*c(mm,mes)*(1.d0/xrc4)*(2.d0*xrc*dbesk2(xrcd)
     1    +5.d0*dbesk3(xrcd))
      facpq=0.d0
      go to 95
c
 3120 fac=fac*c(mm,mes)*(1.d0/xrc5)*((3.d0+4.d0*xrc2)*dbesk2(xrcd)
     1    +16.d0*xrc*dbesk3(xrcd))
      facpq=0.d0
      go to 95
c
 3200 fac=fac*c(mm,mes)*(1.d0/xrc5)*(dbesk2(xrcd)
     1    +2.d0*xrc*dbesk3(xrcd))
      facpq=0.d0
      go to 95
c
 3600 go to (3610,3620),icase
 3610 fac=fac*c(mm,mes)*(1.d0/xrc5)*((5.d0+25.d0*ga2+4.d0*ga2*xrc2)
     1    *xrc*dbesk1(xrcd)+2.d0*(5.d0+25.d0*ga2+(1.d0+8.d0*ga2)
     2    *xrc2)*dbesk2(xrcd))
      facpq=0.d0
      go to 95
c
 3620 fac=fac*c(mm,mes)*(1.d0/xrc5)*((1.d0+6.d0*ga2)*2.d0*xrc
     1    *dbesk1(xrcd)+(5.d0+25.d0*ga2+4.d0*ga2*xrc2)*dbesk2(xrcd))
      facpq=0.d0
      go to 95
c
 3700 go to (3710,3720),icase
 3710 fac=fac*c(mm,mes)*(1.d0/xrc5)*((5.d0-35.d0*ga2-4.d0*ga2*xrc2)*xrc
     1     *dbesk1(xrcd)+2.d0*(5.d0*(1.d0-7.d0*ga2)+(1.d0-10.d0*ga2)
     2     *xrc2)*dbesk2(xrcd))
      facpq=0.d0
      go to 95
c
 3720 fac=fac*c(mm,mes)*(1.d0/xrc5)*(2.d0*(-8.d0+59.d0*ga2+4.d0
     1    *ga2*xrc2)*xrc*dbesk1(xrcd)-(35.d0*(1.d0-7.d0*ga2)+4.d0*
     2       (1.d0-13.d0*ga2)*xrc2)*dbesk2(xrcd))
      facpq=0.d0
      go to 95
c
 3900 fac=fac*c(mm,mes)*(1.d0/xrc6)*((20.d0*(cb2-6.d0*cb3)-4.d0
     1    *(6.d0*cb1-cb2+9.d0*cb3)*xrc2-2.d0*(2.d0*cb1+cb3)*xrc4)
     2    *xrc*dbesk0(xrcd)+(20.d0*(cb2-6.d0*cb3)-2.d0*(12.d0*cb1
     3    -7.d0*cb2+48.d0*cb3)*xrc2-(16.d0*cb1-cb2+10.d0*cb3)*xrc4)
     4    *dbesk1(xrcd))
      facpq=0.d0
      go to 95
c
 4000 go to (4010,4020),icase
 4010 fac=fac*c(mm,mes)*yukx2/xrc6*(24.d0+48.d0*xrc
     1   +43.d0*xrc2+22.d0*xrc3+7.d0*xrc4+4.d0*ga2*
     2   (6.d0+12.d0*xrc+10.d0*xrc2+4.d0*xrc3+xrc4))
      facpq=0.d0
      go to 95
c
 4020 fac= fac*c(mm,mes)*yukx2/xrc7*((120.d0+ 240.d0*xrc+213.d0*xrc2
     1	   +106.d0*xrc3+ 32.d0*xrc4+8.d0*xrc5)*
     2     (dlog(4.d0*xrc)+gammae)-(120.d0- 240.d0*xrc
     3	   + 213.d0*xrc2- 106.d0*xrc3+ 32.d0*xrc4-8.d0*xrc5)*yukx4
     4	    *eix(-4.d0*xrc)-4.d0*xrc*(96.d0+72.d0*xrc+38.d0*xrc2
     5       +7.d0*xrc3))-2.d0*fac*c(mm,mes)*(aibarm1(xrcd))/xrc
      facpq=0.d0
      go to 95
c
 4100 go to (4110,4120,4130),icase
 4110 fac=fac*c(mm,mes)*(yukx2/xrc7)*((15.d0+30.d0*xrc+24.d0*xrc2
     1    +8.d0*xrc3)*(dlog(xrcq)+gammae)+(-15.d0+30.d0*xrc-24.d0*xrc2
     2    +8.d0*xrc3)*yukx4*eix(-xrcq)-4.d0*xrc*(15.d0+15.d0*xrc
     3    +8.d0*xrc2+2.d0*xrc3)-8.d0*ga2*xrc*(3.d0+6.d0*xrc+5.d0*xrc2
     4    +2.d0*xrc3))
      facpq=0.d0
      go to 95
c
 4120 fac=fac*c(mm,mes)*(yukx2/xrc6)*(3.d0+6.d0*xrc+4.d0*xrc2+xrc3)
      facpq=0.d0
      go to 95
c
 4130 fac=fac*c(mm,mes)*(yukx2/xrc7)*(-324.d0*xrc-228.d0*xrc2-48.d0
     1    *xrc3+5.d0*(21.d0+42.d0*xrc+30.d0*xrc2+4.d0*xrc3)*(dlog(xrcq)
     2     +gammae)+5.d0*yukx4*eix(-xrcq)*(-21.d0+42.d0*xrc-30.d0*xrc2
     3    +4.d0*xrc3))+24.d0*fac*c(mm,mes)*(1.d0/xrc3)
     4    *aibarm1(xrcd)
      facpq=0.d0
      go to 95
c
 4200 go to (4210,4220),icase
 4210 fac=fac*c(mm,mes)*(1.d0/xrc4)*(2.d0*xrc*dbesk2(xrcd)+5.d0
     1       *dbesk3(xrcd))
      facpq=0.d0
      go to 95
c
 4220 fac=fac*c(mm,mes)*(1.d0/xrc5)*((3.d0+4.d0*xrc2)*dbesk2(xrcd)
     1     +16.d0*xrc*dbesk3(xrcd))
      facpq=0.d0
      go to 95
c
 4300 go to (4310,4320,4330,4340),icase
 4310 fac=fac*c(mm,mes)*(1.d0/xrc7)*((30.d0+89.d0*xrc2-8.d0*xrc4
     1    +ga2*(300.d0+926.d0*xrc2-32.d0*xrc4)+ga4*(750.d0+2405.d0
     2	  *xrc2+76.d0*xrc4))*dbesk0(xrcd)+(137.d0+8.d0*xrc2+8.d0*xrc4+
     3       2.d0*ga2*(685.d0+106.d0*xrc2+16.d0*xrc4)+ga4*(3425.d0
     4      +860.d0*xrc2+32.d0*xrc4))*xrc*dbesk1(xrcd))-16.d0*
     5      fac*c(mm,mes)*(1.d0+2.d0*ga2)**2*
     6      (aitilm1(xrcd)/xrc)
      facpq=0.d0
      go to 95
c
 4320 fac=fac*c(mm,mes)*(-((2.d0*ga2*xrc*dbesk1(xrcd)+(1.d0+5.d0*ga2)
     1    *dbesk2(xrcd))/xrc3)*2.d0*cd5+(((5.d0+ga2*(25.d0+2.d0*xrc2))
     2    *xrc*dbesk1(xrcd)+(10.d0+xrc2+ga2*(50.d0+11.d0*xrc2))
     3    *dbesk2(xrcd))/xrc5)*cd12)
      facpq=0.d0
      go to 95
c
 4330 fac=fac*c(mm,mes)*(1.d0/xrc5)*((25.d0+ga2*(190.d0-4.d0*xrc2)
     1    +ga4*(325.d0+4.d0
     2    *xrc2))*xrc*dbesk1(xrcd)+2.d0*(25.d0-xrc2+ga2*(190.d0+11.d0
     3    *xrc2)+ga4*(325.d0+44.d0*xrc2))*dbesk2(xrcd))
      facpq=0.d0
      go to 95
c
 4340 fac=fac*c(mm,mes)*((2.d0*ga2*xrc*dbesk2(xrcd)+(3.d0+7.d0*ga2)
     1    *dbesk3(xrcd))*cd3/xrc4)
      facpq=0.d0
      go to 95
c
 5100 fac=fac*c(mm,mes)*(yukx2/xr6)*(24.d0+48.d0*xrc+46.d0*xrc2
     1    +28.d0*xrc3+10.d0*xrc4+xrc5)
      facpq=0.d0
      go to 95
c
 5200 fac=fac*c(mm,mes)*(yukx2/xr6)*(ga2*(48.d0+96.d0*xrc+82.d0*xrc2
     1    +36.d0*xrc3+10.d0*xrc4+3.d0*xrc5)-4.d0*(6.d0+12.d0*xrc
     2    +10.d0*xrc2+4.d0*xrc3+xrc4))
      facpq=0.d0
      go to 95
c
 5300 go to (5310,5320),icase
 5310 fac=fac*c(mm,mes)*(yukx2/xr6)*(24.d0+48.d0*xrc+43.d0*xrc2
     1    +22.d0*xrc3+6.d0*xrc4)
      facpq=0.d0
      go to 95
c
 5320 fac=fac*c(mm,mes)*(yukx2/xr6)*(96.d0+192.d0*xrc+152.d0*xrc2
     1    +62.d0*xrc3+12.d0*xrc4)
      facpq=0.d0
      go to 95
c
 5400 go to (5410,5420),icase
 5410 fac=fac*c(mm,mes)*(yukx2/xr6)*(2.d0*ga2*(2.d0+xrc)*(6.d0
     1    +9.d0*xrc+6.d0*xrc2+2.d0*xrc3)-8.d0*(1.d0+xrc)*(3.d0
     2    +3.d0*xrc+2.d0*xrc2))
      facpq=0.d0
      go to 95
c
 5420 fac=fac*c(mm,mes)*(yukx2/xr6)*(2.d0*ga2*(2.d0+xrc)*(12.d0
     1    +18.d0*xrc+9.d0*xrc2+2.d0*xrc3)-16.d0*(1.d0+xrc)*(3.d0
     2    +3.d0*xrc+xrc2))
      facpq=0.d0
      go to 95
c
c
 8000 continue
c
c
      vvv(mes)=fac
c
c
 8500 continue
c
c
      do 8605 iop=1,30
      vpi(iop)=0.d0
      vop(iop)=0.d0
      vvop(iop)=0.d0
      vst(iop)=0.d0
      vcd(iop)=0.d0
 8605 continue
c
c
c         cream-off ope
c
      do 8610 mes=1,4
      vpi(mes)=vvv(mes)
 8610 continue
c
c
      do 8615 mes=5,ime
c
      if (mgg(mes).le.7) then
c
      if (mgg(mes).eq.1.and..not.indt(mes)) then
      vop(1)=vop(1)+vvv(mes)
      go to 8615
      end if
      if (mgg(mes).eq.1.and.indt(mes)) then
      vop(2)=vop(2)+vvv(mes)
      go to 8615
      end if
      if (mgg(mes).eq.2.and..not.indt(mes)) then
      vop(3)=vop(3)+vvv(mes)
      go to 8615
      end if
      if (mgg(mes).eq.2.and.indt(mes)) then
      vop(4)=vop(4)+vvv(mes)
      go to 8615
      end if
      if (mgg(mes).eq.3.and..not.indt(mes)) then
      vop(5)=vop(5)+vvv(mes)
      go to 8615
      end if
      if (mgg(mes).eq.3.and.indt(mes)) then
      vop(6)=vop(6)+vvv(mes)
      go to 8615
      end if
      if (mgg(mes).eq.4.and..not.indt(mes)) then
      vop(7)=vop(7)+vvv(mes)
      go to 8615
      end if
      if (mgg(mes).eq.4.and.indt(mes)) then
      vop(8)=vop(8)+vvv(mes)
      go to 8615
      end if
      if (mgg(mes).eq.5.and..not.indt(mes)) then
      vop(9)=vop(9)+vvv(mes)
      go to 8615
      end if
      if (mgg(mes).eq.5.and.indt(mes)) then
      vop(10)=vop(10)+vvv(mes)
      go to 8615
      end if
      if (mgg(mes).eq.6.and..not.indt(mes)) then
      vop(11)=vop(11)+vvv(mes)
      go to 8615
      end if
      if (mgg(mes).eq.6.and.indt(mes)) then
      vop(12)=vop(12)+vvv(mes)
      go to 8615
      end if
      if (mgg(mes).eq.7.and..not.indt(mes)) then
      vop(13)=vop(13)+vvv(mes)
      go to 8615
      end if
      if (mgg(mes).eq.7.and.indt(mes)) then
      vop(14)=vop(14)+vvv(mes)
      go to 8615
      end if
c
      else
            if (mgg(mes).le.21) then
c              ST contributions
            vst(mgg(mes)-7)=vst(mgg(mes)-7)+vvv(mes)
            else
            vcd(mgg(mes)-21)=vvv(mes)
            end if
      end if
c
 8615 continue
c
c
c        the charge-independent LO c01 term
c
      vcd(4)=(vcd(1)+vcd(2)+vcd(3))/3.d0
      vst(2)=vst(2)+vcd(4)
c
c
      vvop(1)=0.0625d0*( vst(1)+3.d0*vst(2)+3.d0*vst(3)+9.d0*vst(4))
      vvop(2)=0.0625d0*(-vst(1)+1.d0*vst(2)-3.d0*vst(3)+3.d0*vst(4))
      vvop(3)=0.0625d0*(-vst(1)-3.d0*vst(2)+1.d0*vst(3)+3.d0*vst(4))
      vvop(4)=0.0625d0*( vst(1)-1.d0*vst(2)-1.d0*vst(3)+1.d0*vst(4))
c
      vvop(5)=0.25d0*( vst(5)+3.d0*vst(6))
      vvop(6)=0.25d0*(-vst(5)+1.d0*vst(6))
c
      vvop(7)=0.25d0*( vst(7)+3.d0*vst(8))
      vvop(8)=0.25d0*(-vst(7)+1.d0*vst(8))
c
      vvop(9)=0.0625d0*( vst(9)+3.d0*vst(10)+3.d0*vst(11)+9.d0*vst(12))
      vvop(10)=.0625d0*(-vst(9)+1.d0*vst(10)-3.d0*vst(11)+3.d0*vst(12))
      vvop(11)=.0625d0*(-vst(9)-3.d0*vst(10)+1.d0*vst(11)+3.d0*vst(12))
      vvop(12)=.0625d0*( vst(9)-1.d0*vst(10)-1.d0*vst(11)+1.d0*vst(12))
c
      vvop(13)=0.25d0*( vst(13)+3.d0*vst(14))
      vvop(14)=0.25d0*(-vst(13)+1.d0*vst(14))
c
c
      do 8625 iop=1,19
      vnn(iop)=0.0d0
 8625 continue
c
c
      do 8635 iop=1,14
      vnn(iop)=vop(iop)+vvop(iop)
 8635 continue
c
c
c        charge dependence from contacts
      vnn(15)=(0.5*(vcd(1)+vcd(3))-vcd(2))/6.d0/4.d0
      vnn(16)=-vnn(15)
c
c        charge asymmetry from contacts
      vnn(18)=(vcd(1)-vcd(3))/4.d0/4.d0
      vnn(19)=-vnn(18)
c
c
c        add charge-dependence due to ope
c
      vnn(16)=(vpi(1)-vpi(3))/3.d0+vnn(16)
      vnn(17)=(vpi(2)-vpi(4))/3.d0
c
c
c        the charge-independent ope
c
      vpi(1)=(vpi(1)+2.d0*vpi(3))/3.d0
      vpi(2)=(vpi(2)+2.d0*vpi(4))/3.d0
      vnn(4)=vnn(4)+vpi(1)
      vnn(6)=vnn(6)+vpi(2)
c
c
      if (lpot.lt.100) then
      return
      else
c
c     v8' projection
c     --------------
c
      vnn(1)=vnn(1)+1.25d0*vnn( 9)+0.75d0*vnn(10)+0.75d0*vnn(11)
     x             +2.25d0*vnn(12)+0.75d0*vnn(13)+0.75d0*vnn(14)
      vnn(2)=vnn(2)+0.25d0*vnn( 9)+0.75d0*vnn(10)+0.75d0*vnn(11)
     x             -0.75d0*vnn(12)+0.25d0*vnn(13)+0.25d0*vnn(14)
      vnn(3)=vnn(3)+0.25d0*vnn( 9)+0.75d0*vnn(10)+0.75d0*vnn(11)
     x             -0.75d0*vnn(12)+0.25d0*vnn(13)+0.25d0*vnn(14)
      vnn(4)=vnn(4)+0.25d0*vnn( 9)-0.25d0*vnn(10)-0.25d0*vnn(11)
     x             +1.25d0*vnn(12)+1.d0/12.d0*vnn(13)
     x             +1.d0/12.d0*vnn(14)
      vnn(5)=vnn(5)-0.3125d0*vnn(13)-0.3125d0*vnn(14)
      vnn(6)=vnn(6)-5.d0/48.d0*vnn(13)-5.d0/48.d0*vnn(14)
      vnn(7)=vnn(7)-0.50d0*vnn( 9)+1.50d0*vnn(10)-0.50d0*vnn(11)
     x             +1.50d0*vnn(12)-1.125d0*vnn(13)+1.875d0*vnn(14)
      vnn(8)=vnn(8)+0.50d0*vnn( 9)-1.50d0*vnn(10)+0.50d0*vnn(11)
     x             -1.50d0*vnn(12)+0.625d0*vnn(13)-2.375d0*vnn(14)
c
      do 8645 iop=9,19
      vnn(iop)=0.0d0
 8645 continue
c
      end if
c
      return
c
c
c         more errors and warnings
c         ------------------------
c
c
 9001 write (kwrite,19001) lpot
19001 format (///' error in idahoop: lpot =',i5,
     1'   does not exist in this program.'/' execution terminated.'
     2////)
      go to 9999
c
 9002 write (kwrite,19002) itype,icase
19002 format (///' error in idahoop: cut/fun type',i5,' and icase',
     1i5,' does not exist in this program.'/' execution terminated.'
     2////)
c
 9999 stop
      end
c
c
      subroutine chrpar(lpot)
c
c
c        chrpar reads, writes and stores the parameters for idahoop
c
c
      implicit real*8(a-h,o-z)
c
      common /cwrite/kwrite
c
      common /xcchr/c(20,300),wn,tlamb,ga,fpi,
     1            cb1,cb2,cb3,cb4,cd12,cd3,cd5,cd145,
     2            ic(20,300),mgg(300),ime,
     3            indt(300)
c
      logical indt
      dimension cc(4),cca(4)
      dimension cst(30),cvw(30)
      dimension tab(4,136),ttab(4,136)
      dimension tabct(12,4,58)
      character*4 name(3)
      character*4 opsym(24)/'c   ','ss  ','s12 ',
     1'ls  ','l2  ','l2ss','ls2 ',
     2'c00 ','c01 ','c10 ','c11 ',
     3's120','s121','ls0 ','ls1 ',
     4'l200','l201','l210','l211',
     5'ls20','ls21','c01p','c01n','c01m'/
      character*4 end/'end '/
      character*4 cut/'cut '/,cuta/'cuta'/,fun/'fun '/
      character*4 ntab(3,136)
      logical indca
      logical index/.false./
c
      data hbarc/197.32698d0/
      data pi/3.141592653589793d0/
      save
c
c
c
c
c         parameter tables
c         ----------------
c         ----------------
c
c
c         identification tables
c         ---------------------
c
c
      data ((ntab(i,j),i=1,3),j=1,78)/
     1 'cuta','ll  ','  ',
     2 'ss  ',' ope','p ',
     3 'fun ','    ','  ',
     4 's12 ',' ope','p ',
     5 'fun ','    ','  ',
     6 'ss  ',' ope','p ',
     7 'fun ','    ','  ',
     8 's12 ',' ope','p ',
     9 'fun ','    ','  ',
     * 'cuta','ll  ','  ',
     1 'c   ',' tpn','1 ',
     2 'fun ','    ','  ',
     3 'ss  ',' tpn','1 ',
     4 'fun ','    ','  ',
     5 's12 ',' tpn','1 ',
     6 'fun ','    ','  ',
     7 'c   ',' tpn','2 ',
     8 'fun ','    ','  ',
     9 'ss  ',' tpn','2 ',
     * 'fun ','    ','  ',
     1 's12 ',' tpn','2 ',
     2 'fun ','    ','  ',
     3 'c   ',' tpn','2m',
     4 'fun ','    ','  ',
     5 'c   ',' tpn','2m',
     6 'fun ','    ','  ',
     7 'ss  ',' tpn','2m',
     8 'fun ','    ','  ',
     9 's12 ',' tpn','2m',
     * 'fun ','    ','  ',
     1 'ss  ',' tpn','2m',
     2 'fun ','    ','  ',
     3 's12 ',' tpn','2m',
     4 'fun ','    ','  ',
     5 'ls  ',' tpn','2m',
     6 'fun ','    ','  ',
     7 'ls  ',' tpn','2m',
     8 'fun ','    ','  ',
     9 'c   ',' tpn','3 ',
     * 'fun ','    ','  ',
     1 'ss  ',' tpn','3 ',
     2 'fun ','    ','  ',
     3 's12 ',' tpn','3 ',
     4 'fun ','    ','  ',
     5 'c   ',' tpc','m ',
     6 'fun ','    ','  ',
     7 'c   ',' tpc','m ',
     8 'fun ','    ','  ',
     9 'ss  ',' tpc','m ',
     * 'fun ','    ','  ',
     1 's12 ',' tpc','m ',
     2 'fun ','    ','  ',
     3 'ls  ',' tpc','m ',
     4 'fun ','    ','  ',
     5 'ls  ',' tpc','m ',
     6 'fun ','    ','  ',
     7 'c   ',' tpn','32',
     8 'fun ','    ','  ',
     9 'c   ',' tpn','32',
     * 'fun ','    ','  ',
     1 'c   ',' tpn','32',
     2 'fun ','    ','  ',
     3 'c   ',' tpn','32',
     4 'fun ','    ','  ',
     5 'c   ',' tpn','32',
     6 'fun ','    ','  ',
     7 'c   ',' tpn','32',
     8 'fun ','    ','  ',
     9 'ss  ',' tpn','32',
     * 'fun ','    ','  ',
     1 's12 ',' tpn','32',
     2 'fun ','    ','  ',
     3 'ss  ',' tpn','32',
     4 'fun ','    ','  ',
     5 's12 ',' tpn','32',
     6 'fun ','    ','  ',
     7 's12 ',' tpn','32',
     8 'fun ','    ','  '/
c
      data ((ntab(i,j),i=1,3),j=79,136)/
     9 'cuta','ll  ','  ',
     * 'c01n',' 1s0','  ',
     1 'fun ','    ','  ',
     2 'c01 ',' 1s0','  ',
     3 'fun ','    ','  ',
     4 'c01 ',' 1s0','  ',
     5 'fun ','    ','  ',
     6 'l201','  1s','0 ',
     7 'fun ','    ','  ',
     8 'c00 ',' 1p1','  ',
     9 'fun ','    ','  ',
     * 'c00 ',' 1p1','  ',
     1 'fun ','    ','  ',
     2 'c00 ',' 1p1','  ',
     3 'fun ','    ','  ',
     4 'l200','  1p','1 ',
     5 'fun ','    ','  ',
     6 'c10 ',' 3s1','  ',
     7 'fun ','    ','  ',
     8 'c10 ',' 3s1','  ',
     9 'fun ','    ','  ',
     * 'c10 ',' 3s1','  ',
     1 'fun ','    ','  ',
     2 's120',' 3s1','  ',
     3 'fun ','    ','  ',
     4 's120',' 3s1','  ',
     5 'fun ','    ','  ',
     6 'ls0 ',' 3s1','  ',
     7 'fun ','    ','  ',
     8 'ls0 ',' 3s1','  ',
     9 'fun ','    ','  ',
     * 'ls20','  3s','1 ',
     1 'fun ','    ','  ',
     2 'l210','  3s','1 ',
     3 'fun ','    ','  ',
     4 'c11 ',' 3pj','  ',
     5 'fun ','    ','  ',
     6 'c11 ',' 3pj','  ',
     7 'fun ','    ','  ',
     8 'c11 ',' 3pj','  ',
     9 'fun ','    ','  ',
     * 's121',' 3pj','  ',
     1 'fun ','    ','  ',
     2 's121',' 3pj','  ',
     3 'fun ','    ','  ',
     4 'ls1 ',' 3pj','  ',
     5 'fun ','    ','  ',
     6 'ls1 ',' 3pj','  ',
     7 'fun ','    ','  ',
     8 'ls21','  3p','j ',
     9 'fun ','    ','  ',
     * 'l211','  3p','j ',
     1 'fun ','    ','  ',
     2 'c01p',' 1s0','  ',
     3 'fun ','    ','  ',
     4 'c01m',' 1s0','  ',
     5 'fun ','    ','  ',
     6 'end ','para','m.'/
c
c
c         numerial parameters
c         -------------------
c
c
      data ((tab(i,j),i=1,4),j=1,78)/
     1     7.00000000d0,   0.00d0,    5.0000d0,   1.00d0,
     2     1.29000000d0,  92.40d0,  134.9766d0,   0.00d0,
     3     1.00000000d0,   0.00d0,    0.0000d0,   0.00d0,
     4     1.29000000d0,  92.40d0,  134.9766d0,   0.00d0,
     5     2.00000000d0,   0.00d0,    0.0000d0,   0.00d0,
     6     1.29000000d0,  92.40d0,  139.5702d0,   0.00d0,
     7     1.00000000d0,   0.00d0,    0.0000d0,   0.00d0,
     8     1.29000000d0,  92.40d0,  139.5702d0,   0.00d0,
     9     2.00000000d0,   0.00d0,    0.0000d0,   0.00d0,
     *     8.00000000d0,   0.00d0,    5.0000d0,   1.00d0,
     1     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     2    13.00000000d0,   0.00d0,    0.0000d0,   0.00d0,
     3     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4    14.00000000d0,   1.00d0,    0.0000d0,   0.00d0,
     5     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6    14.00000000d0,   2.00d0,    0.0000d0,   0.00d0,
     7     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8    25.00000000d0,   0.00d0,    0.0000d0,   0.00d0,
     9     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     *    28.00000000d0,   1.00d0,    0.0000d0,   0.00d0,
     1     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     2    28.00000000d0,   2.00d0,    0.0000d0,   0.00d0,
     3     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4    51.00000000d0,   0.00d0,    0.0000d0,   0.00d0,
     5     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     6    52.00000000d0,   0.00d0,    0.0000d0,   0.00d0,
     7     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8    53.00000000d0,   1.00d0,    0.0000d0,   0.00d0,
     9     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *    53.00000000d0,   2.00d0,    0.0000d0,   0.00d0,
     1     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     2    54.00000000d0,   1.00d0,    0.0000d0,   0.00d0,
     3     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     4    54.00000000d0,   2.00d0,    0.0000d0,   0.00d0,
     5     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6    19.00000000d0,   0.00d0,    0.0000d0,   0.00d0,
     7     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     8    20.00000000d0,   0.00d0,    0.0000d0,   0.00d0,
     9     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *    30.00000000d0,   0.00d0,    0.0000d0,   0.00d0,
     1     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     2    31.00000000d0,   1.00d0,    0.0000d0,   0.00d0,
     3     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     4    31.00000000d0,   2.00d0,    0.0000d0,   0.00d0,
     5     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6    39.00000000d0,   0.00d0,    0.0000d0,   0.00d0,
     7     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     8    36.00000000d0,   1.00d0,    0.0000d0,   0.00d0,
     9     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     *    37.00000000d0,   1.00d0,    0.0000d0,   0.00d0,
     1     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     2    37.00000000d0,   2.00d0,    0.0000d0,   0.00d0,
     3     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4    32.00000000d0,   0.00d0,    0.0000d0,   0.00d0,
     5     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     6    36.00000000d0,   2.00d0,    0.0000d0,   0.00d0,
     7     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8    40.00000000d0,   1.00d0,    0.0000d0,   0.00d0,
     9     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *    40.00000000d0,   2.00d0,    0.0000d0,   0.00d0,
     1     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     2    43.00000000d0,   1.00d0,    0.0000d0,   0.00d0,
     3     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     4    43.00000000d0,   2.00d0,    0.0000d0,   0.00d0,
     5     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     6    43.00000000d0,   3.00d0,    0.0000d0,   0.00d0,
     7     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     8    43.00000000d0,   4.00d0,    0.0000d0,   0.00d0,
     9     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *    42.00000000d0,   1.00d0,    0.0000d0,   0.00d0,
     1     1.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2    42.00000000d0,   2.00d0,    0.0000d0,   0.00d0,
     3     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     4    41.00000000d0,   1.00d0,    0.0000d0,   0.00d0,
     5     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     6    41.00000000d0,   2.00d0,    0.0000d0,   0.00d0,
     7     1.00000000d0,   0.00d0,  138.0390d0,   1.00d0,
     8    41.00000000d0,   3.00d0,    0.0000d0,   0.00d0/
c
      data ((tab(i,j),i=1,4),j=79,136)/
     9     0.00000000d0,   0.00d0,    2.0000d0,   1.00d0,
     *     3.58962000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     2    -0.08952900d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     4     0.03600500d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     6    -0.03939300d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     8    11.05911600d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     *    -0.20801800d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     2    -0.00001800d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     4    -0.00002000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     6     2.07316200d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     8    -0.20617600d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     *     0.02129600d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     2    -0.02000700d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  14.00d0,    1.0000d0,   0.70d0,
     4    -0.00000100d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  15.00d0,    1.0000d0,   0.70d0,
     6    -0.50012000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  16.00d0,    1.0000d0,   0.70d0,
     8     0.09631500d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  17.00d0,    1.0000d0,   0.70d0,
     *    -0.00005000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     2    -0.08963000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     4     9.06127400d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     6    -0.23120400d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     8     0.04384700d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     *    -0.00062400d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  14.00d0,    1.0000d0,   0.70d0,
     2     0.00017600d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  15.00d0,    1.0000d0,   0.70d0,
     4    -0.93560000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  16.00d0,    1.0000d0,   0.70d0,
     6    -0.02270700d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  17.00d0,    1.0000d0,   0.70d0,
     8    -0.07670300d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     *     0.03168000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     2     3.67115000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     4     3.64207200d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     6     0.00000000d0,   0.00d0,    0.0000d0,   0.00d0/
c
c
      data ((tabct(2,i,j),i=1,4),j=1,58)/
     9     0.00000000d0,   0.00d0,    2.0000d0,   1.00d0,
     *     1.28494660d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     2     0.13615500d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     4     0.02588600d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     6    -0.04744100d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     8    10.25968000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     *     0.00130100d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     2    -0.00002600d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     4    -0.00002500d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     6     0.37095800d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     8    -0.06413800d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     *     0.02941300d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     2    -0.02003100d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  14.00d0,    1.0000d0,   0.72d0,
     4    -0.00000100d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  15.00d0,    1.0000d0,   0.72d0,
     6    -0.50077900d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  16.00d0,    1.0000d0,   0.72d0,
     8     0.07944200d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  17.00d0,    1.0000d0,   0.72d0,
     *    -0.00003200d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     2    -0.08180500d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     4     5.32851800d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     6    -0.14852000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     8     0.03712800d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     *    -0.00403700d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  14.00d0,    1.0000d0,   0.72d0,
     2    -0.00234700d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  15.00d0,    1.0000d0,   0.72d0,
     4    -0.93977200d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  16.00d0,    1.0000d0,   0.72d0,
     6    -0.01131500d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  17.00d0,    1.0000d0,   0.72d0,
     8    -0.06758900d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     *    -0.00636700d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     2     1.37494000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     4     1.34288800d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     6     0.00000000d0,   0.00d0,    0.0000d0,   0.00d0/
c
c
      data ((tabct(3,i,j),i=1,4),j=1,58)/
     9     0.00000000d0,   0.00d0,    2.0000d0,   1.00d0,
     *    -0.00357370d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     2     0.23877700d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     4     0.03020800d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     6    -0.04646700d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     8    10.94008900d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     *     0.30415700d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     2    -0.00002200d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     4    -0.00003100d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     6    -0.88013000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     8     0.03337800d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     *     0.03757900d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     2    -0.01998300d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  14.00d0,    1.0000d0,   0.75d0,
     4     0.00001400d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  15.00d0,    1.0000d0,   0.75d0,
     6    -0.49720300d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  16.00d0,    1.0000d0,   0.75d0,
     8     0.08328300d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  17.00d0,    1.0000d0,   0.75d0,
     *     0.00077500d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     2    -0.08922200d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     4     4.18430000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     6    -0.08835500d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     8     0.02690600d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     *    -0.00746800d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  14.00d0,    1.0000d0,   0.75d0,
     2    -0.00304300d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  15.00d0,    1.0000d0,   0.75d0,
     4    -0.96172000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  16.00d0,    1.0000d0,   0.75d0,
     6    -0.02242400d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  17.00d0,    1.0000d0,   0.75d0,
     8    -0.05703700d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     *     0.01222800d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     2     0.08355000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     4     0.05583400d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     6     0.00000000d0,   0.00d0,    0.0000d0,   0.00d0/
c
c
      data ((tabct(4,i,j),i=1,4),j=1,58)/
     3     0.00000000d0,   0.00d0,    2.0000d0,   1.00d0,
     4     0.18650650d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     6     0.26070100d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     2    17.37579600d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     4     1.14987300d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     *     0.03592100d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     2     0.18242100d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     6     0.00003300d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  14.00d0,    1.0000d0,   0.70d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  15.00d0,    1.0000d0,   0.70d0,
     *     0.00003200d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  16.00d0,    1.0000d0,   0.70d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  17.00d0,    1.0000d0,   0.70d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     8     5.19812900d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     *     0.07585500d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     4     0.00601200d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  14.00d0,    1.0000d0,   0.70d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  15.00d0,    1.0000d0,   0.70d0,
     8    -1.31461200d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  16.00d0,    1.0000d0,   0.70d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  17.00d0,    1.0000d0,   0.70d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     6     0.26101500d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     8     0.23525000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     *     0.00000000d0,   0.00d0,    0.0000d0,   0.00d0/
c
c
      data ((tabct(5,i,j),i=1,4),j=1,58)/
     3     0.00000000d0,   0.00d0,    2.0000d0,   1.00d0,
     4    -0.80906950d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     6     0.37480400d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     2    17.09651000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     4     1.10425800d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     *    -1.25966640d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     2     0.25515300d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     6     0.00005200d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  14.00d0,    1.0000d0,   0.72d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  15.00d0,    1.0000d0,   0.72d0,
     *    -0.27298500d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  16.00d0,    1.0000d0,   0.72d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  17.00d0,    1.0000d0,   0.72d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     8     3.15881100d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     *     0.03604100d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     4     0.00239100d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  14.00d0,    1.0000d0,   0.72d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  15.00d0,    1.0000d0,   0.72d0,
     8    -1.12995500d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  16.00d0,    1.0000d0,   0.72d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  17.00d0,    1.0000d0,   0.72d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     6    -0.73224900d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     8    -0.75673500d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     *     0.00000000d0,   0.00d0,    0.0000d0,   0.00d0/
c
c
      data ((tabct(6,i,j),i=1,4),j=1,58)/
     3     0.00000000d0,   0.00d0,    2.0000d0,   1.00d0,
     4    -1.52604650d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     6     0.47858400d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     2    20.30045200d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     4     1.52014400d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     *    -2.11088810d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     2     0.37943400d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     6     0.00003200d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  14.00d0,    1.0000d0,   0.75d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  15.00d0,    1.0000d0,   0.75d0,
     *     0.00003100d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  16.00d0,    1.0000d0,   0.75d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  17.00d0,    1.0000d0,   0.75d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     8     2.43712300d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     *     0.08738800d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     4    -0.00366100d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  14.00d0,    1.0000d0,   0.75d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  15.00d0,    1.0000d0,   0.75d0,
     8    -1.14828800d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  16.00d0,    1.0000d0,   0.75d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  17.00d0,    1.0000d0,   0.75d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     6    -1.44536600d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     8    -1.46997300d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     *     0.00000000d0,   0.00d0,    0.0000d0,   0.00d0/
c
c
      data ((tabct(7,i,j),i=1,4),j=1,58)/
     7     0.00000000d0,   0.00d0,    2.0000d0,   1.00d0,
     8    -1.30789200d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     *     0.84552200d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     6    29.62293000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     8     2.06478400d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     4    -1.06169900d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     6     0.67123300d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     *     0.00257400d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  14.00d0,    1.0000d0,   0.70d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  15.00d0,    1.0000d0,   0.70d0,
     4     0.00286100d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  16.00d0,    1.0000d0,   0.70d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  17.00d0,    1.0000d0,   0.70d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     2     4.59197200d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     4     0.62462300d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     8     0.02385000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  14.00d0,    1.0000d0,   0.70d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  15.00d0,    1.0000d0,   0.70d0,
     2    -1.19010200d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  16.00d0,    1.0000d0,   0.70d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  17.00d0,    1.0000d0,   0.70d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     *    -1.21858200d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     2    -1.25084000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     4     0.00000000d0,   0.00d0,    0.0000d0,   0.00d0/
c
c
      data ((tabct(8,i,j),i=1,4),j=1,58)/
     7     0.00000000d0,   0.00d0,    2.0000d0,   1.00d0,
     8    -1.84376100d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     *     0.84348400d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     6    24.45501000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     8     1.85884700d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     4    -2.19819300d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     6     0.62734100d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     *    -0.00018900d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  14.00d0,    1.0000d0,   0.72d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  15.00d0,    1.0000d0,   0.72d0,
     4    -0.10940800d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  16.00d0,    1.0000d0,   0.72d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  17.00d0,    1.0000d0,   0.72d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     2     3.69163300d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     4     0.56754200d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     8     0.01754300d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  14.00d0,    1.0000d0,   0.72d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  15.00d0,    1.0000d0,   0.72d0,
     2    -1.17695500d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  16.00d0,    1.0000d0,   0.72d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  17.00d0,    1.0000d0,   0.72d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     *    -1.75453900d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     2    -1.78439460d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     4     0.00000000d0,   0.00d0,    0.0000d0,   0.00d0/
c
c
      data ((tabct(9,i,j),i=1,4),j=1,58)/
     7     0.00000000d0,   0.00d0,    2.0000d0,   1.00d0,
     8    -2.28683300d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     *     0.87079100d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     6    19.05570000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     8     1.63208700d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     4    -3.07944800d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     6     0.62145400d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     *    -0.00019200d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  14.00d0,    1.0000d0,   0.75d0,
     2     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  15.00d0,    1.0000d0,   0.75d0,
     4    -0.68374800d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  16.00d0,    1.0000d0,   0.75d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  17.00d0,    1.0000d0,   0.75d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     2     3.84533600d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     4     0.64111100d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     8     0.00592200d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  14.00d0,    1.0000d0,   0.75d0,
     *     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  15.00d0,    1.0000d0,   0.75d0,
     2    -1.27680100d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  16.00d0,    1.0000d0,   0.75d0,
     4     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     5     3.00000000d0,  17.00d0,    1.0000d0,   0.75d0,
     6     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     7     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     8     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     9     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     *    -2.19351600d0,   0.00d0,  138.0390d0,   0.00d0,
     1     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     2    -2.22393500d0,   0.00d0,  138.0390d0,   0.00d0,
     3     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     4     0.00000000d0,   0.00d0,    0.0000d0,   0.00d0/
c
c
      data ((tabct(10,i,j),i=1,4),j=1,58)/
     *     0.00000000d0,   0.00d0,    2.0000d0,   1.00d0,
     1    -1.82650880d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     7     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     9     2.60000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     7    -0.92665590d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     9     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  14.00d0,    1.0000d0,   0.70d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  15.00d0,    1.0000d0,   0.70d0,
     7     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  16.00d0,    1.0000d0,   0.70d0,
     9     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  17.00d0,    1.0000d0,   0.70d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     7     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  12.00d0,    1.0000d0,   0.70d0,
     9     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  13.00d0,    1.0000d0,   0.70d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  14.00d0,    1.0000d0,   0.70d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  15.00d0,    1.0000d0,   0.70d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  16.00d0,    1.0000d0,   0.70d0,
     7     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  17.00d0,    1.0000d0,   0.70d0,
     9     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  18.00d0,    1.0000d0,   0.70d0,
     3    -1.85330000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     5    -1.81399000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  11.00d0,    1.0000d0,   0.70d0,
     7     0.00000000d0,   0.00d0,    0.0000d0,   0.00d0/
c
c
      data ((tabct(11,i,j),i=1,4),j=1,58)/
     *     0.00000000d0,   0.00d0,    2.0000d0,   1.00d0,
     1    -1.90320700d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     7     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     9     2.65000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     7    -1.48127000d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     9     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  14.00d0,    1.0000d0,   0.72d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  15.00d0,    1.0000d0,   0.72d0,
     7     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  16.00d0,    1.0000d0,   0.72d0,
     9     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  17.00d0,    1.0000d0,   0.72d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     7     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  12.00d0,    1.0000d0,   0.72d0,
     9     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  13.00d0,    1.0000d0,   0.72d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  14.00d0,    1.0000d0,   0.72d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  15.00d0,    1.0000d0,   0.72d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  16.00d0,    1.0000d0,   0.72d0,
     7     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  17.00d0,    1.0000d0,   0.72d0,
     9     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  18.00d0,    1.0000d0,   0.72d0,
     3    -1.92730000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     5    -1.88746500d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  11.00d0,    1.0000d0,   0.72d0,
     7     0.00000000d0,   0.00d0,    0.0000d0,   0.00d0/
c
c
      data ((tabct(12,i,j),i=1,4),j=1,58)/
     *     0.00000000d0,   0.00d0,    2.0000d0,   1.00d0,
     1    -2.00152800d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     7     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     9     2.63000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     7    -1.91181020d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     9     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  14.00d0,    1.0000d0,   0.75d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  15.00d0,    1.0000d0,   0.75d0,
     7     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  16.00d0,    1.0000d0,   0.75d0,
     9     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  17.00d0,    1.0000d0,   0.75d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     7     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  12.00d0,    1.0000d0,   0.75d0,
     9     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  13.00d0,    1.0000d0,   0.75d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  14.00d0,    1.0000d0,   0.75d0,
     3     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  15.00d0,    1.0000d0,   0.75d0,
     5     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  16.00d0,    1.0000d0,   0.75d0,
     7     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     8     3.00000000d0,  17.00d0,    1.0000d0,   0.75d0,
     9     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     *     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     1     0.00000000d0,   0.00d0,  138.0390d0,   0.00d0,
     2     3.00000000d0,  18.00d0,    1.0000d0,   0.75d0,
     3    -2.02270000d0,   0.00d0,  138.0390d0,   0.00d0,
     4     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     5    -1.98234600d0,   0.00d0,  138.0390d0,   0.00d0,
     6     3.00000000d0,  11.00d0,    1.0000d0,   0.75d0,
     7     0.00000000d0,   0.00d0,    0.0000d0,   0.00d0/
c
c
c        this has been the end of tables
c
c
10000 format (2a4,a2,15a4)
10003 format (' ',34(1h-))
10006 format (' ',2a4,a2,f12.8,f10.6,1x,f10.5,f7.1)
10007 format (' ',2a4,a2,f4.1,1x,2f10.5,f13.5)
10010 format (2a4,a2,4f10.8)
10020 format (' ',2a4,a2,15a4)
10021 format (' ',2a4,a2,4f10.6)
10022 format (' ',2a4,a2,f10.4,f10.2)
c
c
c
c
      if (index) go to 50
      index=.true.
c
c
      ist=0
c
      ime=0
      imee=300
      ile=20
      mge=24
c
c         set all parameters and indices to zero or .false.
c
      do 1 imm=1,imee
      mgg(imm)=0
      indt(imm)=.false.
c
      do 1 il=1,ile
      c(il,imm)=0.d0
    1 ic(il,imm)=0
c
      indca=.false.
c
      sqrpi=dsqrt(pi)
      pi2=pi*pi
      pi3=pi2*pi
      pi4=pi3*pi
      pi5=pi4*pi
      pi32=dsqrt(pi3)
c
c
c
   50 continue
c
c         nucleon mass
c
      wn=938.9183d0
c
c         ga and fpi
c
      ga=1.29d0
      fpi=92.4d0
c
      ga2 = ga * ga
      ga4 = ga2 * ga2
      fpi2 = fpi * fpi
      fpi4 = fpi2 * fpi2
      fpi6 = fpi4 * fpi2
c
c
c          the pi-N LECs
c          -------------
c
      if (lpot.le.3) then
c
c         N3LO LECs
c
c         the c_i LECs
      cc(1)=-1.07d0
      cc(2)=3.20d0
      cc(3)=-5.32d0
      cc(4)=3.56d0
c
      cb1=cc(1)*1.d-3
      cb2=cc(2)*1.d-3
      cb3=cc(3)*1.d-3
      cb4=cc(4)*1.d-3
c
c         the d_i LECs
      cc(1)=1.04d0
      cc(2)=-0.48d0
      cc(3)=0.14d0
      cc(4)=-1.90d0
c
      cd12=cc(1)*1.d-6
      cd3=cc(2)*1.d-6
      cd5=cc(3)*1.d-6
      cd145=cc(4)*1.d-6
c
      else
c
c         NNLO LECs
c
c         the c_i LECs
      cc(1)=-0.74d0
      cc(2)=0.0d0
      cc(3)=-3.61d0
      cc(4)=2.44d0
c
      cb1=cc(1)*1.d-3
      cb2=cc(2)*1.d-3
      cb3=cc(3)*1.d-3
      cb4=cc(4)*1.d-3
c
c         the d_i LECs
      cc(1)=0.0d0
      cc(2)=0.0d0
      cc(3)=0.0d0
      cc(4)=0.0d0
c
      cd12=cc(1)*1.d-6
      cd3=cc(2)*1.d-6
      cd5=cc(3)*1.d-6
      cd145=cc(4)*1.d-6
c
      end if
c
c
c         prepare tables
c
      do 55 ii=1,136
      do 55 i=1,4
   55 ttab(i,ii)=tab(i,ii)
c
      icut=mod(lpot,3)
      if (icut.eq.0) icut=3
c
      go to (51,52,53),icut
   51 continue
c         Rpi=1.0 is default
      go to 54
   52 ttab(4,1)=1.1d0
      ttab(4,10)=1.1d0
      go to 54
   53 ttab(4,1)=1.2d0
      ttab(4,10)=1.2d0
   54 continue
c
      if (lpot.ne.1) then
      do 56 ii=1,58
      do 56 i=1,4
   56 ttab(i,ii+78)=tabct(lpot,i,ii)
      end if
c
c
c         get parameters from tables, line by line
c         ----------------------------------------
c         ----------------------------------------
c
c
      line=0
c
   61 line=line+1
c
c
c         NNLO
      if (lpot.ge.4.and.lpot.le.6.and.line.eq.23) line=79
c         NLO
      if (lpot.ge.7.and.lpot.le.9.and.line.eq.17) line=79
c         LO
      if (lpot.ge.10.and.lpot.le.12.and.line.eq.10) line=79
c
c
      do i=1,4
      if (i.le.3) then
      name(i)=ntab(i,line)
      end if
      cc(i)=ttab(i,line)
      end do

c
c         check if data card just read contains cutoff or function parameters
c
      if (name(1).eq.cut.or.name(1).eq.fun) go to 70
c
c         check if end
c
      if (name(1).eq.end) go to 8500
c
c         check if there is a 'cutall'
c
      if (name(1).eq.cuta) then
c**** write (kwrite,10007) name,cc
      indca=.true.
      do i=1,4
      cca(i)=cc(i)
      end do
      go to 61
      end if
c
c
c
c        write parameters which are no cutoff or function parameters
c        -----------------------------------------------------------
c
c
c**** write (kwrite,10006) name,cc
c
c
c        find out type of contribution
c
      do 63 mg=1,mge
      if (name(1).eq.opsym(mg)) go to 64
   63 continue
      go to 9000
c
c
c
c         store contribution parameters which are no cutoff or
c         ----------------------------------------------------
c         function parameters
c         -------------------
c
   64 ime=ime+1
      if (ime.gt.imee) go to 9011
      mgg(ime)=mg
c
      if (mg.gt.7) then
      ist=ist+1
      cst(ist)=cc(1)
      end if
c
c
c         store coupling constant
      c(1,ime)=cc(1)
      if (cc(2).ne.0.d0) then
      c(1,ime)=(cc(1)/cc(2))**2
      end if
c         store meson mass
      c(4,ime)=cc(3)
c         test and store isospin as logical parameter
      icc=cc(4)
      if (icc.lt.0.or.icc.gt.1) go to 9004
      if (icc.eq.1) indt(ime)=.true.
      mi=1
      mm=5
c
c        check if there is a `cutall' cutoff
c
      if (indca) then
      do i=1,4
      cc(i)=cca(i)
      end do
      go to 72
      else
      go to 61
      end if
c
c
c
c         write cutoff or function parameters
c         -----------------------------------
c
c
   70 continue
c**** write (kwrite,10007) name,cc
   72 continue
c
c
c
c
c          store cutoff or function parameters
c          -----------------------------------
c
      itype=cc(1)
      if (itype.eq.0) go to 8000
      if (itype.lt.0.or.itype.gt.54) go to 9002
      ic(mi,ime)=itype
      icase=cc(2)
      ic(mi+1,ime)=icase
      c4=c(4,ime)
      c47=c4*c4*c4*c4*c4*c4*c4
c
c
      go to(100,200,300,9002,9002,9002,700,800,900,9002,
     1 9002,9002,1300,1400,1500,1600,1700,1800,1900,2000,
     2 9002,9002,9002,9002,2500,2600,9002,2800,2900,3000,
     3 3100,3200,9002,9002,9002,3600,3700,9002,3900,4000,
     4 4100,4200,4300,9002,9002,9002,9002,9002,9002,9002,
     5 5100,5200,5300,5400),itype
c
c
c        for 1PE
c        -------
  100 c(mm,ime)=c4*c4/(48.d0*pi)
      go to 95
c
c        for 1PE
c        -------
  200 c(mm,ime)=1.d0/(48.d0*pi)
      go to 95
c
c
c
  300 c(mm+1,ime)=cc(4)
      nexp=cc(3)
      ic(mi+2,ime)=nexp
      if (nexp.lt.1.or.nexp.gt.2) go to 9003
      go to (320,9002),nexp
  320 if (icase.lt.1.or.icase.gt.18) go to 9002
      go to (321,322,9002,9002,325,326,327,9002,9002,9002,
     1        331,332,333,334,335,336,337,338),icase
c
  321 c(mm,ime)=1.d0
      go to 95
c
  322 c(mm,ime)=1.d0
      go to 95
c
  325 c(mm,ime)=1.d0
      go to 95
c
  326 c(mm,ime)=1.d0
      go to 95
c
  327 c(mm,ime)=1.d0
      go to 95
c
  331 c(mm,ime)=hbarc/pi32
      go to 95
  332 c(mm,ime)=hbarc/pi32
      go to 95
  333 c(mm,ime)=hbarc/pi32
      go to 95
  334 c(mm,ime)=hbarc/pi32
      go to 95
  335 c(mm,ime)=hbarc/pi32
      go to 95
  336 c(mm,ime)=hbarc/pi32
      go to 95
  337 c(mm,ime)=hbarc/pi32
      go to 95
  338 c(mm,ime)=hbarc/pi32
      go to 95
c
c
  700 c(mm,ime)=cc(3)
      c(mm+1,ime)=cc(4)
      go to 95
c
  800 c(mm,ime)=cc(3)
      c(mm+1,ime)=cc(4)
      go to 95
c
c        Piarulli's pi ff
  900 c(mm,ime)=cc(3)
      c(mm+1,ime)=cc(4)
      go to 95
c
 1300 c(mm,ime)=c4/(128.d0*pi3*fpi4)
      go to 95
c
 1400 if (icase.lt.1.or.icase.gt.2) go to 9002
      go to (1410,1420),icase
 1410 c(mm,ime)=ga4*c4/(32.d0*pi3*fpi4)
      go to 95
c
 1420 c(mm,ime)=ga4*c4/(128.d0*pi3*fpi4)
      go to 95
c
 1500 c(mm,ime)=3.d0*ga2/(32.d0*pi2*fpi4)
      go to 95
c
 1600 c(mm,ime)=c4*ga2/(512.d0*pi2*fpi4*wn)
      go to 95
c
 1700 if (icase.lt.1.or.icase.gt.2) go to 9002
      go to (1710,1720),icase
 1710 c(mm,ime)=-3.d0*ga4*c4/(512.d0*pi2*fpi4*wn)
      go to 95
c
 1720 c(mm,ime)=3.d0*ga4*c4/(1024.d0*pi2*fpi4*wn)
      go to 95
c
 1800 if (icase.lt.1.or.icase.gt.2) go to 9002
      go to (1810,1820),icase
 1810 c(mm,ime)=ga2/(48.d0*pi2*fpi4)
      go to 95
c
 1820 c(mm,ime)=ga2/(48.d0*pi2*fpi4)
      go to 95
c
 1900 c(mm,ime)=-3.d0*ga4/(64.d0*pi2*wn*fpi4)
      go to 95
c
 2000 c(mm,ime)=ga2*(ga2-1.d0)/(32.d0*pi2*wn*fpi4)
      go to 95
c
 2500 c(mm,ime)=3.d0*ga2/(32.d0*pi2*fpi4)
      go to 95
c
 2600 c(mm,ime)=3.d0*ga2/(32.d0*pi2*fpi4)
      go to 95
c
 2800 if (icase.lt.1.or.icase.gt.2) go to 9002
      go to (2810,2820),icase
 2810 c(mm,ime)=ga2/(48.d0*pi2*fpi4)
      go to 95
c
 2820 c(mm,ime)=ga2/(48.d0*pi2*fpi4)
      go to 95
c
 2900 if (icase.lt.1.or.icase.gt.2) go to 9002
      go to (2910,2920),icase
 2910 c(mm,ime)=ga2/(48.d0*pi2*fpi4)
      go to 95
c
 2920 c(mm,ime)=ga2/(48.d0*pi2*fpi4)
      go to 95
c
 3000 c(mm,ime)=-3.d0*c47/(32.d0*pi3*fpi4)
      go to 95
c
 3100 if (icase.lt.1.or.icase.gt.2) go to 9002
      go to (3110,3120),icase
 3110 c(mm,ime)=cb4*cb4*c47/(24.d0*pi3*fpi4)
      go to 95

 3120 c(mm,ime)=-cb4*cb4*c47/(96.d0*pi3*fpi4)
      go to 95
c
 3200 c(mm,ime)=3.d0*cb2*ga2*c47/(8.d0*pi3*wn*fpi4)
      go to 95
c
 3600 if (icase.lt.1.or.icase.gt.2) go to 9002
      go to (3610,3620),icase
 3610 c(mm,ime)=cb4*c47/(32.d0*pi3*wn*fpi4)
      go to 95
c
 3620 c(mm,ime)=-cb4*c47/(16.d0*pi3*wn*fpi4)
      go to 95
c
 3700 if (icase.lt.1.or.icase.gt.2) go to 9002
      go to (3710,3720),icase
 3710 c(mm,ime)=cb4*c47/(48.d0*pi3*wn*fpi4)
      go to 95

 3720 c(mm,ime)=cb4*c47/(192.d0*pi3*wn*fpi4)
      go to 95
c
 3900 c(mm,ime)=3.d0*ga2*c47/(32.d0*pi3*wn*fpi4)
      go to 95
c
 4000 if (icase.lt.1.or.icase.gt.2) go to 9002
      go to (4010,4020),icase
c
 4010 c(mm,ime)=3.d0*c47*ga4/(2048.d0*pi3*fpi6)
      go to 95
c
 4020 c(mm,ime)=-3.d0*ga4*c47/(8192.d0*pi3*fpi6)
      go to 95
c
 4100 if (icase.lt.1.or.icase.gt.3) go to 9002
      go to (4110,4120,4130),icase
 4110 c(mm,ime)=(c47*ga4)/(6144.d0*pi3*fpi6)
      go to 95
c
 4120 c(mm,ime)=(c47*ga4*(1.d0+2.d0*ga2))/(1536.d0*pi3*fpi6)
      go to 95
c
 4130 c(mm,ime)=-(c47*ga4)/(49152.d0*pi3*fpi6)
      go to 95
c
 4200 if (icase.lt.1.or.icase.gt.2) go to 9002
      go to (4210,4220),icase
 4210 c(mm,ime)=-(c47*ga2*cd145)/(8.d0*pi3*fpi4)
      go to 95
c
 4220 c(mm,ime)=(c47*ga2*cd145)/(32.d0*pi3*fpi4)
      go to 95
c
 4300 if (icase.lt.1.or.icase.gt.4) go to 9002
      go to (4310,4320,4330,4340),icase
 4310 c(mm,ime)=-c47/(9216.d0*pi5*fpi6)
      go to 95
c
 4320 c(mm,ime)=-c47/(8.d0*pi3*fpi4)
      go to 95
c
 4330 c(mm,ime)=c47/(9216.d0*pi5*fpi6)
      go to 95
c
 4340 c(mm,ime)=-c47/(16.d0*pi3*fpi4)
      go to 95
c
 5100 c(mm,ime)=3*ga4/(1024.d0*pi2*fpi4*wn)
      go to 95
c
 5200 c(mm,ime)=ga2/(512.d0*pi2*fpi4*wn)
      go to 95
c
 5300 if (icase.lt.1.or.icase.gt.2) go to 9002
      go to (5310,5320),icase
 5310 c(mm,ime)=-ga4/(512.d0*pi2*fpi4*wn)
      go to 95
c
 5320 c(mm,ime)=ga4/(2048.d0*pi2*fpi4*wn)
      go to 95
c
 5400 if (icase.lt.1.or.icase.gt.2) go to 9002
      go to (5410,5420),icase
 5410 c(mm,ime)=-ga2/(1536.d0*pi2*fpi4*wn)
      go to 95
c
 5420 c(mm,ime)=ga2/(3072.d0*pi2*fpi4*wn)
      go to 95
c
c
   95 mi=mi+3
      mm=mm+3
      if(mi.gt.ile.or.mm.gt.ile) go to 9011
 8000 continue
      go to 61
c
c
c         write end
c         ---------
c
c
 8500 continue
c**** write (kwrite,10020) name
c**** write (kwrite,10003)
c**** write (kwrite,10003)
c
c
      c01np=cst(1)
      c01pp=cst(27)
      c01nn=cst(28)
c
      c01cd=(0.5d0*(c01pp+c01nn)-c01np)/6.d0
      c01ca=(c01pp-c01nn)/4.d0
      cvw(27)=c01cd/4.d0
      cvw(28)=-cvw(27)
      cvw(29)=c01ca/4.d0
      cvw(30)=-cvw(29)
c
c        the charge-independent c01
      cst(1)=(c01pp+c01nn+c01np)/3.d0
c
      cvw(1)=(3.d0*cst(1)+cst(5)+3.d0*cst(9)+9.d0*cst(18))/16.d0
      cvw(2)=(1.d0*cst(1)-cst(5)-3.d0*cst(9)+3.d0*cst(18))/16.d0
      cvw(3)=(-3.d0*cst(1)-1.d0*cst(5)+cst(9)+3.d0*cst(18))/16.d0
      cvw(4)=(-1.d0*cst(1)+1.d0*cst(5)-cst(9)+1.d0*cst(18))/16.d0
c
      cvw(5)=(3.d0*cst(2)+cst(6)+3.d0*cst(10)+9.d0*cst(19))/16.d0
      cvw(6)=(1.d0*cst(2)-cst(6)-3.d0*cst(10)+3.d0*cst(19))/16.d0
      cvw(7)=(-3.d0*cst(2)-1.d0*cst(6)+cst(10)+3.d0*cst(19))/16.d0
      cvw(8)=(-1.d0*cst(2)+1.d0*cst(6)-cst(10)+1.d0*cst(19))/16.d0
c
      cvw(9)=(3.d0*cst(3)+cst(7)+3.d0*cst(11)+9.d0*cst(20))/16.d0
      cvw(10)=(1.d0*cst(3)-cst(7)-3.d0*cst(11)+3.d0*cst(20))/16.d0
      cvw(11)=(-3.d0*cst(3)-1.d0*cst(7)+cst(11)+3.d0*cst(20))/16.d0
      cvw(12)=(-1.d0*cst(3)+1.d0*cst(7)-cst(11)+1.d0*cst(20))/16.d0
c
      cvw(13)=(cst(12)+3.d0*cst(21))/4.d0
      cvw(14)=(-cst(12)+cst(21))/4.d0
c
      cvw(15)=(cst(13)+3.d0*cst(22))/4.d0
      cvw(16)=(-cst(13)+cst(22))/4.d0
c
      cvw(17)=(cst(14)+3.d0*cst(23))/4.d0
      cvw(18)=(-cst(14)+cst(23))/4.d0
c
      cvw(19)=(cst(15)+3.d0*cst(24))/4.d0
      cvw(20)=(-cst(15)+cst(24))/4.d0
c
      cvw(21)=(cst(16)+3.d0*cst(25))/4.d0
      cvw(22)=(-cst(16)+cst(25))/4.d0
c
      cvw(23)=(3.d0*cst(4)+cst(8)+3.d0*cst(17)+9.d0*cst(26))/16.d0
      cvw(24)=(1.d0*cst(4)-cst(8)-3.d0*cst(17)+3.d0*cst(26))/16.d0
      cvw(25)=(-3.d0*cst(4)-1.d0*cst(8)+cst(17)+3.d0*cst(26))/16.d0
      cvw(26)=(-1.d0*cst(4)+1.d0*cst(8)-cst(17)+1.d0*cst(26))/16.d0
c
c
c
c**** write (kwrite,10030) (cvw(i),i=1,4),
c****1                     (cvw(i2),i2=5,8),
c****3                     (cvw(i4),i4=13,14),
c****5                     (cvw(i6),i6=17,18),
c****2                     (cvw(i3),i3=9,12),
c****4                     (cvw(i5),i5=15,16),
c****6                     (cvw(i7),i7=19,20),
c****7                     (cvw(i8),i8=21,22),
c****8                     (cvw(i9),i9=23,26),
c****9                     (cvw(i0),i0=27,30)
10030 format(//
     *       ' Q0: Cc  Ct  Cs   Cst ',4e16.8/
     1       ' Q2: C1  C2  C3   C4  ',4e16.8/
     3       ' Q2: C5      C6       ',2e16.8/
     5       ' Q2: C7      C8       ',2e16.8/
     2       ' Q4: D1  D2  D3   D4  ',4e16.8/
     4       ' Q4: D5      D6       ',2e16.8/
     6       ' Q4: D7      D8       ',2e16.8/
     7       ' Q4: D9      D10      ',2e16.8/
     8       ' Q4: D11,D12,D13, D14 ',4e16.8/
     9       ' Q0: CT12,CsT12,Ctz,Cstz ',4e16.8//)
c
c
      return
c
c
c         errors
c         ------
c         ------
c
c
 9000 write (kwrite,19000) name(1)
19000 format(/////' error in chrpar: contribution   ',a4,'   does not
     1 exist in this program.'/' execution terminated.'////)
      go to 9999
c
c
 9002 write (kwrite,19002) cc(1),cc(2)
19002 format (/////' error in chrpar: cut/fun type',f10.4,'  and icase',
     1f10.4,' does not exist in this program.'/' execution terminated.'
     2////)
      go to 9999
c
c
 9003 write (kwrite,19003) cc(3)
19003 format (/////' error in chrpar:   nexp  has the non-permissible
     1value',f10.4,'  .'/' execution terminated.'////)
      go to 9999
c
c
 9004 write (kwrite,19004) cc(4)
19004 format (/////' error in chrpar: isospin has the non-permissible
     1value',f10.4,'  .'/' execution terminated.'////)
      go to 9999
c
c
 9011 write (kwrite,19011)
19011 format (/////' error in chrpar: too many contributions with respec
     1t to the dimensions given'/' in this program. execution terminated
     2.'////)
      go to 9999
c
c
 9999 stop
      end
c
c
      subroutine opera
c
c        operator matrices
c
c        The first index (k) of opmat(k,i) runs from 1 to 5
c        and refers, for fixed j
c        and using the notation v(s;l,l'),
c        to the following matrix elements:
c        v(0;j,j) , v(1;j,j) , v(1;j-1,j-1) , v(1;j+1,j+1) ,
c        v(1;j-1,j+1)
c
c        The second index (i) runs from 1 to 7 and denotes
c        the various operators, with
c        1   c    (central)
c        2   ss   (sigma1 dot sigma2)
c        3   s12  (tensor)
c        4   ls   (spin-orbit)
c        5   l2   (L square)
c        6   l2ss (L square times sigma1 dot sigma2)
c        7   ls2  (quadratic spin-orbit)
c
      implicit real*8 (a-h,o-z)
      common/copera/opmat(5,7),jj
c
c
      x = dfloat(jj)
      x1 = x + 1.d0
      x2 = 2.d0 * x + 1.d0
c
c
      do 8 i=1,7
      do 8 k=1,5
    8 opmat(k,i) = 0.d0
c
c        central and spin-spin operators
      do 10 i=1,2
      do 10 k=1,4
   10 opmat(k,i) = 1.d0
      opmat(1,2) = -3.d0
c
c        spin-orbit operator
      opmat(2,4) = -1.d0
      opmat(3,4) = x - 1.d0
      opmat(4,4) = - ( x + 2.d0 )
c
c        tensor operator
      opmat(2,3) = 2.d0
      opmat(3,3) = -2.d0 * opmat(3,4) / x2
      opmat(4,3) =  2.d0 * opmat(4,4) / x2
      opmat(5,3) = dsqrt(x*x1) * 6.d0 / x2
c
c
c        l**2 operator
      opmat(1,5)=x*x1
      opmat(2,5)=x*x1
      opmat(3,5)=x*opmat(3,4)
      opmat(4,5)=x1*(x+2.d0)
c
c
c        l2ss  (l square times sigma1 dot sigma2)
      do 15 k=1,4
   15 opmat(k,6)=opmat(k,5)*opmat(k,2)
c
c
c        ls**2 operator
      opmat(2,7) = opmat(2,4) * opmat(2,4)
      opmat(3,7) = opmat(3,4) * opmat(3,4)
      opmat(4,7) = opmat(4,4) * opmat(4,4)
c
c
c        the case j=0
c
      if (jj.eq.0) then
      do 25 i=1,7
      opmat(2,i) = 0.d0
      opmat(3,i) = 0.d0
   25 opmat(5,i) = 0.d0
      end if
c
c
      return
      end
c
c
c $Id: besk0.F,v 1.1.1.1 1996/02/15 17:49:09 mclareni Exp $
c
c $Log: besk0.F,v $
c Revision 1.1.1.1  1996/02/15 17:49:09  mclareni
c Kernlib
c
c
      function dbesk0(dx)
      double precision x,y,r,a,a0,a1,a2,b,b0,b1,b2,t(10)
      double precision u0,u1,u2,u3,u4,u5,u6,u7,u8,u9
      double precision f,f1,f2,f3,c,c0,pi1,ce,eps,h,alfa,d
      double precision zero,one,two,four,five,six,seven,eight,nine,half
      double precision c1(0:14),c2(0:15),c3(0:12)
      double precision dbesk0,dx

      data zero /0.0d0/, one /1.0d0/, two /2.0d0/
      data four /4.0d0/, five /5.0d0/, six /6.0d0/, seven /7.0d0/
      data eight /8.0d0/, nine /9.0d0/, half /0.5d0/

      data t /16.0d0,368.0d0,43.0d0,75.0d0,400.0d0,40.0d0,
     1        48.0d0,12.0d0,20.0d0,28.0d0/

      data pi1 /1.25331 41373 155d0/, ce /0.57721 56649 0153d0/
      data eps /1.0d-14/

      data c1( 0) /0.12773 34398 1218d3/
      data c1( 1) /0.19049 43201 7274d3/
      data c1( 2) /0.82489 03274 4024d2/
      data c1( 3) /0.22274 81924 2462d2/
      data c1( 4) /0.40116 73760 1793d1/
      data c1( 5) /0.50949 33654 3998d0/
      data c1( 6) /0.04771 87487 9817d0/
      data c1( 7) /0.00341 63317 6601d0/
      data c1( 8) /0.00019 24693 5969d0/
      data c1( 9) /0.00000 87383 1550d0/
      data c1(10) /0.00000 03260 9105d0/
      data c1(11) /0.00000 00101 6973d0/
      data c1(12) /0.00000 00002 6883d0/
      data c1(13) /0.00000 00000 0610d0/
      data c1(14) /0.00000 00000 0012d0/

      data c2( 0) /0.24027 70596 4072d3/
      data c2( 1) /0.36947 40739 7287d3/
      data c2( 2) /0.16997 34116 9840d3/
      data c2( 3) /0.49020 46377 7263d2/
      data c2( 4) /0.93884 97325 2684d1/
      data c2( 5) /0.12594 79763 6677d1/
      data c2( 6) /0.12377 69641 1492d0/
      data c2( 7) /0.00924 43098 6287d0/
      data c2( 8) /0.00054 06238 9649d0/
      data c2( 9) /0.00002 53737 9603d0/
      data c2(10) /0.00000 09754 7830d0/
      data c2(11) /0.00000 00312 4957d0/
      data c2(12) /0.00000 00008 4643d0/
      data c2(13) /0.00000 00000 1963d0/
      data c2(14) /0.00000 00000 0039d0/
      data c2(15) /0.00000 00000 0001d0/

      data c3( 0) /+0.98840 81742 3083d0/
      data c3( 1) /-0.01131 05046 4693d0/
      data c3( 2) /+0.00026 95326 1276d0/
      data c3( 3) /-0.00001 11066 8520d0/
      data c3( 4) /+0.00000 06325 7511d0/
      data c3( 5) /-0.00000 00450 4734d0/
      data c3( 6) /+0.00000 00037 9300d0/
      data c3( 7) /-0.00000 00003 6455d0/
      data c3( 8) /+0.00000 00000 3904d0/
      data c3( 9) /-0.00000 00000 0458d0/
      data c3(10) /+0.00000 00000 0058d0/
      data c3(11) /-0.00000 00000 0008d0/
      data c3(12) /+0.00000 00000 0001d0/

      round(d)  =  sngl(d+(d-dble(sngl(d))))

      x=dx

    9 if(x .le. zero) then
       print *,'dbsk0 out of range'
       stop
      endif
      if(x .lt. half) then
       y=x/eight
       h=two*y**2-one
       alfa=-two*h
       b1=zero
       b2=zero
       do 1 i = 14,0,-1
       b0=c1(i)-alfa*b1-b2
       b2=b1
    1  b1=b0
       r=b0-h*b2
       b1=zero
       b2=zero
       do 2 i = 15,0,-1
       b0=c2(i)-alfa*b1-b2
       b2=b1
    2  b1=b0
       b1=-(ce+log(half*x))*r+b0-h*b2
      else if(x .gt. five) then
       r=one/x
       y=five*r
       h=two*y-one
       alfa=-two*h
       b1=zero
       b2=zero
       do 3 i = 12,0,-1
       b0=c3(i)-alfa*b1-b2
       b2=b1
    3  b1=b0
       b1=pi1*sqrt(r)*(b0-h*b2)
       b1=exp(-x)*b1
      else
       y=(t(1)*x)**2
       a0=one
       a1=(t(1)*x+seven)/nine
       a2=(y+t(2)*x+t(3))/t(4)
       b0=one
       b1=(t(1)*x+nine)/nine
       b2=(y+t(5)*x+t(4))/t(4)
       u1=one
       u4=t(6)
       u5=t(7)
       c=zero
       f=two
    4  c0=c
       f=f+one
       u0=t(8)*f**2-one
       u1=u1+two
       u2=u1+two
       u3=u1+four
       u4=u4+t(9)
       u5=u5+t(10)
       u6=one/u3**2
       u7=u2*u6
       u8=-u7/u1
       u9=t(1)*u7*x
       f1=u9-(u0-u4)*u8
       f2=u9-(u0-u5)*u6
       f3=-u8*(u3-six)**2
       a=f1*a2+f2*a1+f3*a0
       b=f1*b2+f2*b1+f3*b0
       c=a/b
       if(abs((c0-c)/c) .ge. eps) then
        a0=a1
        a1=a2
        a2=a
        b0=b1
        b1=b2
        b2=b
        go to 4
       endif
       b1=pi1*c/sqrt(x)
       b1=exp(-x)*b1
      endif

      dbesk0=b1

      return
      end
c
c
c $Id: besk1.F,v 1.1.1.1 1996/02/15 17:49:09 mclareni Exp $
c
c $Log: besk1.F,v $
c Revision 1.1.1.1  1996/02/15 17:49:09  mclareni
c Kernlib
c
c
      function dbesk1(dx)
      double precision x,y,r,a,a0,a1,a2,b,b0,b1,b2,t(12)
      double precision u0,u1,u2,u3,u4,u5,u6,u7,u8,u9
      double precision f,f1,f2,f3,c,c0,pi1,ce,eps,h,alfa,d
      double precision zero,one,two,three,four,five,six,eight,half
      double precision c1(0:14),c2(0:14),c3(0:11)
      double precision dbesk1,dx

      data zero /0.0d0/, one /1.0d0/, two /2.0d0/, three /3.0d0/
      data four /4.0d0/, five /5.0d0/, six /6.0d0/, eight /8.0d0/
      data half /0.5d0/

      data t /16.0d0,3.2d0,2.2d0,432.0d0,131.0d0,35.0d0,336.0d0,
     1        40.0d0,48.0d0,12.0d0,20.0d0,28.0d0/

      data pi1 /1.25331 41373 155d0/, ce /0.57721 56649 0153d0/
      data eps /1.0d-14/

      data c1( 0) /0.22060 14269 2352d3/
      data c1( 1) /0.12535 42668 3715d3/
      data c1( 2) /0.42865 23409 3128d2/
      data c1( 3) /0.94530 05229 4349d1/
      data c1( 4) /0.14296 57709 0762d1/
      data c1( 5) /0.15592 42954 7626d0/
      data c1( 6) /0.01276 80490 8173d0/
      data c1( 7) /0.00081 08879 0069d0/
      data c1( 8) /0.00004 10104 6194d0/
      data c1( 9) /0.00000 16880 4220d0/
      data c1(10) /0.00000 00575 8695d0/
      data c1(11) /0.00000 00016 5345d0/
      data c1(12) /0.00000 00000 4048d0/
      data c1(13) /0.00000 00000 0085d0/
      data c1(14) /0.00000 00000 0002d0/

      data c2( 0) /0.41888 94461 6640d3/
      data c2( 1) /0.24989 55490 4287d3/
      data c2( 2) /0.91180 31933 8742d2/
      data c2( 3) /0.21444 99505 3962d2/
      data c2( 4) /0.34384 15392 8805d1/
      data c2( 5) /0.39484 60929 4094d0/
      data c2( 6) /0.03382 87455 2688d0/
      data c2( 7) /0.00223 57203 3417d0/
      data c2( 8) /0.00011 71310 2246d0/
      data c2( 9) /0.00000 49754 2712d0/
      data c2(10) /0.00000 01746 0493d0/
      data c2(11) /0.00000 00051 4329d0/
      data c2(12) /0.00000 00001 2890d0/
      data c2(13) /0.00000 00000 0278d0/
      data c2(14) /0.00000 00000 0005d0/

      data c3( 0) /+1.03595 08587 724d0/
      data c3( 1) /+0.03546 52912 433d0/
      data c3( 2) /-0.00046 84750 282d0/
      data c3( 3) /+0.00001 61850 638d0/
      data c3( 4) /-0.00000 08451 720d0/
      data c3( 5) /+0.00000 00571 322d0/
      data c3( 6) /-0.00000 00046 456d0/
      data c3( 7) /+0.00000 00004 354d0/
      data c3( 8) /-0.00000 00000 458d0/
      data c3( 9) /+0.00000 00000 053d0/
      data c3(10) /-0.00000 00000 007d0/
      data c3(11) /+0.00000 00000 001d0/

      round(d)  =  sngl(d+(d-dble(sngl(d))))

      x=dx

    9 if(x .le. zero) then
       print *,'dbesk1 out of range'
       stop
      endif
      if(x .lt. half) then
       y=x/eight
       h=two*y**2-one
       alfa=-two*h
       b1=zero
       b2=zero
       do 1 i = 14,0,-1
       b0=c1(i)-alfa*b1-b2
       b2=b1
    1  b1=b0
       r=y*(b0-b2)
       b1=zero
       b2=zero
       do 2 i = 14,0,-1
       b0=c2(i)-alfa*b1-b2
       b2=b1
    2  b1=b0
       b1=(ce+log(half*x))*r+one/x-y*(b0-b2)
      else if(x .gt. five) then
       r=one/x
       y=five*r
       h=two*y-one
       alfa=-two*h
       b1=zero
       b2=zero
       do 3 i = 11,0,-1
       b0=c3(i)-alfa*b1-b2
       b2=b1
    3  b1=b0
       b1=pi1*sqrt(r)*(b0-h*b2)
       b1=exp(-x)*b1
      else
       y=(t(1)*x)**2
       a0=one
       a1=t(2)*x+t(3)
       a2=(y+t(4)*x+t(5))/t(6)
       b0=one
       b1=t(2)*x+one
       b2=(y+t(7)*x+t(6))/t(6)
       u1=one
       u4=t(8)
       u5=t(9)
       c=zero
       f=two
    4  c0=c
       f=f+one
       u0=t(10)*f**2+three
       u1=u1+two
       u2=u1+two
       u3=u1+four
       u4=u4+t(11)
       u5=u5+t(12)
       u6=one/(u3**2-four)
       u7=u2*u6
       u8=-u7/u1
       u9=t(1)*u7*x
       f1=u9-(u0-u4)*u8
       f2=u9-(u0-u5)*u6
       f3=u8*(four-(u3-six)**2)
       a=f1*a2+f2*a1+f3*a0
       b=f1*b2+f2*b1+f3*b0
       c=a/b
       if(abs((c0-c)/c) .ge. eps) then
        a0=a1
        a1=a2
        a2=a
        b0=b1
        b1=b2
        b2=b
        go to 4
       endif
       b1=pi1*c/sqrt(x)
       b1=exp(-x)*b1
      endif

      dbesk1=b1

      return

      end

	function dbesk2(dx)
	   implicit real*8 (a-h, o-z)

	dbesk2= dbesk0(dx)+(2.d0/dx)*dbesk1(dx)
                        return
                        end

	function dbesk3(dx)
	   implicit real*8 (a-h, o-z)

	dbesk3= dbesk1(dx)+(4.d0/dx)*dbesk2(dx)
                        return
                        end

c
c
      function oneme(x)
c
c        oneme(x) = 1 - exp(-x)
c
c        NOTE: this code is useful only for
c              0 < x < 1.
c
c
      implicit real*8 (a-h,o-z)
c
c
      i=0
      term=1.d0
      sum0=0.d0
      sum1=1.d0
c
c
   10 i=i+1
      term=term*x/dfloat(i)
      sum0=sum0+term
      sum1=sum1+term
      if(dabs(term/sum0).gt.1.d-16) go to 10
c
c
      oneme=sum0/sum1
c
c
      return
      end
        subroutine gaussid(ax,bx,n,x,w)
        implicit double precision (a-h,o-z)
        parameter(pi=3.1415926535897932d0)
        dimension x(n),w(n)

        eps = 1.d-15
        do 2 i=1,n
        z=dcos(pi*(i-.25d0)/(n+.5d0))
    1   z1=z-alegfid(n,0,z,0)/alegfid(n,1,z,0)*dsqrt((1.d0+z)*(1.d0-z))
        if (dabs(z1-z).gt.eps) then
        z=z1
        goto 1
        end if
        x(i) = z1
    2   w(i) = 2.d0/(alegfid(n,1,z1,0)**2)
c
c        at this point, x and w are the gauss points and weights
c        for the interval (-1,+1); x is descending.
c
c        now, re-scale gauss points and weights for interval (ax,bx)
c        such that x is ascending.
c
      alpha=0.5d0*(ax+bx)
      beta=0.5d0*(bx-ax)
      do 3 j=1,n
      x(j)=alpha-beta*x(j)
    3 w(j)=beta*w(j)
c
        return
        end
cern	  c315	    version    05/03/68 alegfid	     94 	       c
      double precision function alegfid(l,m,y,norm)
c     norm = 0 , unnormalised legendre functions
c     norm = 1 , normalised legendre functions
      implicit double precision (a-h,o-z)
      ma = iabs(m)
      if (l - ma) 1,2,2
    1 alegfid = 0.d0
      return
    2 if (m) 3,4,4
    3 ib = l - ma + 1
      ie = l + ma
      prod = 1.d0
      do 5 i = ib,ie
      a = i
    5 prod = prod*a
      fact = 1.0d0/prod
      go to 6
    4 fact = 1.d0
    6 mm = ma
      if (ma - 1) 7,7,8
    7 p0 = 1.0d0
      if (ma) 9,9,10
    9 p1 = y
      go to 11
   10 p0 = 0.0d0
      p1 = dsqrt(1.00 - y**2)
   11 if (l - 1) 12,13,14
   12 alegfid = p0
      go to 23
   13 alegfid = p1 * fact
      go to 23
   14 pnmn1 = p0
      pn = p1
      am = ma
      k = l - 1
      do 15 n = 1,k
      an = n
      pnpl1 = (1.d0/(an-am+1.d0))*
     1        ((2.d0*an+1.d0)*y*pn-(an+am)*pnmn1)
      pnmn1 = pn
   15 pn = pnpl1
   16 alegfid = pnpl1*fact
      go to 23
    8 z2 = (1.d0-y**2)
      am = ma
      ham = am/2.d0
      zham = z2**ham
      if (ma.eq.(l-1)) go to 17
      nb = l+1
      ne = 2*l
      prod = 1.d0
      do 18 ni = nb,ne
      ai = ni
   18 prod = prod*ai
      denom = 2.d0**l
      dnn = prod/denom
      if (ma.eq.l) go to 19
   17 ne = 2*l-1
      prod = 1.d0
      do 20 ni = l,ne
      ai = ni
   20 prod = prod*ai
      denom = 2.d0**(l-1)
      dnm1n = y*prod/denom
      if (ma.eq.(l-1)) go to 21
      me = l-1-ma
      do 22 mn = 1,me
      an = l
      am = l-1-mn
      dnm2n=(1.d0/((an-am)*(an+am+1.d0)))*
     *      (2.d0*(am+1.d0)*y*dnm1n-z2*dnn)
      dnn = dnm1n
   22 dnm1n = dnm2n
      alegfid = dnm2n*zham*fact
      go to 23
   19 alegfid = dnn*zham*fact
      go to 23
   21 alegfid = dnm1n*zham*fact
   23 if (norm.eq.1) go to 24
      return
   24 b = l
      if (m) 25,26,25
   25 ms = -m/mm
      jb = l - mm + 1
      je = l + mm
      prod = 1.d0
      do 27 i = jb,je
      a = i
   27 prod = prod*a
      factor = 0.5d0*(2.d0*b+1.d0)*(prod**ms)
      go to 28
   26 factor = 0.5d0*(2.d0*b+1.d0)
   28 factor = dsqrt(factor)
      alegfid = alegfid*factor
      return
      end
	 function eix(z)
c
c        eix function
c        note: this version of eix is only valid for negative arguments.
c
	 implicit real*8 (a-h,o-z)
         dimension x(100),w(100)
         data pih/1.570796326794897d0/
         logical index/.false./
         save
c
c
         if (z.ge.0.d0) then
         write (6,20000)
20000    format (' error: function eix has a non-negative argument.'/
     1   ' execution terminated.')
         stop
         end if
c
c
         if (index) go to 10
         index=.true.
c
         n=50
         c=1.d0
c****    write (6,10000) n,c
10000    format (i6,f10.5)
c
         call gaussid(0.d0,1.d0,n,x,w)
c
         do i=1,n
         xx=pih*x(i)
         x(i)=c*dtan(xx)
         dc=1.d0/dcos(xx)
         w(i)=c*pih*dc*dc*w(i)
         end do
c
c
   10    zz=-z
c
         aint=0.d0
         do i=1,n
         xx=x(i)+zz
         term=dexp(-xx)/xx*w(i)
         aint=aint+term
         end do
c
         eix=-aint
c
         return
         end
c
c
	 function aibarm1(z)
	 implicit real*8 (a-h,o-z)
         dimension x(1000),w(1000)
         data pih/1.570796326794897d0/
         logical index/.false./
         save
c
c
         if (index) go to 10
         index=.true.
c
         n=200
         c=0.01d0
c****    write (6,10000) n,c
10000    format (i6,f10.5)
c
         call gaussid(0.d0,1.d0,n,x,w)
c
         do i=1,n
         xx=pih*x(i)
         x(i)=c*dtan(xx)+1.d0
         dc=1.d0/dcos(xx)
         w(i)=c*pih*dc*dc*w(i)
         end do
c
c
   10    aint=0.d0
         do i=1,n
         xx=x(i)
         term=dexp(-z*xx)/xx*dlog((xx+1.d0)/(xx-1.d0))*w(i)
         aint=aint+term
         end do
c
         aibarm1=aint
c
         return
         end
c
c
	 function aitilm1(z)
	 implicit real*8 (a-h,o-z)
         dimension x(1000),w(1000)
         data pih/1.570796326794897d0/
         logical index/.false./
         save
c
c
         if (index) go to 10
         index=.true.
c
         n=200
         c=1.d0
c****    write (6,10000) n,c
10000    format (i6,f10.5)
c
         call gaussid(0.d0,1.d0,n,x,w)
c
         do i=1,n
         xx=pih*x(i)
         x(i)=c*dtan(xx)+1.d0
         dc=1.d0/dcos(xx)
         w(i)=c*pih*dc*dc*w(i)
         end do
c
c
   10    aint=0.d0
         do i=1,n
         xx=x(i)
         term=dexp(-z*xx)/xx*dlog(xx+dsqrt(xx*xx-1.d0))*w(i)
         aint=aint+term
         end do
c
         aitilm1=aint
c
         return
         end

      end module