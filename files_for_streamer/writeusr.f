      subroutine writeusr(USER_UNIT,version,avedirc,aveflxu,aveflxd,
     &                avedfdt,swflxd,netflx,dtdt,
     &                cldfrac,surffrac,nstrshort,nstrlong,NCOEF,
     &                surfalb,emis,gemis,nsurfs,whichsurfs,
     &                nover,novertype,overtype,overfrac,
     &                ntype,nlev,zen,ntypetotal,nclouds,
     &                ctopnum,cbottnum,zcthick,USEOLD,
     &                alt,press,temps,wv,ozone,rh,aerosol, 
     &                FRWV,FRO3,FRRHHZ,FRCO2,FRO2,FRWV2,
     &                year,month,day,hour,lat,lon,
     &                STDPROF,AERMOD,AERVERT,cldtauband,cldtemp,
     &                cldpress,cldtau,cldthick,tsurf,rsurf,h2ocol,o3,
     &                tauhaze,albtotal,cldfrc,nrec,satalb,sattb,solflux,
     &                rectitle,infile,IZD,ITD,IWV,IO3,ICTOP,ICTHK,IWAVE,
     &                DESCRIP,descripfile,WEIGHTBND,weightfile,GASABS,
     &                ALBTYPE,EMISSTYPE,CLDFORCE,nphases,cldphase1,
     &                cldre1,cldwc1,cldphase2,cldre2,cldwc2,cldwp,
     &                pcthick,OUTLEVS,INFRD,SOLAR,RAYLISHRT,FLUXES,
     &                wvstart,wvend,bstart,bend,channel,ntheta,theta,
     &                nphi,phi,averad,alb06_adj,bbalb_adj,bandalb,
     &                USERCLOUD,cloudfile,brdffile)
c-----------------------------------------------------------------------------
c  Writes whatever you want!  You have to put in the WRITE statments for
c  whatever variables and in whatever format you desire.  Output must go
c  to the unit USER_UNIT, which corresponds to the file that you (should have)
c  specified in the options ('userfile').  The variables that are available
c  are described briefly below.
c
c  The WRITE statements must take the following form:
c
c       write(USER_UNIT,<fmnt>) <variable-list>
c  
c  Output variable descriptions are given in the file OUTVARS.DOC.
c
c  PARAMETERS (CONSTANTS) USED IN DIMENSIONING ARRAYS: see 'prmstmnt.inc'
c
c  Called by:
c    streamer
c
c  Calls:
c    specint
c-----------------------------------------------------------------------------
      implicit none

      include 'prmstmnt.inc'

      integer    i,j,ntypetotal,nclouds,name
      integer    cbottnum(MAXTYPES-1),ctopnum(MAXTYPES-1)
      integer    nlev,ntype(MAXTYPES),ALBTYPE,EMISSTYPE
      integer    nover,novertype(MCLD),overtype(MCLD,MCLD)
      integer    nphases,cldphase1(MCLD),cldphase2(MCLD)
      integer    OUTLEVS,AERVERT,channel
      integer    month,year,day,nrec,bstart,bend,NCOEF
      integer    ntheta,nphi,nstrshort,nstrlong,STDPROF,AERMOD
      integer    nsurfs,whichsurfs(MAXSURFS)
      integer    IZD,ITD,IWV,IO3,ICTOP,ICTHK,IWAVE,USER_UNIT
      real*8     alt(MLEV),press(MLEV),temps(MLEV),wv(MLEV)
      real*8     ozone(MLEV),rh(MLEV),aerosol(MLEV)
      real*8     h2ocol,o3,tauhaze
      real*8     cldre1(MCLD),cldwc1(MCLD)
      real*8     cldre2(MCLD),cldwc2(MCLD),cldwp(MCLD)
      real*8     pcthick(MAXTYPES-1),zcthick(MAXTYPES-1)
      real*8     cldfrac(MAXTYPES),overfrac(MCLD),cldthick(MAXTYPES)
      real*8     cldtemp(MAXTYPES),cldpress(MAXTYPES),cldtau(MAXTYPES)
      real*8     avedirc(MLEV,2), aveflxu(MLEV,2), dtdt(MLEV)
      real*8     avedfdt(MLEV,2), swflxd(MLEV), netflx(MLEV)
      real*8     aveflxd(MLEV,2),albtotal(MAXTYPES+1)
      real*8     lat,lon,rsurf,tsurf,wvst,wven,wvstart,wvend
      real*8     surffrac(MAXSURFS),zen,raddeg,bandalb
      real*8     cldfrc(MLEV,2),hour,alb06_adj,bbalb_adj
      real*8     cldtauband(NSPEC_TOT,MCLD)
      real*8     surfalb(NSPEC_IR+1:NSPEC_TOT),emis,gemis(NSPEC_IR+1)
      real*8     FRWV,FRO3,FRRHHZ,FRCO2,FRO2,FRWV2,solflux
      real*8     phi(MPHI), theta(MANG), averad(MANG,MPHI,MLEV)
      real*8     dtheta(MANG),satalb(MANG,MPHI),sattb(MANG,MPHI)
      logical    INFRD,SOLAR,FLUXES,GASABS,USERCLOUD
      logical    CLDFORCE,WEIGHTBND,USEOLD,DESCRIP,RAYLISHRT
      character      phasetxt(2)*3, proflabel(NSTANDS)*30
      character      aerlabel(NAERMODS)*30
      character*60   infile,weightfile,version,descripfile
      character*80   rectitle,brdffile,cloudfile

      integer  nbzen,nbtheta,nbphi,nbwave,brdftype
      real*8   bzens(MBZEN),bthetas(MBTHETA),bphis(MBPHI)
      real*8   bwaves(NSPEC_TOT),brdfs(MBPHI,MBTHETA,MBZEN,NSPEC_TOT)
      real*8   brdfparams(MBPARAMS),brdfscaler
      
      common /brdfdat/ brdfparams,bzens,bthetas,bphis,bwaves,brdfs,
     &                 brdfscaler,nbzen,nbtheta,nbphi,nbwave,brdftype

      data raddeg /57.29578/
      data phasetxt /'Liq', 'Ice'/
      data proflabel /'Tropical','Mid-latitude Summer',
     &                'Mid-latitude Winter',
     &                'Subarctic Summer', 'Subarctic Winter',
     &                'Arctic Summer', 'Arctic Winter'/
      data aerlabel /'Tropospheric','Rural','Urban','Maritime','Arctic',
     &               'Smoke'/
c------------------------------------------------------------------------
c  What was the spectral interval for calculations?  'wvst' is the lowest
c  wave number; 'wven' is the highest.  

      call specint(3,wvst,wven,bstart,bend)

c  Convert polar angles to degrees if radiances were computed.  If the cosine
c  was negative then make sure the degrees are negative as well.

      if (.not. FLUXES) then
         do i=1,ntheta
            dtheta(i) = acos(theta(i)) * raddeg
            if (dtheta(i) .gt. 90.) dtheta(i) = dtheta(i) - 180.
         end do
      end if

c************************* USER WRITING STARTS HERE **************************

c  The following code is an example.  Replace it with your own WRITEs.
c23456 marker for line continuation (f77)
      name = ntype(1)
      write(USER_UNIT,*) rectitle
      write(USER_UNIT,10) 'z','p','T','wv','O3','RH','Aer.',
     $ 'LW_up','LW_dn'
10    format ('',a7,a8,a7,a7,a7,a7,a7,a7,a7)

	  do i=1,MLEV
      write(USER_UNIT,11)  alt(i),press(i),temps(i),wv(i),ozone(i), 
     $ rh(i),aerosol(i),aveflxu(i,2),aveflxd(i,2)
11    format ('',f7.2,f8.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2,f7.2)
      enddo

c       write(USER_UNIT,'(a)') 'This message is from "writeusr.f".  '//
c     &   'Did you mean to use CPRINT?  Is it in the right place?'

      return      
      end
