
;------------------------------------------------------
pro apply_beam,tube;,repeats,shape,dure,del,gap
  common ebeam, amp,index,cutoff,dure,shape,h00,del,eprof,sym,repeats
  ;may wanna get clever with 3 phases, fast ramp, 2nd/slowramp/flat, decay
                                ;if(tube.time lt del) then return
  gap=0.0
  if shape eq 1 then begin
     amplitude=0                ;(1.0-2*abs(tube.time-0.5*dure-del)/dure)>0.0
     ;cgplot,0,0,/nodata,xr=[0,50],yr=[-2,2]
     for  i=1, repeats do begin
        repamp=(1.0-2*abs(tube.time-i*(0.5*dure+del)-(i-1)*(0.5*dure+gap))/dure)>0.0
        ;cgoplot,tube.time,repamp
        amplitude=amplitude+repamp
     endfor
  endif
  ;cgoplot,tube.time,amplitude,color='red'
 ; stop
  if shape eq 0 then amplitude=1.0/dure
  if tube.time gt repeats*(dure+gap+del) then begin
     tube.heat=h00
     return
  endif
  ;stop
  ;amp=
  eheat_calc,tube,eheat,amplitude
  ;stop
  half=tube.n/2
  if sym eq 1 then tube.heat=[h00[0:half]+eheat[0:half],h00[half+1:tube.n-1]] else tube.heat=h00+eheat
 
end
;---------------------,----------------------------------
pro eheat_calc,tube,eheat,amplitude,debug=debug
  common ebeam, amp,index,cutoff,dure,shape,h00,del,eprof,sym,repeats

;----------------- from ionconstants in preft
a_he = 10.0d0^( -1.075 );    solar He abundance = 10.925
mph = 1.0 + 4*a_he;          mass per H nucleus [ amu ]
pph = 2.0 + 3*a_he;          particles per H nucleus
eph = 1.0 + 2*a_he;          electrons per H nucleus

ipamu=(pph-eph)/mph
;-----------------------------

;let us take the midpoint, and calculate column from there to either end
col=fltarr(tube.n)
midn=tube.n/2
eden=tube.rho*tube.epamu/1.67d-8
col=total(eden*tube.dl_e*1d8,/cumulative)
col=abs(col-col[midn])
col[midn]=eden[midn]*5d6;adhoc giving a 50km slab of material at the 0 point to avoid huge spike from betafunc
iden=ipamu/1.67d-8*tube.rho
delta=abs(index)
ece=cutoff*1.6022d-9
avge=ece*(1-delta)/(2-delta)
L=65.1+alog(avge)-0.5*alog(eden)
fc=amp
e=4.803e-10
k=2.*!PI*e^4.
gam=l   
u0=1.
bet=2.
col_st=u0*ece^2/(3.*gam*k)
xc=col/col_st
mask=where(xc lt 1.0)
beta_fun=fltarr(tube.n)+1.0
if xc[0] ne -1 then beta_fun[mask]=Ibeta(delta/2.,1./3.,xc[mask])
beta_fun=beta(delta/2.,1./3.)*beta_fun 
;emslie78
;I(N,t) is ergs/cm^-3/s
;	I(N,t)=1/2 K n gam (delta-2) Beta(delta/2,2/(4+bet)) Flux(t)/E1^2((2 + bet/2)gam K N/(mu0 E1^2))^(-delta/2)
;	K=2 pi fundcharge^4;
;	n=local number density;
;	gam=LAM
;	LAM=ln(mo v^2 min(eta,rl)/Z z e^2)
;	rl is larmor radius
;	eta is velocity/plasma osc freq
;	mo=Mm/(M+m)
;	v=velocity
;	z=electron charge
;	Z=nuclear charge
;	e=esu
;	delta= parameter
;	bet=2
;	flux=energy flux at ergs/cm/cm/s
;	E1=cutoff energy 
;	N=column density, int from -inf to x of n(xp)dxp
;	mu0=1(beamed electrons)
;	I=1*pi*e^4*numberdensity of target*LAM*(delta-2)
;Beta(delta,1/3) Flux /E1^2
;(3*LAM*2*pi*e^4*column/(mu0*1^2))^(-delta/2)
I=1*!pi*e^4*eden*gam*(delta-2.0)*beta_fun*(amp*amplitude)/ece^2*(col/col_st)^(-delta/2)
;;erg/cm^3/s,dont trust using the complete beta function everywhere
;like in emslie 78, so we are using the beta function approach from
;1994 hawley
eheat=i*tube.dl_e ;(yeaaaaa why arent we using 1d8*tube.dl_e here?) Suspicious now about the lack of tube.b to get me into cm from mx...and why my length is in megameters Length is left in Mm as heat is stored as 1d8 erg/cm^2/s, so no mx conversion and leaving in Mm is fine.

eprof=eheat
if keyword_set(debug) then stop
end
;------------------------------------------------------------------------

pro init_ebeammulti, amp_, index_, cutoff_,dure_,shape_,del_,tube,sym_,repeats_
  common ebeam, amp,index,cutoff,dure,shape,h00,del,eprof,sym,repeats
  ;the point of this module is to load the user supplied parameters of the ebeam into the common module
  ;we have amplitude, spectral index,low nergy cutoff, duration, and shape with 0 being flat and 1 triangle
  amp=amp_;currently this is total e flux in beam.
  index=index_
  cutoff=cutoff_
  dure=dure_
  del=del_
  shape=shape_
  h00=tube.heat
  sym=sym_
  repeats=repeats_
                                ;there needs to be normalization here,
                                ;cus thats unregulated compared to the heating
end
