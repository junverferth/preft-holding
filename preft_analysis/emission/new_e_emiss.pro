function new_e_emiss,tarr,beamparams,lowerbound,upperbound,resolution,phi=phi,order=order,energies=energies,times=times,erg=erg,goes=goes,multi=multi
;this function will take in the tube array, and beam params and
;calculate the proper nonthermal brem from the beam, then interpolate
;it to the stac variables and output that in 3 bands, 3-6 6-8 8-10 kev
  if ~keyword_set(phi) then phi=1d17 ;default to 1d18 Mx of flux
  ;upperbound=100.0;50.0
  ;lowerbound=0.0;0.1
                                ;resolution=1.0                ;0.1
  ;stop
  resolution2=resolution
  energies=findgen((upperbound-lowerbound)/resolution,increment=resolution2,start=lowerbound);kev
  
  energy=energies
  ;resolution[-1]=0.0
  if ~keyword_set(order) then order=100
  tstart=(where(min(abs(tarr.time-min(beamparams.del))) eq abs(tarr.time-min(beamparams.del)))-5>0)[0]
  tstop=(where(min(abs(tarr.time-max(beamparams.dure+beamparams.del))) eq abs(tarr.time-max(beamparams.dure+beamparams.del)))+5)[0]<(n_elements(tarr)-1)
  t=n_elements(tarr)
                                ;stop
  n=n_elements(tarr[0].t)       ;tarr[0].n
  len=fltarr(n_elements(tarr))
  len=tarr.l[-1]
  maxflux=beamparams.amp         ;beamflux in erg/cm2/s
  maxfluxkev=maxflux/1.602d-9    ;flux in kev/cm^2/s
  
  gauss_laguerre_quadr,p,weight,order,0
  
  time=tarr.time                ;time array
  ;energy=dindgen(energies)*resolution+lowerbound
  cutoff=double(beamparams.cutoff)
  dura=beamparams.dure;/multi
  if beamparams.shape eq 1 then amp=(1.0-2*abs(tarr.time-0.5*dura-beamparams.del)/dura)>0.0
  if beamparams.shape eq 0 then amp=fltarr(t)+1.0/beamparams.dure
  if keyword_set(multi) then begin
     amp[21:40]=amp[1:20]
     amp[41:60]=amp[1:20]
     amp[61:80]=amp[1:20]
  endif
  ;stop
  thickbr=dblarr(t,n,n_elements(energy))
  totbr=dblarr(t,n)
  thickbr2=dblarr(t,n,n_elements(energy))
  totbr2=dblarr(t,n)
  ionfrac=fltarr(n_elements(tarr),n_elements(tarr[0].t))
  ;if ion eq 1 then begin
     for i=0,n_elements(tarr)-1 do begin
                                ;stop
        ionfrac[i,*]=get_ieq(tarr[i].t*1d6,'h_2')
     endfor
  ;endif else ionfrac=ionfrac+1.0
  
  m_e=5.10998950d2
  delta=beamparams.index
  avg_e=cutoff*(delta-1.)/(delta-2.)
  l2=25.1+alog(avg_e)
                                ;stop
                                ;time loop
                                ;stop
  for k=tstart,tstop do begin 
     flux=maxfluxkev*amp[k]
     eden=tarr[k].rho*tarr[k].epamu/1.67d-8
     bc=(shift(tarr[k].b,-1)+tarr[k].b)*0.5
     area=phi/bc
     dl=shift(tarr[k].l,-1)-tarr[k].l
     dl[-1]=dl[-2]
                                ;ionfrac=
     l1=66.+1.5*alog(avg_e)-.5*alog(eden)
                                ;cell loop
     for j=0,n-1 do begin
        
                                ;loop over photon_e
        for i=0,n_elements(energy)-1 do begin
           finecon=-0.045848400982856752
           if energy[i] lt cutoff then x=p+cutoff else x=p+energy[i]
           cross_section2=1.1289762657670394d-21 /(energy[i]*x*(x+2*m_e)) $
                          *alog((1.0+sqrt(((x-energy[i])*(x-energy[i]+2*m_e))/(x*(x+2*m_e)))) $
                                /(1.0-sqrt(((x-energy[i])*(x-energy[i]+2*m_e))/(x*(x+2*m_e))))) $
                          *sqrt(x*(x+2*m_e))*(x-energy[i]+m_e) $
                          /(sqrt((x-energy[i])*(x-energy[i]+2*m_e))*(x+m_e)) $
                          *(1.0-exp(finecon*(x+m_e)/sqrt(x*(x+2*m_e)))) $
                          /(1.0-exp(finecon*(x-energy[i]+m_e)/sqrt((x-energy[i])*(x-energy[i]+2*m_e))))
                                ; 1.1289762657670394d-21 = 16 <Z^2> r_0^2 alpha m_e^2 / 3
                                ; 1021.9979 = 2 m_e
                                ; 0.045848400982856752 = 2 pi alpha
                                ; alpha = fine structure constant  
                                ;integral over the xsection
           total_x_sec2=total(weight*x^(2.-delta)*exp(p)*cross_section2)
           const=2.7294774653513575d-9
           thickbr2[k,j,i]= area[j]*flux*(delta-2.0)*total_x_sec2/(cutoff^(2.-delta)*(delta-1.)*(ionfrac[k,j]*l1[j]+(1-ionfrac[k,j])*l2))*const*dl[j]/len[k]
           ;; 2.7294774653513575d-9 = 1 / (8 pi^2 AU^2 e^4 (1.602d-9)^2)
           ;; 1.602d-9 erg to keV conversion
                                ;stop       
        endfor
                                ;stop
     endfor
                                ;stop
  endfor
                                ;stop
  ;three=(where(min(abs(energy-3.)) eq abs(energy-3.)))[0]
  ;six=(where(min(abs(energy-6.)) eq abs(energy-6.)))[0]
  ;eight=(where(min(abs(energy-8.)) eq abs(energy-8.)))[0]
  ;ten=(where(min(abs(energy-10.)) eq abs(energy-10.)))[0]
  if keyword_set(erg) then begin
     for i=0,n_elements(energy)-1 do begin
        phot_erg=1.6022e-9*energy[i]
        thickbr2[*,*,i]=thickbr2[*,*,i]*phot_erg
     endfor
  endif
  ;gl1=(where(min(abs(energy-24.8)) eq abs(energy-24.8)))[0]
  ;gl2=(where(min(abs(energy-3.1)) eq abs(energy-3.1)))[0]
  ;gs1=(where(min(abs(energy-12.4)) eq abs(energy-12.4)))[0]
  ;gs2=(where(min(abs(energy-1.55)) eq abs(energy-1.55)))[0]
  ;stop
  if keyword_set(times) then begin
     nontherm=total(thickbr2,2,/nan)
     if keyword_set(goes) then nontherm=total(nontherm*resolution,2,/nan)
  endif else begin
     nontherm=total(thickbr2*resolution,3,/nan)
  endelse
  ;totbr2=total(thickbr2*resolution,3,/nan)
  ;br236=total(thickbr2[*,*,three:six-1]*resolution,3,/nan)
  ;br268=total(thickbr2[*,*,six:eight-1]*resolution,3,/nan)
  ;br2810=total(thickbr2[*,*,eight:ten-1]*resolution,3,/nan)
  ;goeslong=total(thickbr2[*,*,gl2:gl1-1]*resolution,3,/nan)
  ;goesshort=total(thickbr2[*,*,gs2:gs1-1]*resolution,3,/nan)
                                ;stop
  ;nontherm={units:"photons/cm^2/s(/kev in brem, int over kev in the others)",energies:energy,brem:thickbr2,totalemiss:totbr2,b36:br236,b68:br268,b810:br2810,phi:phi,goeslong:goeslong,goesshort:goesshort}
  ;stop
  return,nontherm
end
