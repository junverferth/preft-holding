; a set of programs which will update the heating profile of a tube
; over time.

pro apply_heater, tube

  common heater_params, h0, hprof, dur_h, hmax,delay,symm
  if(tube.time lt delay) then return
  if( tube.time gt 1.0*dur_h+delay ) then begin
     tube.heat = h0
     return
  endif
  dur=dur_h
  amp = ( 1.0 - 2*abs( tube.time - 0.5*dur-delay )/dur_h ) > 0.0
  ;gauamp= 1*exp(-0.5*((tube.time-0.5*dur_h)/2.021587)^2)
  half=tube.n/2
  if symm eq 1 then tube.heat=[h0[0:half]+amp*hmax*hprof[0:half],h0[half+1:tube.n-1]] else tube.heat = h0 + amp*hmax*hprof ;*tube.b/tube.b[0] 
  ;tube.heat=h0+gauamp*hmax*hprof
  ;print,amp,hmax,max(hprof)
return
end

; ---------------------------------------------

pro init_heater, tube, etot, dur=dur, len=len, l0=l0,delay_=delay_,shape=shape,symm_=symm_
common heater_params, h0, hprof, dur_h, hmax,delay,symm
symm=symm_
h0 = tube.heat;  save the initial heat to be augmented
calc_len, tube
if( not keyword_set( delay_ ) ) then delay_ = 0.0;delay before we kickoff the heating
if( not keyword_set( l0 ) ) then l0 = 0.5*max( tube.l );  midpoint
if( not keyword_set( len ) ) then len = 0.5*max( tube.l );  heat center 50%
if( not keyword_set( dur ) ) then dur = 10.0;  10 seconds of heating
dur_h = dur
delay=delay_
;h2prof = ( 1.0 - 2*abs( tube.l - l0 )/len ) > 0.0 ;  normalized tent
hprof=dblarr(tube.n)
if shape eq 0 then hprof=(1/sqrt(2*!pi)/l0)*exp(-0.5*((tube.l-len)/l0)^2) ;gaussian profile
if shape eq 1 then begin
   low=2.0
   high=5.0
   feet=[where(tube.l le (high) and tube.l ge (low)),where(tube.l ge (tube.l[tube.n-1]-high) and tube.l le (tube.l[tube.n-1]-low))]
   hprof[feet]=1.0                                      ;foot point centric
endif
if shape eq 2 then hprof=replicate(1.0,tube.n) ;uniform heating
if shape eq 3 then hprof=((1/sqrt(2*!pi)/l0)*exp(-0.5*((tube.l-len)/l0)^2))+((1/sqrt(2*!pi)/l0)*exp(-0.5*((tube.l-(tube.l[tube.n-1]-len))/l0)^2));double gaussian focussed in lower leg
;h2prof=h2prof*tube.dl
hprof = hprof*tube.dl;   will store heat in each cell
;stop
htot = total( hprof )
hprof = hprof/htot;  normalized: integrates to 1
hmax = 2*etot/dur ;  for tent time profile

;h2tot=total(h2prof)
;h2prof=h2prof/h2tot
;h2max=2*etot/dur
;stop
return
end



