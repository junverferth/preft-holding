; radiative losses from fit to Chitanti

function lam_rlf, t

common rad_loss_arrays, t0, t1, f0, p

lam = t*0d0

nr =  n_elements( t0 )
FOR i=0, nr-1 DO BEGIN
  j = where( ( t GE t0(i) ) AND ( t LE t1(i) ) )
  IF( j(0) GT -1 ) THEN lam(j) = f0(i)*t(j)^(p(i))
ENDFOR

return, lam
END

pro init_new_lam, tmin=tmin, quiet=quiet

common rad_loss_arrays, t0, t1, f0, p
;lt0 = [4.0,4.97,5.67,6.18,6.55,6.90,7.63]
lt0=[4.0,4.159,4.28,4.43,4.6,4.89,5.4,5.58,6.17,6.53,7.0,7.37,7.65,8.49]
;lt1 = [4.97,5.67,6.18,6.55,6.90,7.63,9.00];
lt1=[4.159,4.28,4.43,4.6,4.89,5.4,5.58,6.17,6.53,7.0,7.37,7.65,8.49,9.0]
;lf0 = [-30.9626,-16.0555,-21.7212,-12.4522,-24.4609,-15.2604,-26.7077] ;
lf0=[-62.187008608171269,-22.05021674392341,-13.231078912235443,-22.686566554570273,-29.853769394826024,-22.253305359309092,-13.895763199774283,-22.788435356504870,-8.8522077601966540,-23.119514414288613,-13.750612193652493,-22.864702228073543,-25.402215669649991,-26.560021938841437]
;p = [2.0,-1.0,0.0,-1.5,1./3.,-1.0,0.50] ;
p=[9.7143440753862897,0.075691172282918412,-1.9859310965066082,0.14672919872864534,1.7065680016256219,0.15336378886634547,-1.3948801749272668,0.19792218698292527,-2.0609269737824487,0.12342338626463684,-1.2146239204378744,0.022602225056695568,0.35422035987900935,0.49060123799525213]
;stop
if( keyword_set( tmin ) ) then begin

   if alog10(tmin) gt lt0[1] then begin
      lowend=(where(alog10(tmin) lt lt0))
      lt0=[alog10(tmin),lt0[lowend]]
      lt1=[lt1[lowend[0]-1],lt1[lowend]]
      lf0=[lf0[lowend[0]-1],lf0[lowend]]
      p=[p[lowend[0]-1],p[lowend]]
   endif else lt0[0] = alog10( tmin )

endif

if( not keyword_set( quiet ) ) then begin
   print, '; PPL radiative loss function - version 3'
   print, '; CHIANTI v10.0.2 schmelz_ext'
   print, ';  tmin = ', 10.d0^(lt0[0])
endif
;stop
t0 = 10.0d0^(lt0)
t1 = 10.0d0^(lt1)
f0 = 10.0d0^(lf0)
p=double(p)
;help,t0,t1,f0
;print,t0
;print,t1
;print,f0
return
end

