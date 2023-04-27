; radiative losses from fit to Chitanti

function lam_rlf, t

common rad_loss_arrays, t0, t1, f0, p

lam = t*0.0

nr =  n_elements( t0 )
FOR i=0, nr-1 DO BEGIN
  j = where( ( t GE t0(i) ) AND ( t LE t1(i) ) )
  IF( j(0) GT -1 ) THEN lam(j) = f0(i)*t(j)^(p(i))
ENDFOR

return, lam
END

pro init_lam_rlf, tmin=tmin, quiet=quiet

common rad_loss_arrays, t0, t1, f0, p
;new arrays match the HYDRAD power loss RLF
lt0 = [  4.0, 4.896, 5.419, 5.563, 6.183, 6.563, 6.978, 7.467 ];[4.0,4.97,5.67,6.18,6.55,6.90,7.63];
lt1 = [  4.896, 5.419, 5.563, 6.183, 6.563, 6.978, 7.467, 9.000 ];[4.97,5.67,6.18,6.55,6.90,7.63,9.00];
lf0 = [  -29.411, -21.927, -10.565, -22.849,  -8.679, -23.867, -13.248, -25.105 ];[-30.9626,-16.0555,-21.7212,-12.4522,-24.4609,-15.2604,-26.7077];
p = [    1.659,   0.131,  -1.966,   0.242,  -2.050,   0.264,  -1.257,   0.331 ];[2.0,-1.0,0.0,-1.5,1./3.,-1.0,0.50];

if( keyword_set( tmin ) ) then lt0[0] = alog10( tmin )
if( not keyword_set( quiet ) ) then begin
  print, '; PPL radiative loss function - version 2'
  print, ';  tmin = ', 10.d0^(lt0[0])
endif

t0 = 10.0d0^(lt0)
t1 = 10.0d0^(lt1)
f0 = 10.0d0^(lf0)

return
end

