; radiative losses from fit to Chitanti

function lam_rlf, t

return, 0.0*t
END

pro set_zero_rad_loss, quiet=quiet

if( not keyword_set( quiet ) ) then print, '; NO RADIATIVE LOSSES'

return
end

