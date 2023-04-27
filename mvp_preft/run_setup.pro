preft,/verbos                 
forward_function lam_rlf
forward_function rest_force
forward_function tube_erg
common tube_length, tota,chro,coro
bnot=141.0
rec=45.00
tminr=0.02
tminl=0.02
tmax=2.0
n=200
chrl=3.80
chrr=3.80
len=7.0
cor_lenl=10.0
cor_lenr=10.0
ztop=20.0*sin(rec/!radeg)
zbot=0.0
pscale=1.0

flat=0
post=0

stage1dt=1.0
stage1tt=10.0

init_lam_rlf, tmin=1.1d6*tminl
;set_zero_rad_loss
init_gscs,rec,bnot,ztop,zbot

tube=init_rtv_field_var(n=n,chrl=chrl,tminl=tminl,tmax=tmax,len=len,cor_lenl=cor_lenl,cor_lenr=cor_lenr,post=post,ypos=ypos,chrr=chrr,tminr=tminr,pscale=pscale,flat=flat)

field_at_points,tube
drag_const=replicate(0.0,tube.n)
tube.drag_const=drag_const

set_heat_for_equilib, tube
;tube.heat=0.0
tube.visc_max=200.0
tube.inv_hflf=6.
tarr=[tube]
dt=stage1dt
ttot=stage1tt
nst=ceil(ttot/dt)
steptime=0.0
elapsedtime=0.0
stop
print,'Starting run at t=0.0 Tmax=',max(tube.t)
for i=1, nst do begin
   ;cgplot,tube.l,tube.heat,/ylog,yrange=[1d-10,1d3]
   adv_tube, tube, dt, max=10000000L
   tarr = [ tarr, tube ]
   print,'step#:',i,' time=',tube.time,'     Max Temp:',max(tube.t)
   save,tarr,filename='run.sav'
endfor
end
