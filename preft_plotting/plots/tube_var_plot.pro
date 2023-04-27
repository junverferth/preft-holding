pro tube_var_plot, tarr,name, fp=fp, vmax=vmax, vfp=vfp, denmax=denmax, tmax=tmax, $
  lhalf=lhalf,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,IC=IC,FC=FC,thick=thick

if( not keyword_set( fp ) ) then fp=0.6*max(tarr[0].l)
if( not keyword_set( denmax ) ) then denmax=1.0d13
if( not keyword_set( tmax ) ) then tmax=1.0d8
if( ~keyword_set(xmin))then xmin= -3.0
if( ~keyword_set(xmax))then xmax= 50.0
if( ~keyword_set(ymin))then ymin= -1.0
if( ~keyword_set(ymax))then ymax= 10.0
;jups_open,name,900,1000
;calc_len, tarr[0]
ink='black'
background='white'
!p.thick=1
!p.charsize=1.0
!p.charthick=1
lh = max(tarr[0].l)                
if n_elements(fC) gt 0 then fh = max(fc.l)
if n_elements(IC) gt 0 then ih = max(ic.l)
if( not keyword_set( lhalf ) ) then lhalf=lh
;print,'check2'
ypos=[0.07,.315,.325,.6,.65,.95]
xpos=[.09,.5,.5,.91]
;ypos=[0.06,0.25, 0.44, 0.63, 0.81,0.86,0.96 ]
n_e = tarr[0].rho*tarr[0].epamu/1.67d-8
vpar = total( tarr[0].v*tarr[0].tv_e, 1 )
if( not keyword_set( vmax ) ) then vmax = max( vpar )
if( not keyword_set( vfp ) ) then vfp = 1.5
leftedge=tarr[0].l[0]-lh/2.0
mid=lh*.6+leftedge
colors=['black','blue','dark green','orange','red'] ;,'purple'
lines=[0,2,4,0,5]


cgplot, tarr[0].l-lh/2.0, n_e, xr=[leftedge,mid], xst=1, pos=[xpos[0], ypos[2], xpos[1], ypos[3]], ytit='n!de!n [cm!u-3!n]', yr=[1.0d8,denmax],color=ink,background=backgound,/ylog,thick=thick,xtickf='no_tick_label',charsize=0.7 ;1.5;0.7
cgoplot,ic.l-ic.l[-1]/2.0,ic.rho*ic.epamu/1.67d-8,lines=5,color='red'
cgoplot, [0,0], [ 1.0d6, denmax], lines=1

cgplot, tarr[0].l-lh/2.0, 1.0d6*tarr[0].t, xr=[leftedge,mid], xst=1, /noerase, pos=[xpos[2], ypos[0], xpos[3], ypos[1]], yr=[1.0e4,tmax],color=ink,background=backgound,/ylog,xtit='Distance from Loop Center [ Mm ]',thick=thick,ytickf='no_tick_label',charsize=0.7 ;1.5;0.7
cgoplot,ic.l-ic.l[-1]/2.0,ic.t*1d6,lines=5,color='red'
cgoplot, [0,0], [.1,1d10 ], lines=1
cgaxis, /yax, yst=1, ytit='T [K]',charsize=0.7;1.5;0.7

cgplot, tarr[0].l-lh/2.0, vpar, xr=[leftedge,mid], xst=1, /noerase, pos=[xpos[0], ypos[0], xpos[1], ypos[1]], yr=vfp*[-1.0,1.0], ytit='v [Mm/s]',color=ink,background=backgound,thick=thick,xtit='Distance from Loop Center [Mm]',charsize=0.7 ;1.5;0.7
cgoplot,ic.l-ic.l[-1]/2.0,total(ic.v*ic.tv_e,1),lines=5,color='red'
cgoplot, [0,0], [ -4,4 ], lines=1

cgplot, tarr[0].l-lh/2.0, tarr[0].p, xr=[leftedge,mid], xst=1, /noerase, pos=[xpos[2], ypos[2], xpos[3], ypos[3]],ytickf='no_tick_label', yr=[.1,100],color=ink,background=backgound,/ylog,thick=thick,xtickf='no_tick_label',charsize=0.7 ;1.5;0.7
cgoplot,ic.l-ic.l[-1]/2.0,ic.p,lines=5,color='red'
cgoplot, [0,0], [ 0.1,100.0 ], lines=1
cgaxis, /yax, yst=1, ytit='p [erg/cm!u3!n]',charsize=0.7;1.5;0.7
;cgtext,.5,.98,string(tarr.time),/normal,charsize=1.,alignment=0.5
xmin=1.1*tarr[0].x[0,0]

cgplot,tarr[0].x[0,*],tarr[0].x[2,*],yrange=[-10.0,50.],xrange=[-60,60],xtitle='X [ Mm ]',ytitle='Z [ Mm ]',title='',/noerase,pos=[xpos[0],ypos[4],xpos[3],ypos[5]],color=ink,background=backgound,thick=thick,/iso,charsize=0.7 ;1.5;0.7
cgoplot,ic.x[0,*],ic.x[2,*],lines=5,color='red'
;print,lh,leftedge,mid

; for i=1,n_elements(tarr)-1 do begin
;    lh = max(tarr[i].l)
;    ;leftedge=tarr[i].l[0]-lh/2.0
;    ;mid=lh*.6+leftedge
   
;    n_e = tarr[i].rho*tarr[i].epamu/1.67d-8
;    n_e=smooth(n_e,5)
;    vpar = total( tarr[i].v*tarr[i].tv_e, 1 )
   
;    cgplot, tarr[i].l-lh/2.0, n_e, xr=[leftedge,mid], xst=1, pos=[xpos[0], ypos[2], xpos[1], ypos[3]], yr=[1.0d8,denmax],background=backgound,/ylog,lines=lines[i],thick=thick,xtickf='no_tick_label', /noerase,color=colors[i],charsize=0.7;1.5;0.7

;    cgplot, tarr[i].l-lh/2.0, 1.0d6*tarr[i].t, xr=[leftedge,mid], xst=1, /noerase, pos=[xpos[2], ypos[0], xpos[3], ypos[1]], yr=[1.0e4,tmax],background=backgound,/ylog,lines=lines[i],thick=thick,ytickf='no_tick_label',color=colors[i],charsize=0.7;1.5;0.7

;    cgplot, tarr[i].l-lh/2.0, vpar, xr=[leftedge,mid], xst=1, /noerase, pos=[xpos[0], ypos[0], xpos[1], ypos[1]], yr=vfp*[-0.50,1.0],background=backgound,lines=lines[i],thick=thick,color=colors[i],charsize=0.7;1.5;0.7

;    pres=smooth(tarr[i].p,5)
;    cgplot, tarr[i].l-lh/2.0, pres, xr=[leftedge,mid], xst=1, /noerase, pos=[xpos[2], ypos[2], xpos[3], ypos[3]],ytickf='no_tick_label', yr=[.1,100],background=backgound,/ylog,lines=lines[i],thick=thick,xtickf='no_tick_label',color=colors[i],charsize=0.7;1.5;0.7

; print,lh,leftedge,mid
;    cgplot,tarr[i].x[0,*],tarr[i].x[2,*],/noerase,pos=[xpos[0],ypos[4],xpos[3],ypos[5]],color=colors[i],lines=lines[i],thick=thick,/iso,yrange=[-10,50],xrange=[-60,60],charsize=0.7;1.5;0.7
; endfor

;cgtext,.18,.92,'(a)',/normal,charsize=1.;50
;cgtext,.11,.56,'(b)',/normal,charsize=1.;50
;cgtext,.52,.56,'(c)',/normal,charsize=1.;50
;cgtext,.11,.28,'(d)',/normal,charsize=1.;50
;cgtext,.52,.28,'(e)',/normal,charsize=1.;50



;jups_close
;write_png,'flat_timesteps.png',tvrd(/true)
;write_png,'peak_timesteps.png',tvrd(/true)
;write_png,'high_timesteps.png',tvrd(/true)
end
