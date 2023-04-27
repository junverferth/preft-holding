
function init_rtv_field_var, n=n,chrr=chrr,chrl=chrl,tminr=tminr,tmax=tmax, len=len,cor_lenl=cor_lenl,tangent=tangent,epsilon=epsilon,ypos=ypos,post=post,offset=offset,cor_lenr=cor_lenr,flat=flat,tminl=tminl,chrom=chrom,exp_l=exp_l,exp_r=exp_r,test=test,corshape=corshape,pscale=pscale
common cs ,rec_ang,ztop,zbot,b0;,depthl,depthr,bpeakl,bpeakr,x0l,x0r,bminl,bminr,nozl,nozr,cor_rise,cdl,cdr
if( not keyword_set( chr ) ) then chr = 2.0
if( not keyword_set( tmin ) ) then tmin = 0.03
if( not keyword_set( tmax ) ) then tmax = 2.0
if( not keyword_set( n ) ) then n = 400
if( not keyword_set( epsilon ) ) then epsilon = 0.00001
if( not keyword_set( len ) ) then len = 3.0*(ztop-zbot)
if( not keyword_set( tangent ) ) then tangent = [1.0,0,0]
if( not keyword_set( cor_lenl ) ) then cor_lenl = 0.0
if( not keyword_set( cor_lenr ) ) then cor_lenr = cor_lenl
if( not keyword_set( offset ) ) then offset = 0.0
if n_elements(post) lt 0 then post = 0
forward_function lam_rlf
common constr,chrn,dif_le,lengt,peak,cp_rho
;with this redone version i would use the ypos keyword, instead of
;tweaking ztop and zobt so that you have the amount in the cs that you want
tot_len=cor_lenl+cor_lenr+len
cp_rho=1.0

;chr=chrl
;tmin=tminl
;tube=rtv_tube_constr(tot_len,chr=chr,tmax=tmax,tmin=tmin,n=n)
;stop
tube=rtv_tube_test(tot_len,chrr=chrr,tmax=tmax,tminr=tminr,n=n,tminl=tminl,chrl=chrl,chrom=chrom,exp_l=exp_l,exp_r=exp_r,corshape=corshape,pscale=pscale)
test=tube


n=tube.n
db=fltarr(3,1)
b=0.0
x=fltarr(3,1)
n=1
v=fltarr(3,1)
tv_e=fltarr(3,1)
time=0
ptube={n:n,x:x,v:v,b:b,db:db,tv_e:tv_e,time:time}
if ~keyword_set(flat) then begin
   ;dl=len/float(n)
   calc_tv,tube
   calc_len,tube
   ;find midpoint
   nh=tube.n/2
   cel=max(where(tube.t[0:nh] le tminl))
   cer=min(where(tube.t[nh:(tube.n-1)] le tminr))+nh-1
   xel=tube.x[0,cel]
   xer=tube.x[0,cer]
   mid_point=((xer-xel)/2.0+xel)
   mid_point=min(abs(tube.x[0,*]-mid_point),index)
   mid_point=index
   tube.x[0,*]=tube.x[0,*]-tube.x[0,mid_point]
   tube.x[2,*]=(ztop-zbot)/2.0
   dz=0
   dx=0
   if offset ne 0 then begin
      nooffpos=tube.x[*,mid_point]
      oofpos=min(abs(tube.x[0,*]-offset),index)
      mid_point=index
   endif
    if post eq 0 then begin
       i = mid_point+1
       tangent=[cos(rec_ang/!radeg),0,-sin(rec_ang/!radeg)]
                                 ;rightside integration
; stop
       while ((tube.x[2,i-1]+dz) le ztop) do begin
          ptube.x=tube.x[*,i-1]
          field_at_points,ptube,bz=bz,bx=bx,by=by
          dz=tube.dl_e[i-1]*bz/sqrt(by^2+bz^2);ptube.b
          dx=tube.dl_e[i-1]*by/sqrt(by^2+bz^2);ptube.b
          tangent=[dx/sqrt(dx^2+dz^2),0,dz/sqrt(dx^2+dz^2)]
          tube.x[*,i]=tube.x[*,i-1]+[dx,0,dz]
          ;stop
          i=i+1
       endwhile
; stop      
       n=tube.n ;try to make the last piece have its correct length but still sit up on the edge
       tube.x[2,(i):(n-1)]=ztop
       tube.x[0,i]=sqrt(tube.dl_e[i-1]^2-(tube.x[2,i]-tube.x[2,i-1])^2)+tube.x[0,i-1]
       for k=i,n-1 do begin
          tube.x[0,k]=tube.dl_e[k-1]+tube.x[0,(k-1)]
       endfor
    ;stop                            ;leftside integration
       i=mid_point-1
       tangent=[-cos(rec_ang/!radeg),0,-sin(rec_ang/!radeg)]
       while((tube.x[2,i+1]+dz) ge zbot) do begin
          ptube.x=tube.x[*,i+1]
      ;stop
          field_at_points,ptube,bz=bz,bx=bx,by=by
          dz=-tube.dl_e[i-1]*bz/sqrt(by^2+bz^2);ptube.b
          dx=-tube.dl_e[i-1]*by/sqrt(by^2+bz^2);ptube.b
          tangent=[dx/sqrt(dx^2+dz^2),0,dz/sqrt(dx^2+dz^2)]
          tube.x[*,i]=tube.x[*,i+1]+[dx,0,dz]
          i=i-1
       endwhile
;stop
       k=i
       tube.x[2,0:(i)]=zbot
       tube.x[0,i]=-sqrt(tube.dl_e[i]^2-(tube.x[2,i+1]-tube.x[2,i])^2)+tube.x[0,i+1]

       while (k ge 0) do begin
          tube.x[0,k]=tube.x[0,(k+1)]-tube.dl_e[k]
          k=k-1
       endwhile
;stop
    endif
   
   if post eq 1 then begin
      if(keyword_set(ypos)) then begin
         tube.x[2,*]=ypos
         ;stop
         i = mid_point+1
         tangent=[cos(rec_ang/!radeg),0,-sin(rec_ang/!radeg)]
           ;stop                     ;rightside integration
         while ((tube.x[2,i-1]+dz) ge zbot) do begin
            ptube.x=tube.x[*,i-1]
            field_at_points,ptube,bz=bz,bx=bx,by=by
            dz=-tube.dl_e[i-1]*bz/sqrt(by^2+bz^2);ptube.b
            dx=tube.dl_e[i-1]*by/sqrt(by^2+bz^2);ptube.b
            tangent=[dx/sqrt(dx^2+dz^2),0,dz/sqrt(dx^2+dz^2)]
            tube.x[*,i]=tube.x[*,i-1]+[dx,0,dz]
            ;stop
            i=i+1
         endwhile
         n=tube.n
         tube.x[2,(i):(n-1)]=zbot
         tube.x[0,i]=tube.x[0,i-1]+sqrt(tube.dl_e[i-1]^2-(tube.x[2,i-1]-zbot)^2)
     ;stop
         for k=i,n-1 do begin
            tube.x[0,k]=tube.dl_e[k-1]+tube.x[0,(k-1)]
         endfor
                                ;leftside integration
         for i=0,mid_point do begin
            tube.x[0,mid_point-i]=-tube.x[0,mid_point+i]
            tube.x[2,mid_point-i]=tube.x[2,mid_point+i]
         endfor
      endif
   
      if(~keyword_set(ypos)) then begin
         i = mid_point+1
         tangent=[cos(rec_ang/!radeg),0,-sin(rec_ang/!radeg)]
                                ;rightside integration
         while ((tube.x[2,i-1]+dz) ge zbot) do begin
            ptube.x=tube.x[*,i-1]
            field_at_points,ptube,bz=bz,bx=bx,by=by
            dz=-tube.dl_e[i-1]*bz/sqrt(by^2+bz^2);ptube.b
            dx=tube.dl_e[i-1]*by/sqrt(by^2+bz^2);ptube.b
            tangent=[dx/sqrt(dx^2+dz^2),0,dz/sqrt(dx^2+dz^2)]
            tube.x[*,i]=tube.x[*,i-1]+[dx,0,dz]
            i=i+1
         endwhile
         n=tube.n
                             ;   stop
         tube.x[2,(i):(n-1)]=zbot
         tube.x[0,i]=tube.x[0,i-1]+sqrt(tube.dl_e[i-1]^2-(tube.x[2,i-1]-zbot)^2)
         ;i=i+1
 
         for k=i+1,n-1 do begin
            tube.x[0,k]=tube.dl_e[k-1]+tube.x[0,(k-1)]
         endfor
         ;stop                   ;leftside integration

         for i=0,mid_point do begin
            tube.x[0,mid_point-i]=-tube.x[0,mid_point+i]
            tube.x[2,mid_point-i]=tube.x[2,mid_point+i]
         endfor
         
      endif
      ;stop
      round_rx_bend,tube,mid_point,nb=5
   endif
   peak=where( tube.x[2,*] ge max(tube.x[2,*]))
   peak=peak[0]
endif
;====================================================================================================================================================
if keyword_set(flat) then begin
   tube.x[0,*]=tube.x[0,*]-tube.x[0,((tube.n-1)/2.0)]
   ;stop
   return,tube
endif
;`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-`-
;stop
calc_tv,tube
calc_len,tube
;r=0.01*(1.38/1.67/tube.mpp)
;tube.rho=tube.p/r/tube.t
set_rho,tube
calc_prho,tube
;stop
lengt=max(tube.l)
cp_rho=min(tube.rho)
;save,/variables,filename='init_rtv_save.sav'
;stop
return, tube
end
; --------------------------------------------------------------

