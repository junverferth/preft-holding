;=========================
function rest_force,tube,dv
common cs ,rec_ang,ztop,zbot,b0
common constr, chrn,dif_le,lengt,peak,cp_rho
if ~keyword_set(cp_rho) then cp_rho=10.
;damp/(2*sprint cnst)>1.0 gives overdamp =1.0 crit damp <1.0 underdamp
spr_cnst = 2.0*((b0/15.0)+cp_rho/15);give some density bit too
damp = 20.0*((b0/5.0)+cp_rho/5.0);give some density bit too
a=[(where(tube.x[2,*] lt zbot,count1))];,(where(tube.x[2,*] gt ztop,count2))]
dv_rest=dblarr(3,tube.n)
vector=dblarr(3,tube.n)

;its here that we need a smarter determination of the vector
;we could do a full cross product, we would need to do it by element
;to make sure that we got the upward facing one. we could just do
;both crossproducts, and then run a where the z component is <0 to
;replace from one to the other?
;if we do it that way, do we need to ensure that we dont have weird
;boundary from flip flopping? shouldnt be two bad though.
vector[0,*]=tube.tv_e[2,*]
vector[1,*]=0.0
vector[2,*]=-tube.tv_e[0,*]

;basic code for swapping from the one plane perp to the other, if the
;direction is bad

alpha=where(vector[2,*] gt 0)
;if alpha[0] ge 0 then begin
   vector[0,alpha]=-tube.tv_e[2,alpha]
   vector[2,alpha]=tube.tv_e[0,alpha]
;endif

;count=[count1,count2]
if (count1 ne 0) then begin;or (count2 ne 0) then begin
   dv_rest[*,a]=vector[*,a]*([1,1,1]#( spr_cnst^2*(tube.x[2,a]-zbot) - damp*total(vector[*,a]*tube.v[*,a],1))) ;for a dense loop, it snags on itself in the x direction
   ;dv_rest[*,*]=0.0;testing for the b-ramp
;temp fix?

endif
;stop
return, dv_rest
end
; -----------------------------
pro field_at_points, tube,dv=dv,chrom=chrom,corleft=corleft,corright=corright,corshape=corshape,bz=bz,bx=bx,flip=flip,by=by;,chrom=chrom,by=by
;  the field strength and its gradients at n points
  common cs,rec_ang,ztop,zbot,b0,depthl,depthr,bpeakl,bpeakr,x0l,x0r,bminl,bminr,nozl,nozr,cor_rise,cdl,cdr
  common fieldshape,fieldlen,fieldstr,fieldderiv,cs_scale,ysc


  if ~keyword_set(corshape) then corshape = 0
  if ~keyword_set(corleft) then corleft=10.0
  if ~keyword_set(corright) then corright=10.0
  B_y=b0*cos(rec_ang/!radeg)
  B_z=b0*sin(rec_ang/!radeg)
  
  bz=sqrt(B_y^2)
  tube.db[0,*] = 0.0
  tube.db[1,*] = 0.0

;cdr cdl are the phys size of the legs
p=0.7
GF=B_y
y=tube.x[2,*]/ysc+p;sqrt(1+2*p*p);+p
by=cs_scale*b0*sin(rec_ang/!radeg)*float(2*sqrt( y^2 - p^2 )/( y^2 + 1.0 )/sqrt( 1.0 + p^2 ))
tube.b=sqrt(GF^2+by^2)
tube.db[2,*]=by^2*y*(1+2*p^2-y^2)/((y^2-p^2)*(1+y^2)*tube.b)/ysc

;stop
c=where(y le p,count)
if count ne 0 then begin
   tube.db[2,c] = 0.0
   tube.b[c] = GF
endif

return
end
; -----------------------------
pro init_fitfield,rec_ang_,ztop_,zbot_,b0_,depthl_,depthr_,bpeakl_,bpeakr_,x0l_,x0r_,bminr_,bminl_,nozl_,nozr_,cor_rise_,cdl_,cdr_,ysc_,ypos=ypos
  common cs, rec_ang,ztop,zbot,b0,depthl,depthr,bpeakl,bpeakr,x0l,x0r,Bminl,Bminr,nozl,nozr,cor_rise,cdl,cdr
common fieldshape,fieldlen,fieldstr,fieldderiv,cs_scale,ysc
nx = 1001
ny = 1001
p = 0.7
ysc = ysc_
xmax = 2.0
ymax = ypos+p
xmin = 0.0
ymin = zbot_ + 0.01
x = xmin + (xmax-xmin)*findgen(nx)/float(nx-1)
y = ymin + (ymax-ymin)*findgen(ny)/float(ny-1)
xx = x # replicate( 1.0, ny )
yy = replicate( 1.0, nx ) # y
ii = complex( 0.0, 1.0 )
zz = xx + ii*yy
a = -alog( ( sqrt( zz^2 + p^2 ) + zz*sqrt(1+p^2) )/( p*sqrt(1-zz^2 ) ) )
a = a + alog( ( zz + sqrt( zz^2 + p^2 ) )/p )/sqrt(1+p^2)
f = float( 2*ii*a )
fcrit = !pi*( 1.0 - 1.0/sqrt(p^2+1) )
; the separatrix
contour, f[1:*,*], x[1:*], y, path_info=path_info, path_xy=path_xy, $
         /path_data, lev=[ fcrit ]

jj = where( ( path_xy[1,*] lt p ) and ( path_xy[0,*] lt 1.0 ) )
x_sx = reform( path_xy[0,jj] )
y_sx = reform( path_xy[1,jj] )
z_sx = x_sx + ii*y_sx
bxy = 2*ii*sqrt( z_sx^2 + p^2 )/( z_sx^2 -1 )/sqrt(1+ p^2 )
bx = imaginary( bxy )
by = float( bxy )
b = sqrt( bx^2 + by^2 )
nb = n_elements( b )
len = fltarr( nb )
dl=fltarr(nb)
for i=1, nb-1 do begin
  dl[i] = sqrt( ( x_sx[i]-x_sx[i-1] )^2 + ( y_sx[i]-y_sx[i-1] )^2 )
  len[i] = len[i-1] + dl[i]
endfor
nft = 15
ft = lsfit( len[(nb-nft):(nb-2)], b[(nb-nft):(nb-2)]^2 )
l0 = -ft[0]/ft[1]
;len = l0 - len

b_cs=b0_*sin(rec_ang_/!radeg)*float(2*sqrt( (y)^2 - p^2 )/( (y)^2 + 1.0 )/sqrt( 1.0 + p^2 ))
bmat=b0_*sin(rec_ang_/!radeg)*float(2*sqrt( (ypos/ysc+p)^2 - p^2 )/( (ypos/ysc+p)^2 + 1.0 )/sqrt( 1.0 + p^2 ))
;b0_=
gf=b0_*cos(rec_ang_/!radeg)
dy=b_cs/sqrt(b_cs^2+gf^2)
dx=b/sqrt(b^2+gf^2)
dz=gf/sqrt(b^2+gf^2)
dz_cs=gf/sqrt(b_cs^2+gf^2)

;bmatch=


db=b-shift(b,1)
db[0]=db[1]
dl_b=len-shift(len,1)
dl_b[0]=dl_b[1]

cs_scale=b0_*sin(rec_ang_/!radeg)/bmat
fieldlen=len
fieldstr=b
fieldderiv=db/dl_b

stop
rec_ang=rec_ang_
ztop=ztop_
zbot=zbot_
b0=b0_
depthl=depthl_


cs_scale=b0_*sin(rec_ang_/!radeg)/bmat
fieldlen=len
fieldstr=b
fieldderiv=db/dl_b

stop
rec_ang=rec_ang_
ztop=ztop_
zbot=zbot_
b0=b0_
depthl=depthl_
depthr=depthr_
bpeakl=bpeakl_
bpeakr=bpeakr_
x0l=x0l_
x0r=x0r_
bminl=bminl_
bminr=bminr_
nozl=nozl_
nozr=nozr_
cor_rise=cor_rise_
cdr=cdr_
cdl=cdl_
return
end
