;=========================
function rest_force,tube,dv
  common cs, rec_ang,ztop,zbot,b0,depthl,depthr,bpeakl,bpeakr,x0l,x0r,Bminl,Bminr,nozl,nozr,cor_rise,cdl,cdr
  common constr, chrn,dif_le,lengt,peak,cp_rho
  if ~keyword_set(cp_rho) then cp_rho=10.
;damp/(2*sprint cnst)>1.0 gives overdamp =1.0 crit damp <1.0 underdamp
  spr_cnst = 2.0*((b0/15.0)+cp_rho/15)  ;give some density bit too
  damp = 20.0*((b0/5.0)+cp_rho/5.0)     ;give some density bit too
  a=[(where(tube.x[2,*] lt zbot,count1))] ;,(where(tube.x[2,*] gt ztop,count2))]
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
  if (count1 ne 0) then begin   ;or (count2 ne 0) then begin
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
  common fieldshape,fieldlen,fieldstr,fieldderiv,cs_scale,ysc,a,p


  if ~keyword_set(corshape) then corshape = 0
  if ~keyword_set(corleft) then corleft=10.0
  if ~keyword_set(corright) then corright=10.0

  B_y=25.;b0*cos(rec_ang/!radeg)


  tube.db[0,*] = 0.0
  tube.db[1,*] = 0.0

  bz=B_y
;y=tube.x[2,*]/ysc+p;sqrt(1+2*p*p);+p
;b_z=cs_scale*b0*sin(rec_ang/!radeg)*float(2*sqrt( y^2 - p^2 )/( y^2 + 1.0 )/sqrt( 1.0 + p^2 ))
;tube.b=sqrt(GF^2+b_z^2)
;tube.db[2,*]=b_z^2*y*(1+2*p^2-y^2)/((y^2-p^2)*(1+y^2)*tube.b)/ysc
  z=tube.x[2,*]+p
  scaled=sqrt(z^2-p^2)/(z^2+a^2)
  fac=2*B0*sqrt(p^2+a^2)
  by=fac*scaled
  b=sqrt(by^2+bz^2)
  dby=fac*(z*(a^2+2*p^2-z^2)/(sqrt(z^2-p^2)*(z^2+a^2)^2))
  db=by*dby/b
  tube.b=b
  tube.db[2,*]=db


  ;increment=5d-4
  ;zp=z+increment
  ;zm=z-increment
  ;scaledp=sqrt(zp^2-p^2)/(zp^2+a^2)
  ;scaledm=sqrt(zm^2-p^2)/(zm^2+a^2)
  ;byp=fac*scaledp
  ;bym=fac*scaledm
  ;bp=sqrt(byp^2+bz^2)
  ;bm=sqrt(bym^2+bz^2)
  ;dbcheck=(bp-bm)/(2*increment)

  ;stop
  c=where(z le p,count)
  if count ne 0 then begin
     tube.db[2,c] = 0.0
     tube.b[c] = bz
  endif
;  stop
  return
end

; -----------------------------
pro init_pinch,rec_ang_,ztop_,zbot_,b0_,depthl_,depthr_,bpeakl_,bpeakr_,x0l_,x0r_,bminr_,bminl_,nozl_,nozr_,cor_rise_,cdl_,cdr_,p_,a_,ypos=ypos
  common cs, rec_ang,ztop,zbot,b0,depthl,depthr,bpeakl,bpeakr,x0l,x0r,Bminl,Bminr,nozl,nozr,cor_rise,cdl,cdr
  common fieldshape,fieldlen,fieldstr,fieldderiv,cs_scale,ysc,a,p
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
  p=p_
  a=a_
  return
end
;-----------------------------------------------
