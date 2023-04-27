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
pro field_at_points,tube
;  the field strength and its gradients at n points 
  common cs,rec_ang,ztop,zbot,b0
  common fieldshape,fieldlen,fieldstr,fieldderiv,cs_scale
  by=b0
  bz=0.0
  tanhscale=1.5
  
  tube.b=b0
  tube.db[0,*]=0.0
  tube.db[1,*]=0.0
  tube.db[2,*]=0.0
  mid=where(tube.x[2,*] le cs_scale and tube.x[2,*] ge -cs_scale)
                                ;sin2
  tube.b[mid]=b0+(sin(!pi*tube.x[2,mid]/2.0/cs_scale)*sin(!pi*tube.x[2,mid]/2.0/cs_scale)-1)*fieldstr
  tube.db[2,mid]=sin(!pi*tube.x[2,mid]/2.0/cs_scale)*cos(!pi*tube.x[2,mid]/2.0/cs_scale)*fieldstr*!pi/cs_scale
    ;stop                           ;tanh?
  ;tube.b[mid]=b0+fieldstr/2.0*(tanh((tube.x[2,mid]-cs_scale/tanhscale)*fieldlen)-tanh((tube.x[2,mid]+cs_scale/tanhscale)*fieldlen))
  ;tube.db[2,mid]=fieldstr/2.0*fieldlen*(1.0/(cosh(fieldlen*(tube.x[2,mid]-cs_scale/tanhscale)))^2-1.0/(cosh(fieldlen*(tube.x[2,mid]+cs_scale/tanhscale)))^2)

  ;stop
end
pro init_alfven,bnot,width,height,ztop_,zbot_,tanhwidth
  common cs, rec_ang,ztop,zbot,b0
  common fieldshape,fieldlen,fieldstr,fieldderiv,cs_scale
  b0=bnot
  cs_scale=width
  fieldstr=height
  fieldlen=tanhwidth
  ztop=ztop_
  zbot=zbot_
  
end
