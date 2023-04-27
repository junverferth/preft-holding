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
pro field_at_points, tube,bz=bz,bx=bx,by=by,chrom=chrom,corshape=corshape,flip=flip
;  the field strength and its gradients at n points 
  common cs,rec_ang,ztop,zbot,b0
  common fieldshape,fieldlen,fieldstr,fieldderiv,cs_scale
By=b0*cos(rec_ang/!radeg)
Bz=b0*sin(rec_ang/!radeg)
ramp=1.24
offset=200.0
tube.db[0,*] = 0.0
tube.db[1,*] = 0.0

z=ztop-zbot

;we need to consider how to keep the field from moving upward, in
;2018/2021 we have a steep exponential hump to prevent this, is there
;a better approach. the field we have can be approximated as a linear
;drop with R. so the final result should have to drop too. this
;implies that the field should drop uniformly over x, and then also
;have a drop in z. this may not work. backup plan is a uniform field
;in x and then a linear drop in z, with a hill or well to block lengthening

;stop
edge=where(tube.x[2,*] le zbot or tube.x[2,*] gt ztop)
edgea=where(tube.x[2,*] gt ztop)
edgeb=where(tube.x[2,*] lt zbot)
tube.db[2,edge]=0.0
;stop
return
end
; -----------------------------
pro init_arms,rec_ang_,b0_,ztop_,zbot_
  common cs, rec_ang,ztop,zbot,b0
common fieldshape,fieldlen,fieldstr,fieldderiv,cs_scale

;len = l0 - len

rec_ang=rec_ang_

b0=b0_
zbot=zbot_
ztop=ztop_
return
end
;-----------------------------------------------
