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
pro field_at_points, tube,bz=bz,bx=bx,by=by,chrom=chrom,corshape=corshape,flip=flip,flat=flat
;  the field strength and its gradients at n points 
  common cs,rec_ang,ztop,zbot,b0
  common fieldshape,fieldlen,fieldstr,fieldderiv,cs_scale
By=b0*cos(rec_ang/!radeg)
Bz=b0*sin(rec_ang/!radeg)
bstor=tube.b
dbstor=tube.db
tube.db[0,*] = 0.0
tube.db[1,*] = 0.0
if keyword_set(flat) then begin
   tube.b=b0
   return
endif
;stole another common block to be repurposed
;fieldstr is the equiv of zbot in expanse
;field len is the equiv of ztop
;fieldderiv is the proportion of expansion like rec_ang was.

;stop
l1=where(tube.x[0,*] ge (tube.x[0,0]+fieldstr) and tube.x[0,*] le (tube.x[0,tube.n/2]-fieldlen))
m=where(tube.x[0,*] ge (tube.x[0,tube.n/2]-fieldlen) and tube.x[0,*] le (tube.x[0,tube.n/2]+fieldlen))
r1=where(tube.x[0,*] le (tube.x[0,-1]-fieldstr) and tube.x[0,*] ge (tube.x[0,tube.n/2]+fieldlen))
l2=tube.x[0,max(l1)]-tube.x[0,min(l1)]
r2=tube.x[0,min(r1)]-tube.x[0,max(r1)]
tube.b=b0
if l1[0] ne -1 then tube.b[l1]=b0*(1-fieldderiv*sin((tube.x[0,l1]-tube.x[0,l1[0]])*!pi/2.0/l2)^2)
if r1[0] ne -1 then tube.b[r1]=b0*(1-fieldderiv*sin((tube.x[0,r1[-1]]-tube.x[0,r1])*!pi/2.0/r2)^2)
if m[0] ne -1 then tube.b[m]=tube.b[r1[0]]
if l1[0] ne -1 then tube.db[0,l1]=-b0*fieldderiv*sin((tube.x[0,l1]-tube.x[0,l1[0]])*!pi/2.0/l2)*cos((tube.x[0,l1]-tube.x[0,l1[0]])*!pi/2.0/l2)*!pi/l2
if r1[0] ne -1 then tube.db[0,r1]=b0*fieldderiv*sin((tube.x[0,r1[-1]]-tube.x[0,r1])*!pi/2.0/r2)*cos((tube.x[0,r1[-1]]-tube.x[0,r1])*!pi/2.0/r2)*!pi/r2
if m[0] ne -1 then tube.db[0,m]=0.0

;with the way we formulated b sheet, we need to have something to add
;to it so that its not 0 outside the sheet, else the derivative makes
;no sense
;could fix the derivative or we could just split the guide field into
;a part for the expansion and a part for the sheet.
bstor=tube.b
dbstor=tube.db
bsheet=tube.b
dbsheet=tube.db

z=ztop-zbot
bsheet=sqrt(bz*bz*(1-((tube.x[2,*]-z/2.0)*2/z)^2));
dbsheet[2,*]=-4*bz*bz*(tube.x[2,*]-z/2.0)/(z*z)          ;how do we enforce this only in sheet

edge=where(tube.x[2,*] le zbot or tube.x[2,*] ge ztop)
dbsheet[2,edge]=0.0
bsheet[edge]=0.0
;stop
;we've attempted to fix the derivative by removing the
;multiplication and division of dbsheet by bsheet in its deff and then
;in the tube.db[2,*] line, if this doesnt work, we will split the
;guide field in half and place it into bsheet and bstor equally, and
;that should solve the problem

;did not fix, i will return to this later, for now this remains broken
;where the tube goes into the current sheet, implying that the scaling
;i did makes it broken.
tube.b=sqrt(bstor^2+bsheet^2)
tube.db[0,*]=(bstor*dbstor[0,*])/tube.b
tube.db[1,*]=0.0
tube.db[2,*]=(dbsheet[2,*])/tube.b

;issue occuring in db[0,*] and b when tube.x[2,*] le zbot...at least thats
;easy to track down

;stop


return
end
; -----------------------------
pro init_expgs,rec_ang_,b0_,ztop_,zbot_,base,peak,expansion
  common cs, rec_ang,ztop,zbot,b0
common fieldshape,fieldlen,fieldstr,fieldderiv,cs_scale

;len = l0 - len
fieldderiv=expansion
fieldstr=base
fieldlen=peak
rec_ang=rec_ang_;in this case, let rec_angle be the expansion ratio desired as a fraction, eg .5 is a doubling of area, .99 is a 100 fold)
b0=b0_
zbot=zbot_;lets reuse bot as the chrom distance keep as zero for now
ztop=ztop_;reuse this as the distance from midpoint to top of expanding field right now, place this at the top of the chromosphere?
return
end
;-----------------------------------------------
