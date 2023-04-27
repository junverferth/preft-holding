; set the parameters in a tube to be an RTV, constant pressure,
; equilibrium.  feet of the loop are at tmin located a distance chr from
; each end.  The apex, midway between these
; in physical length, is at tmax.

function rtv_tube_test, len, chrr=chrr, tmax=tmax, tminr=tminr, n=n, $
  debug=debug,tminl=tminl,chrl=chrl,chrom=chrom,exp_l=exp_l,exp_r=exp_r,corshape=corshape,pscale=pscale
common constr,chrn,dif_le,lengt,peak,cp_rho
common cs,rec_ang,ztop,zbot,b0;,depthl,depthr,bpeakl,bpeakr,x0l,x0r,bminl,bminr,nozl,nozr,cor_rise,cdl,cdr

if( not keyword_set(n) ) then n=400;  number of points in tube interior
if( not keyword_set(tminl) ) then tminl=0.03;  30,000K chromosphere
if( not keyword_set(tminr) ) then tminr=tminl
if( not keyword_set(tmax) ) then tmax=2.0
if( not keyword_set(chrl) ) then chrl=2.0;  length of chromopshere in Mm
if( not keyword_set(chrr) ) then chrr=chrl
if( not keyword_set(pscale) ) then pscale=1.02;keyword to allow the coronal pressure to be slid up and down.
tmin_r=tminr
tmin_l=tminl

nh=n/2;may want to make this ratio flexible later on

cs=ion_consts()
chi=cs.mpp*cs.epamu/(1.38d-10);             [ MK / erg ]
kap0=1.0d-6                   ;             [ erg / cm / s / K ]
;stop
;temp grid
fctr=0.9
nb_l=4000
nb_r=4000
nc_l=240
nc_r=240
;mid range
tb_r=(tmin_r/fctr)*exp(findgen(nb_r)*alog(fctr*fctr*tmax/tminr)/float(nb_r-1))
tb_l=(tmin_l/fctr)*exp(findgen(nb_l)*alog(fctr*fctr*tmax/tminl)/float(nb_l-1))
na_r=ceil(2*(tb_r[0]-tmin_r)/(tb_r[1]-tb_r[0]))
na_l=ceil(2*(tb_l[0]-tmin_l)/(tb_l[1]-tb_l[0]))
ta_r=tmin_r+(tb_r[0]-tmin_r)*(findgen(na_r)/float(na_r-1))^2
ta_l=tmin_l+(tb_l[0]-tmin_l)*(findgen(na_l)/float(na_l-1))^2
nc_r=ceil(2*(tmax-tb_r[nb_r-1])/(tb_r[nb_r-1]-tb_r[nb_r-2]))
nc_l=ceil(2*(tmax-tb_l[nb_l-1])/(tb_l[nb_l-1]-tb_l[nb_l-2]))
tc_r=tmax+(tb_r[nb_r-1]-tmax)*(findgen(nc_r)/float(nc_r-1))^2
tc_r=reverse(tc_r)
tc_l=tmax+(tb_l[nb_l-1]-tmax)*(findgen(nc_l)/float(nc_l-1))^2
tc_l=reverse(tc_l)
t_r=[ta_r,tb_r[1:nb_r-2],tc_r]
t_l=[ta_l,tb_l[1:nb_l-2],tc_l]
nt_r=n_elements(t_r)
nt_l=n_elements(t_l)
dlt_r=2*(t_r-shift(t_r,1))/(t_r+shift(t_r,1))
dlt_r[0]=dlt_r[1]
dlt_l=2*(t_l-shift(t_l,1))/(t_l+shift(t_l,1))
dlt_l[0]=dlt_l[1]
;computations
lam_r=lam_rlf(t_r*1.0d6) ;[ erg cm^3 / s ] 
lam_l=lam_rlf(t_l*1.0d6)
f_r=0.1*total(t_r^(1.5)*lam_r*dlt_r/kap0,/cum); [1.0e10 cm^4 K^5]
f_l=0.1*total(t_l^(1.5)*lam_l*dlt_l/kap0,/cum)
rt_r=(f_r-f_r[nt_r-1]*(t_r^(3.5)-tmin_r^(3.5))/(tmax^(3.5)-tmin_r^(3.5))) > 1.0d-22
rt_l=(f_l-f_l[nt_l-1]*(t_l^(3.5)-tmin_l^(3.5))/(tmax^(3.5)-tmin_l^(3.5))) > 1.0d-22
xi_r=100*0.707107*t_r^(2.5)/sqrt(rt_r);[1.0e8 cm^-2]
xi_l=100*0.707107*t_l^(2.5)/sqrt(rt_l)
ig_r=0.5*(xi_r*dlt_r+shift(xi_r*dlt_r,1))
ig_r[0]=0.5*(xi_r[0]*dlt_r[0])
ig_l=0.5*(xi_l*dlt_l+shift(xi_l*dlt_l,1))
ig_l[0]=0.5*(xi_l[0]*dlt_l[0])
col_e_r=total(ig_r,/cum);[1.0e8 cm^-2]
col_e_l=total(ig_l,/cum)
p_r=(2.0/chi/len)*total(t_r*xi_r*dlt_r);[erg/cm^3];strictly speaking, i could tune this value to reduce density?
p_r=pscale*p_r
p_l=(2.0/chi/len)*total(t_l*xi_l*dlt_l)
p_l=pscale*p_l
;stop
;create the grid
col_norm_r=col_e_r/max(col_e_r)
col_norm_l=col_e_l/max(col_e_l)
grid_dm_r=(dindgen(nh)+3.0)^(0.5)
grid_dm_l=(dindgen(nh)+3.0)^(0.5) ;leave in dup in case we switch nh later
grid_col_r=[0.0,total(grid_dm_r[0:(nh-2)],/cum)]
grid_col_l=[0.0,total(grid_dm_l[0:(nh-2)],/cum)]
tgh_r=interpol(t_r,col_norm_r,grid_col_r/grid_col_r[nh-1])
tgh_l=interpol(t_l,col_norm_l,grid_col_l/grid_col_l[nh-1])
dlh_r=tgh_r*grid_dm_r
dlh_l=tgh_l*grid_dm_l
l_tot_r=total(dlh_r)
l_tot_l=total(dlh_l)
;stop
dlh_r=0.5*len*dlh_r/l_tot_r
dlh_l=0.5*len*dlh_l/l_tot_l
;arrays for central portion
dl_c=[dlh_l,reverse(dlh_r)]
t_c=[tgh_l,reverse(tgh_r)]
;chromosphere
if ~keyword_set(exp_r) then chr_exp_r=0.5 else chr_exp_r=exp_r
if ~keyword_set(exp_l) then chr_exp_l=0.5 else chr_exp_l=exp_l
nchr_r=ceil(3.0^(chr_exp_r)*((chr_exp_r+1.0)*chrr/dlh_r[0])^(1.0/(chr_exp_r+1.0)))
nchr_l=ceil(3.0^(chr_exp_l)*((chr_exp_l+1.0)*chrl/dlh_l[0])^(1.0/(chr_exp_l+1.0)))
dl_chr_r=dlh_r[0]*((3.0+findgen(nchr_r))/3.0)^(chr_exp_r)
dl_chr_l=dlh_l[0]*((3.0+findgen(nchr_l))/3.0)^(chr_exp_l)
dl=[reverse(dl_chr_l),dl_c,dl_chr_r]
l=[0.0,total(dl,/cum)]
;stop
n_tot=n_elements(l)
n_t_tot=n_elements(t_c)
t_chr_r=t_c[n_t_tot-1]
t_chr_l=t_c[0]
t_arr=[replicate(t_chr_l,nchr_l),t_c,replicate(t_chr_r,nchr_r)]
z_chr_r=l[(n_tot-nchr_r):(n_tot-1)]-l[n_tot-nchr_r-1]
z_chr_l=l[0:nchr_l-1]-l[nchr_l]
gp=2.74d-4
r=0.01*(1.38/1.67/cs.mpp)
h_ch_r=t_chr_r*r/gp
h_ch_l=t_chr_l*r/gp
p_chr_r=p_r*exp(z_chr_r/h_ch_r)
p_chr_l=p_l*exp(-z_chr_l/h_ch_l)
;stop
tube=make_tube(n_tot);look into making sure that n_tot is odd so symmetry as best i can

tube.x[0,*]=l
tube.x[2,*]=zbot
tube.p=[p_chr_l,replicate(p_l,n_elements(dl_c)),p_chr_r]

tube.t=[t_arr,tmin_r]
tube.p[n_tot-1]=tube.p[n_tot-2]
tube.gpar[0:nchr_l]=-gp
tube.gpar[(n_tot-nchr_r):(n_tot-1)]=gp

tube.rho=tube.p/r/tube.t
;stop
set_rho,tube;,chrom=chrom,corshape=corshape
calc_len,tube
set_heat_for_equilib,tube

chrn=nchr_l
dif_le=dl
lengt=l
;save,/variables,filename='rtv_tube_save.sav'
;stop
return, tube
end
