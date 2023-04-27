pro update_gravity,tube
  common gravell,leng,dx,dz,dl
  common tube_length,tota,chro,coro
  lef=where(tube.l le chro)
  rig=where(tube.l ge tube.l[-1]-chro)
  cent=where(tube.l ge chro and tube.l le tube.l[-1]-chro)
  tube.gpar[lef]=tube.gpar[0]
  tube.gpar[rig]=tube.gpar[-1]
  corlen=total(tube.dl[cent],/cumulative)
  corlen=corlen/corlen[-1]*leng[-1]
  cordx=interpol(dx,leng,corlen)
  cordz=interpol(dz,leng,corlen)
  cortheta=atan(cordz/cordx)
  corgrav=tube.gpar[0]*sin(cortheta)
  tube.gpar[cent]=corgrav

end
pro init_gravity,a,b,tube,pres=pres

  common tube_length,tota,chro,coro
                                ;initialize to arb ellipse, a is
                                ;width, b is height, calc the
                                ;arclength, then calc half length of
                                ;tube minus chro, assume vert for all
                                ;chro. then we scale the length to
                                ;coro. and we can get what the arc
                                ;length is, and then we get the vert
                                ;comp along the arc for time.
                                ;so something like a parameterized
                                ;curve from new transform using a,b
                                ;calc the arc length using ~1k
                                ;segments 
  
;------------ad hoc pressure equalization for grav
  chro=tube.l[(min(where(tube.t[0:tube.n/2] gt tube.t[0])))[0]]
  index=(min(where(tube.t[0:tube.n/2] gt tube.t[0])))[0]
  rindex=(max(where(tube.t[tube.n/2:-1] gt tube.t[-1])))[0]+tube.n/2
  tota=tube.l[-1]
  coro=tota-2*chro
  lef=where(tube.l le chro)
  rig=where(tube.l ge tube.l[-1]-chro)
  cent=where(tube.l ge chro and tube.l le tube.l[-1]-chro)
;------------------------------
  common gravell,leng,dx,dz,dl
  if a eq 0 and b eq 0 then begin
     dz=replicate(0.0,10001)
     dl=dindgen(10001)/10000
     dx=dl
     leng=dindgen(10001)/10000
  endif else begin
     t=findgen(n_elements(tube.l))/(n_elements(tube.l)-1)*!pi
     dt=t-shift(t,1)
     dt[0]=dt[1]
     pts=10001
     arclen=total(sqrt((a^2*sin(t)^2+b^2*cos(t)^2))*dt)
     l=arclen/(pts-1)
     x=-a
     z=0.0
     dx_=[]
     dz_=[]
     dl_=[]
     for i=1,pts-1 do begin
        dzdx1=-x[-1]/(z[-1]+0.001)*b^2/a^2
        theta=atan(dzdx1)
        h=l*cos(theta)
        k1=h*dzdx1
        xg=x[-1]+.5*h        
        zg=z[-1]+.5*k1
        dzdx2=-xg/(zg+0.001)*b^2/a^2
        theta2=atan(dzdx2)
        h2=l*cos(theta2)
        k2=h2*dzdx2
        x=[x,x[-1]+h2] 
        z=[z,z[-1]+k2]
        dx_=[dx_,h2]
        dz_=[dz_,k2]
        dl_=[dl_,sqrt(h2^2+k2^2)]
     endfor
     len_=total(dl_,/cumulative)
                                ;help,dx_,dz_,dl_,len_
                                ;stop
     lim=max(where(z gt 0.0))
     leng=len_[0:lim]
     dl=dl_[0:lim]
     dz=dz_[0:lim]
     dx=dx_[0:lim]
                                ;len is our length, and dx,dz,dl will
                                ;let us calc the tangent at any point
                                ;to multiply against gpar[0]
     
;------------------------------
     
     tube.gpar[lef]=tube.gpar[0]
     tube.gpar[rig]=tube.gpar[-1]
     corlen=total(tube.dl[cent],/cumulative)
     corlen=corlen/corlen[-1]*leng[-1]
     cordx=interpol(dx,leng,corlen)
     cordz=interpol(dz,leng,corlen)
     cortheta=atan(cordz/cordx)
     corgrav=tube.gpar[0]*sin(cortheta)
     tube.gpar[cent]=corgrav
                                ;stop
     
     
     if keyword_set(pres) then begin
        r=0.01*(1.38/1.67/tube.mpp)
        p1=fltarr(tube.n/2-index+2)
        p2=fltarr(rindex-tube.n/2+2)
        p1[0]=tube.p[index]
        p2[0]=tube.p[rindex]
        for i=index,tube.n/2,1 do begin
           hch=tube.t[i-1]*r/tube.gpar[i]
           tube.p[i]=tube.p[i-1]*exp(tube.dl[i-1]/hch)
           p1[i-index+1]=p1[i-index]*exp(tube.dl[i-1]/hch)
        endfor
        for i=1,374 do begin
           tube.p[tube.n/2+i]=tube.p[tube.n/2-i]
        endfor
        tube.rho=tube.p/r/tube.t
        set_rho,tube
     endif
     
  endelse
  
end
