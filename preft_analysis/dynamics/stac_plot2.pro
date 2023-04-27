function stac_plot2,tarr,dx=dx,nx=nx
  if ~keyword_set(dx) then dx=.120
  if ~keyword_set(nx) then nx=170.0/dx;ceil((max(tarr[0].l)-min(tarr[0].l))/dx*1.1)
  x=(findgen(nx)-(nx-1)/2.0)*dx
  n=n_elements(tarr)
  imagepres=fltarr(nx,n)
  imageb=imagepres
  imagerho=imagepres
  imagetemp=imagepres
  imageden=imagepres
  imagevelx=imagepres
  imagevely=imagepres
  minval=fltarr(n_elements(tarr))
  cut=intarr(n_elements(tarr))
  cut2=cut
  
  for i=0,n_elements(tarr)-1 do minval[i]=tarr[i].l[-1]/2.0
  for i=0,n_elements(tarr)-1 do cut[i]=max(where(x le -minval[i]))
  for i=0,n_elements(tarr)-1 do cut2[i]=min(where(x ge minval[i]))
  cut[where(cut eq -1)]=0
  cut2[where(cut2 eq -1)]=nx-1
     for i=0, n-1 do begin
        lh=tarr[i].l[tarr[i].n-1]/2.0
        n_e=tarr[i].rho*tarr[i].epamu/1.67d-8
        imageden[*,i]=interpol(n_e,tarr[i].l-lh,x)
        imageden[0:cut[i],i]=0.0
        imageden[cut2[i]:n_elements(x)-1,i]=0.0
     endfor
     for i=0, n-1 do begin
        lh=tarr[i].l[tarr[i].n-1]/2.0
        imagetemp[*,i]=interpol(tarr[i].t,tarr[i].l-lh,x)
        imagetemp[0:cut[i],i]=0.0
        imagetemp[cut2[i]:n_elements(x)-1,i]=0.0
     endfor
     for i=0, n-1 do begin
        lh=tarr[i].l[tarr[i].n-1]/2.0
        imagepres[*,i]=interpol(tarr[i].p,tarr[i].l-lh,x)
        imagepres[0:cut[i],i]=0.0
        imagepres[cut2[i]:n_elements(x)-1,i]=0.0
     endfor
     for i=0, n-1 do begin
        lh=tarr[i].l[tarr[i].n-1]/2.0
        imageb[*,i]=interpol(tarr[i].b,tarr[i].l-lh,x)
        imageb[0:cut[i],i]=0.0
        imageb[cut2[i]:n_elements(x)-1,i]=0.0
     endfor
     for i=0, n-1 do begin
        lh=tarr[i].l[tarr[i].n-1]/2.0
        imagerho[*,i]=interpol(tarr[i].rho,tarr[i].l-lh,x)
        imagerho[0:cut[i],i]=0.0
        imagerho[cut2[i]:n_elements(x)-1,i]=0.0
     endfor
     for i=0,n-1 do begin
        lh=tarr[i].l[tarr[i].n-1]/2.0
        imagevelx[*,i]=interpol(tarr[i].v[0,*],tarr[i].l-lh,x)
        imagevelx[0:cut[i],i]=0.0
        imagevelx[cut2[i]:n_elements(x)-1,i]=0.0
        imagevely[*,i]=interpol(tarr[i].v[2,*],tarr[i].l-lh,x)
        imagevely[0:cut[i],i]=0.0
        imagevely[cut2[i]:n_elements(x)-1,i]=0.0
     endfor
  ;sample=indgen(n_elements(tarr)-1001)+1000
  sample=0.0;[indgen(101)*10,sample];[0,100,200,300,400,500,600,700,800,900,sample]
  time=tarr[sample].time

                                ;edge check and fix

  output={den:imageden,temp:imagetemp,pres:imagepres,rho:imagerho,b:imageb,vel_x:imagevelx,vel_z:imagevely,time:time,x:x,sample:sample}
  return,output
end
