;preft batch plotting
pro stackrun,stack1,name,eps=eps

xpos=[.08,.26,.3,.48,.52,.7,.74,.92]
ypos=[.1,.9]
cbpos=[.94,.1,.96,.9]

if keyword_set(eps) then jups_open,name+'preft_comp',1600,900

loadct,70
cgimage,stack1.vel_x,xvector=stack1.x,title='vel +- 250km/s',position=[xpos[0],ypos[0],xpos[1],ypos[1]],charsize=0.7,/axes,minv=-.25,maxv=.25

loadct,39
cgimage,alog10(stack1.den),xvector=stack1.x,title='e- den log8.5-11',position=[xpos[2],ypos[0],xpos[3],ypos[1]],charsize=0.7,/axes,minv=8.5,maxv=11,/noerase
cgimage,alog10(stack1.temp*1d6),xvector=stack1.x,title='temp log4-7.5',position=[xpos[4],ypos[0],xpos[5],ypos[1]],charsize=0.7,/axes,minv=4,maxv=7.5,/noerase
cgimage,alog10(stack1.pres),xvector=stack1.x,title='pressure log-1-2',position=[xpos[6],ypos[0],xpos[7],ypos[1]],charsize=0.7,/axes,minv=-1,maxv=2,/noerase

cgtext,.525,.02,'Length from Center [Mm]',/normal,alignment=0.5,charsize=1.0,orientation=0.0
cgtext,.03,.5,'Time [s]',/normal,alignment=0.5,charsize=1.0,orientation=90.0

if keyword_set(eps) then jups_close

end
