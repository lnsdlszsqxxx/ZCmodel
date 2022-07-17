;plot Tsub in ZC based on ZC1987 Appendix
;plot1: linearized Tsub
;plot2: Tsub for diff Hbar
;plot3: Tsub for diff parameters
;liang, 2022.7.15

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  ; Load the NCL file that contains the
begin
pi=3.14159266
b1=1.25
b2=3.0
T1=28.
T2=-40.
Hbar=1.0 ;in the model, Hbar unit is 100 meter, 1.0 is the mean Hbar

;define parameter value ranges.
varb1=new(5,float)
varb2=new(5,float)
varT1=new(5,float)
varT2=new(5,float)
varHbar=new(5,float)
varb1=(/1.05,1.15,1.25,1.35,1.45/)
varb2=(/2,2.5,3,3.5,4/)
varT1=(/24,26,28,30,32/)
varT2=(/-60,-50,-40,-30,-20/)
varHbar=(/0.6,0.8,1.0,1.2,1.4/)


;in 90 yr, h min=-87.6847,max=122.261,therefore here I use -80~80. dx=0.5
nh=321 ;dx=0.5
h = fspan(-80,80,nh)
Tsub = new((/nh/),float)
Tsub1 = new((/5,nh/),float)
Tsub1 = new((/5,nh/),float)
y = new((/2,nh/),float)
y1 = new((/5,nh/),float)
y2 = new((/5,nh/),float)
y3 = new((/5,nh/),float)
y4 = new((/5,nh/),float)
y5 = new((/5,nh/),float)

;plot 1__________________________________________________________________
;define Tsub based on ZC1987 Eq. A13
Tsub(nh/2+1:nh-1)=T1*( tanh(b1*(Hbar+.01*h(nh/2+1:nh-1))) - tanh(b1*Hbar) ) ;h>0
Tsub(0:nh/2-1)=T2*( tanh(b2*(Hbar-.01*h(0:nh/2-1))) - tanh(b2*Hbar) );h<0
Tsub(nh/2)=( Tsub(nh/2-1) + Tsub(nh/2+1) )/2.0 ;h=0, take the mean

;linear regretion, least squares method
rc = regline (h,Tsub)
yReg = rc*h + rc@yintercept
y(0,:)=Tsub
y(1,:)=yReg
;------------------------------------------------------------------------

;plot 2__________________________________________________________________
do i=0,4
Hbar=varHbar(i) ;change Hbar
Tsub(nh/2+1:nh-1)=T1*( tanh(b1*(Hbar+.01*h(nh/2+1:nh-1))) - tanh(b1*Hbar) )
Tsub(0:nh/2-1)=T2*( tanh(b2*(Hbar-.01*h(0:nh/2-1))) - tanh(b2*Hbar) )
Tsub(nh/2)=( Tsub(nh/2-1) + Tsub(nh/2+1) )/2.0
y1(i,:)=Tsub
end  do
Hbar=1.0
;------------------------------------------------------------------------

Hbar=.5
;plot 3__________________________________________________________________
do i=0,4
T1=varT1(i) ;change T1
Tsub(nh/2+1:nh-1)=T1*( tanh(b1*(Hbar+.01*h(nh/2+1:nh-1))) - tanh(b1*Hbar) )
Tsub(0:nh/2-1)=T2*( tanh(b2*(Hbar-.01*h(0:nh/2-1))) - tanh(b2*Hbar) )
Tsub(nh/2)=( Tsub(nh/2-1) + Tsub(nh/2+1) )/2.0
y2(i,:)=Tsub
end  do
T1=28
;------------------------------------------------------------------------

;plot 4__________________________________________________________________
do i=0,4
T2=varT2(i) ;change T2
Tsub(nh/2+1:nh-1)=T1*( tanh(b1*(Hbar+.01*h(nh/2+1:nh-1))) - tanh(b1*Hbar) )
Tsub(0:nh/2-1)=T2*( tanh(b2*(Hbar-.01*h(0:nh/2-1))) - tanh(b2*Hbar) )
Tsub(nh/2)=( Tsub(nh/2-1) + Tsub(nh/2+1) )/2.0
y3(i,:)=Tsub
end  do
T2=-40
;------------------------------------------------------------------------

;plot 5__________________________________________________________________
do i=0,4
b1=varb1(i) ;change b1
Tsub(nh/2+1:nh-1)=T1*( tanh(b1*(Hbar+.01*h(nh/2+1:nh-1))) - tanh(b1*Hbar) )
Tsub(0:nh/2-1)=T2*( tanh(b2*(Hbar-.01*h(0:nh/2-1))) - tanh(b2*Hbar) )
Tsub(nh/2)=( Tsub(nh/2-1) + Tsub(nh/2+1) )/2.0
y4(i,:)=Tsub
end  do
b1=1.25
;------------------------------------------------------------------------

;plot 6__________________________________________________________________
do i=0,4
b2=varb2(i) ;change b2
Tsub(nh/2+1:nh-1)=T1*( tanh(b1*(Hbar+.01*h(nh/2+1:nh-1))) - tanh(b1*Hbar) )
Tsub(0:nh/2-1)=T2*( tanh(b2*(Hbar-.01*h(0:nh/2-1))) - tanh(b2*Hbar) )
Tsub(nh/2)=( Tsub(nh/2-1) + Tsub(nh/2+1) )/2.0
y5(i,:)=Tsub
end  do
b2=3.0
;------------------------------------------------------------------------



wks = gsn_open_wks("x11","gsun01n")  ; Open an X11 workstation.
;wks = gsn_open_wks("pdf","tanh_Tsub")  ; Open an X11 workstation.
plot = new(6,graphic)
res=True
res@gsnDraw   = False  ; do not draw the plot (plots will be paneled)
res@gsnFrame  = False  ; do not advance the frame (plots will be paneled)

res@tiXAxisString = "h (Upper ocean depth anomaly, meter)"
res@tiYAxisString = "T~B~sub~N~ (~S~o~N~C)"
res@tiMainString  = "T~B~sub~N~=T1{tanh[b1(H~B~bar~N~+.01h)]-tanh(b1*H~B~bar~N~)}, h>0~C~T~B~sub~N~=T2{tanh[b2(H~B~bar~N~ -.01h)]-tanh(b2*H~B~bar~N~)}, h<0~C~T1=28~S~o~N~C,T2=-40~S~o~N~C,b1=1.25,b2=3.0,H~B~bar~N~=1.0"
res@xyLineThicknesses = (/3.0,3.0/)
res@xyDashPatterns  = (/0,1/)

plot(0) = gsn_xy(wks,h,y,res)         ; Draw an XY plot with 1 curve.

;this is another way for legend, called simple legend
gres = True
gres@YPosPercent = 90.    ; expressed as %, 0->100
gres@XPosPercent = 5      ; expressed as %, 0->100
lineres = True
lineres@lgLineColors = (/"black","black"/) ; line colors
lineres@lgLineThicknesses = 3.0                        ; line thicknesses
lineres@lgDashIndexes= (/0,1/)
lineres@LineLengthPercent = 9.            ; expressed as %, 0->100
textres = True
textres@lgLabels = (/"ZC model T~B~sub~N~","linearized T~B~sub~N~"/)
textres@lgLabelFontHeights=0.022
plot(0) = simple_legend(wks,plot(0),gres,lineres,textres)
;draw(plot(0))
;frame(wks)
;--------------------------------------------------------

;plot 2_______________________________________________________
res2=True
res2@gsnDraw   = False  ; do not draw the plot (plots will be paneled)
res2@gsnFrame  = False  ; do not advance the frame (plots will be paneled)

res2@tiXAxisString = "h (m)"
res2@tiYAxisString = "T~B~sub~N~ (~S~o~N~C)"
res2@tiMainString  = "T~B~sub~N~ for different Hbar"
res2@xyLineThicknesses = (/3,3,3,3,3/)
res2@xyDashPatterns  = (/2,1,0,3,16/)
plot(1) = gsn_xy(wks,h,y1,res2)         ; Draw an XY plot with 1 curve.

;this is another way for legend, called simple legend
gres2 = True
gres2@YPosPercent = 90.    ; expressed as %, 0->100
gres2@XPosPercent = 5      ; expressed as %, 0->100
lineres2 = True
lineres2@lgLineThicknesses = 3.0                        ; line thicknesses
lineres2@lgDashIndexes= (/2,1,0,3,16/)
lineres2@LineLengthPercent = 9.            ; expressed as %, 0->100
textres2 = True
textres2@lgLabels = "Hbar="+varHbar
textres2@lgLabelFontHeights=0.022
plot(1) = simple_legend(wks,plot(1),gres2,lineres2,textres2)
;draw(plot(1))
;frame(wks)
;--------------------------------------------------------------------

;plot 3_______________________________________________________
res2@tiMainString  = "T~B~sub~N~ for different T1"
plot(2) = gsn_xy(wks,h,y2,res2)         ; Draw an XY plot with 1 curve.
textres2@lgLabels = "T1="+varT1
plot(2) = simple_legend(wks,plot(2),gres2,lineres2,textres2)
;draw(plot(2))
;frame(wks)
;--------------------------------------------------------------------

;plot 4_______________________________________________________
res2@tiMainString  = "T~B~sub~N~ for different T2"
plot(3) = gsn_xy(wks,h,y3,res2)         ; Draw an XY plot with 1 curve.
textres2@lgLabels = "T2="+varT2
plot(3) = simple_legend(wks,plot(3),gres2,lineres2,textres2)
;draw(plot(3))
;frame(wks)
;--------------------------------------------------------------------

;plot 5_______________________________________________________
res2@tiMainString  = "T~B~sub~N~ for different b1"
plot(4) = gsn_xy(wks,h,y4,res2)         ; Draw an XY plot with 1 curve.
textres2@lgLabels = "b1="+varb1
plot(4) = simple_legend(wks,plot(4),gres2,lineres2,textres2)
;draw(plot(4))
;frame(wks)
;--------------------------------------------------------------------

;plot 6_______________________________________________________
res2@tiMainString  = "T~B~sub~N~ for different b2"
plot(5) = gsn_xy(wks,h,y5,res2)         ; Draw an XY plot with 1 curve.
textres2@lgLabels = "b2="+varb2
plot(5) = simple_legend(wks,plot(5),gres2,lineres2,textres2)
;draw(plot(5))
;frame(wks)
;--------------------------------------------------------------------

resP             = True
resP@gsnMaximize = True
;resP@gsnPanelRowSpec = True
;resP@gsnDraw   = False  ; do not draw the plot (plots will be paneled)
;resP@gsnFrame  = False  ; do not advance the frame (plots will be paneled)
gsn_panel(wks,plot(2:5),(/2,2/),resP)

end
