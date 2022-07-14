load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  ; Load the NCL file that contains the
begin
  y1 = new((/5,100/),float)
  y2 = new((/5,100/),float)
b1=new(5,float)
b1=(/1.05,1.15,1.25,1.35,1.45/)
b2=b1-0.95
pi=3.14159266
hbar=.5
x0 = fspan(-4*pi,4*pi,100)/pi
y0=tanh(x0)


x = fspan(0,1,100)
do i=0,4
y1(i,:) = tanh(b1(i)*(hbar+x)) - tanh(b1(i)*hbar)
end do  
do i=0,4
y2(i,:) = tanh(b2(i)*(hbar+x)) - tanh(b2(i)*hbar)
end do  

;wks = gsn_open_wks("x11","gsun01n")  ; Open an X11 workstation.
wks = gsn_open_wks("pdf","tanh")  ; Open an X11 workstation.
plot = new(4,graphic)
res=True
res@gsnDraw   = False  ; do not draw the plot (plots will be paneled)
res@gsnFrame  = False  ; do not advance the frame (plots will be paneled)
res@tiXAxisString = "x"
res@tiYAxisString = "y"
res@tiMainString  = "y=tanh(x)"
plot(0) = gsn_xy(wks,x0,y0,res)         ; Draw an XY plot with 1 curve.
;y0=tanh(x0+0.1)-tanh(x0)
;res@tiMainString  = "y=tanh(x+1)-tanh(x)"
;plot(1) = gsn_xy(wks,x0,y0,res)         ; Draw an XY plot with 1 curve.

res@pmLegendDisplayMode    = "Always"         ; turn on legend
res@pmLegendSide           = "Top"            ; Change location of 
res@pmLegendParallelPosF   = .15               ; move units right
res@pmLegendOrthogonalPosF = -0.5             ; move units down
res@pmLegendWidthF         = 0.2              ; Change width and
res@pmLegendHeightF        = 0.2             ; height of legend.
res@lgPerimOn              = False            ; turn off box around
res@lgLabelFontHeightF     = 7              ; label font height
res@xyLineThicknesses = (/3,3,4,3,3/)
res@xyDashPatterns  = (/2,3,0,12,14/)
res@tiXAxisString = "upper layer depth anomaly"
res@tiYAxisString = "tanh(b1*(hbar+h)) - tanh(b1*hbar)"
res@tiMainString      = "h>0"
res@xyExplicitLegendLabels = ""+b1
plot(1) = gsn_xy(wks,x,y1,res)         ; Draw an XY plot with 1 curve.
res@tiYAxisString = "tanh(b2*(hbar-h)) - tanh(b2*hbar)"
res@tiMainString      = "h<0" 
res@xyExplicitLegendLabels = ""+b2
plot(2) = gsn_xy(wks,x,y2,res)         ; Draw an XY plot with 1 curve.


resP             = True
resP@gsnMaximize = True
resP@gsnPanelRowSpec = True
;resP@gsnDraw   = False  ; do not draw the plot (plots will be paneled)
;resP@gsnFrame  = False  ; do not advance the frame (plots will be paneled)
gsn_panel(wks,plot,(/1,2/),resP)

end
