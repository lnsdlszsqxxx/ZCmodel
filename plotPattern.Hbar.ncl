;plot var patterns for verification in ZC
;H1 is depth anom
;liang, 2022.6.30

begin

yrStart=1
yrEnd=90

nx=34
ny=30
nyear=yrEnd-yrStart+1
nt=(yrEnd-yrStart+1)*12 

U = new((/nt,ny,nx/),"float")
V = new((/nt,ny,nx/),"float")
H = new((/nt,ny,nx/),"float") ;H1 H2 in ZC
T = new((/nt,ny,nx/),"float") ;T is Hbar
TAUx = new((/nt,ny,nx/),"float")
TAUy = new((/nt,ny,nx/),"float")

suffix=yrStart+"-"+yrEnd
histFnH="ZCCPD_H1_34x30_m.bin"

lat = new((/ny/),"float")
lon = new((/nx/),"float")
lon(:)=fspan(101.25,286.875,nx)
lat(:)=fspan(-29,29,ny)

lat2d = new((/ny,nx/),"float")
lon2d = new((/ny,nx/),"float")
do i=0,nx-1
lat2d(:,i)=lat
end do
do i=0,ny-1
lon2d(i,:)=lon
end do

;--------------
;read data
;--------------
H1 = fbindirread(histFnH,0,(/nt,ny,nx/),"float")
H1@lat2d=lat2d
H1@lon2d=lon2d
H1@long_name = "anom upper layer depth"
H1@units     = "unit meter"

print(max(H1))
print(min(H1))
printMinMax(H1,1)
exit

T(0,0,:)=(/1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.65,1.7,1.75,1.75,1.75,1.5,1.25,1.,.932,.864,.796,.728,.66,.592,.524,.5,.5,.5,.5,.5/) ;this is from ssta.f 

do j=1,29
T(0,j,:)=T(0,j-1,:)
end do

T@lat2d=lat2d
T@lon2d=lon2d
T@long_name = "mean upper layer depth"
T@units     = "unit 100 meter"

Tclmt=dim_avg_n_Wrap(T,0)

wks = gsn_open_wks("x11","sst")   ;Open a workstation.
;wks = gsn_open_wks("pdf","Hbar_prescribed") ; Open a workstation.
gsn_define_colormap(wks,"NCV_blue_red")
;gsn_define_colormap(wks,"grads_rainbow")
plot = new(8,graphic)

lat1=-20
lat2=20
lon1=129
lon2=276

res                      = True
res@mpLimitMode = "LatLon"    ; Limit the map view.
res@mpMinLatF            = lat1      ; range to zoom in on
res@mpMaxLatF            = lat2
res@mpMinLonF            = lon1
res@mpMaxLonF            = lon2
res@mpGridAndLimbOn      = True
res@mpGridLineDashPattern  = 2          ; lat/lon lines as dashed
;res@mpProjection = "CylindricalEquidistant"  ; Change the map projection.
res@mpCenterLonF = (lon1+lon2)/2
res@mpCenterLatF = (lat1+lat2)/2
;
;res@vpXF      = 0.1   ; Increase size and change location
;res@vpYF      = 0.92  ; of vector plot.
;res@vpWidthF  = 0.75
;res@vpHeightF = 0.75

;res@mpLimitMode       = "Corners"  ; Zoom in on the plot area.
;res@mpLeftCornerLonF  = lon1
;res@mpRightCornerLonF = lon2
;res@mpLeftCornerLatF  = lat1
;res@mpRightCornerLatF = lat2
;
;;res@mpPerimOn         =   True    ; Turn on map perimeter.
res@gsnAddCyclic = False
res@gsnDraw  = False                          ; don't draw
res@gsnFrame = False                          ; don't advance frame
res@cnFillOn            = True            ; turn on color
;res@cnLinesOn            = True
;;res@cnLineColor          = "white"
res@cnLineLabelsOn       = True

res@lbLabelBarOn        = False ;turn off individual cb's

res@tmXBLabelFontHeightF = 0.012 ;lat lon label size
res@tmYLLabelFontHeightF = 0.012
;res@tmXBTickSpacingF=30 ;tick space
res@tmYLTickSpacingF=15

res@gsnLeftStringFontHeightF=0.02 ;topleft string size
res@gsnRightStringFontHeightF=0.02 ;topright string size

res@cnLevelSelectionMode = "ManualLevels"   
a= .5
b= 1.8
res@cnMinLevelValF       = a
res@cnMaxLevelValF       = b
;res@cnLevelSpacingF      = (b-a)/15.0
res@cnLevelSpacingF      = 0.1
res@lbLabelFontHeightF  = 0.01    ; make labels smaller
;;res@lbOrientation = "vertical"

plot(0) = gsn_csm_contour_map(wks,T(0,:,:),res)

ElNino="El-Nin~H-13V2F35~D~FV-2H3~o"
LaNina="La-Nin~H-13V2F35~D~FV-2H3~a"

resP             = True
resP@gsnMaximize = True
resP@txString  = "~Z120~Prescribed Hbar in ZC model"
resP@gsnPanelLabelBar    = True ; add common colorbar
;resP@lbOrientation       = "vertical" ; vertical label bar
resP@lbLabelFontHeightF  = 0.015              ; make labels smaller
gsn_panel(wks,plot(0),(/1,1/),resP)

;**********************************
; now maximize the plot
;**********************************   
; Maximize these plots for PS or PDF output. If output is
; to an X11 window or an NCGM file, then this routine doesn't
; do anything.
psres = True
;maximize_output(wks,psres)  ; calls draw and frame for you 


;;---------------------
;;save data
;;---------------------
;print("saving...")
;;---Remove file just in case
;filename1 ="./restoringSST_ClmtSST_from_run_clmtwind_DegC.bin"
;system("rm -f " + filename1)
;;fbindirwrite(filename1,sstclmt)
;setfileoption("bin","WriteByteOrder","BigEndian")
;fbinwrite(filename1,sstclmt)
;system("ls -l " + filename1)


end
