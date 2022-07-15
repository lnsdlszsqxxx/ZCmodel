;plot var patterns for verification in ZC
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
H1 = new((/nt,ny,nx/),"float") ;H1 H2 in ZC
T = new((/nt,ny,nx/),"float") ;T is SST
TAUx = new((/nt,ny,nx/),"float")
TAUy = new((/nt,ny,nx/),"float")

suffix=yrStart+"-"+yrEnd
histFnU="ZCHIST_SSTA_34x30_Deg.bin"
histFnV="ZCHIST_SSTA_34x30_Deg.bin"
histFnH1="ZCCPD_H1_34x30_m.bin"
histFnT="ZCHIST_SSTA_34x30_Deg.bin"
histFnTAUx="ZCHIST_TAUxA_34x30_dyne.bin"
histFnTAUy="ZCHIST_SSTA_34x30_Deg.bin"

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
H1 = fbindirread(histFnH1,0,(/nt,ny,nx/),"float")

H1@lat2d=lat2d
H1@lon2d=lon2d
H1@long_name = "Upper ocean depth anom"
H1@units     = "Unit meter"

printMinMax(H1,0)

;Tclmt=dim_avg_n_Wrap(T,0)

wks = gsn_open_wks("x11","sst")   ;Open a workstation.
;wks = gsn_open_wks("pdf","./H1_evolution_yr12-13") ; Open a workstation.
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

res@tmXBLabelFontHeightF = 0.02 ;lat lon label size
res@tmYLLabelFontHeightF = 0.02
;res@tmXBTickSpacingF=30 ;tick space
res@tmYLTickSpacingF=15

res@gsnLeftStringFontHeightF=0.02 ;topleft string size
res@gsnRightStringFontHeightF=0.02 ;topright string size

res@cnLevelSelectionMode = "ManualLevels"   
a= -100.
b= 100.
res@cnMinLevelValF       = a
res@cnMaxLevelValF       = b
;res@cnLevelSpacingF      = (b-a)/15.0
res@cnLevelSpacingF      = 10
res@lbLabelFontHeightF  = 0.01    ; make labels smaller
;;res@lbOrientation = "vertical"


;plot(0) = gsn_csm_contour_map(wks,Tclmt(:,:),res)

iyr12Jan=132 ;year 12, january, this great peak
iyr13Dec=155
icounter=0
month=(/"Jan(0)","Apr(0)","Jul(0)","Oct(0)","Jan(1)","Apr(1)","Jul(1)","Oct(1)"/)
do i=iyr12Jan,iyr13Dec,3
T@long_name = month(icounter)
plot(icounter) = gsn_csm_contour_map(wks,H1(i,:,:),res)
icounter=icounter+1
end do

ElNino="El-Nin~H-13V2F35~D~FV-2H3~o"
LaNina="La-Nin~H-13V2F35~D~FV-2H3~a"

resP             = True
resP@gsnMaximize = True
resP@txString  = "H anomalies evolution during an "+ElNino+" Event in Year 12"
resP@gsnPanelLabelBar    = True ; add common colorbar
;resP@lbOrientation       = "vertical" ; vertical label bar
resP@lbLabelFontHeightF  = 0.015              ; make labels smaller
gsn_panel(wks,(/plot(0),plot(4),plot(1),plot(5),plot(2),plot(6),plot(3),plot(7)/),(/4,2/),resP)

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
