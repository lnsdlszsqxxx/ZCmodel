;plot NINO3.4 index for verification
;compare with obs
;liang, 2020.9.22

begin

yrStart=200
yrEnd=261

nx=166
ny=80
nz=7   ;sigma coordiante
nnz=11 ;for z coordinate
nt=(yrEnd-yrStart+1)*12 ;62 years
nxp1=nx+1
nyp1=ny+1

U = new((/nt,nnz,nyp1,nxp1/),"float")
V = new((/nt,nnz,nyp1,nxp1/),"float")
H = new((/nt, nz,nyp1,nxp1/),"float")
T = new((/nt,nnz,nyp1,nxp1/),"float")
TAUx =  new((/nt,nyp1,nxp1/),"float")
TAUy =  new((/nt,nyp1,nxp1/),"float")
Q =     new((/nt,nyp1,nxp1/),"float")
sst =   new((/nt,nyp1,nxp1/),"float") ;(nyp1,nxp1)=(lat,lon)

suffix=yrStart+"-"+yrEnd
histFnU="history.littleendian.direct.monthly.U.ms-1.z."+suffix
histFnV="history.littleendian.direct.monthly.V.ms-1.z."+suffix
histFnH="history.littleendian.direct.monthly.H.m."+suffix
histFnT="history.littleendian.direct.monthly.T.DegC.z."+suffix
histFnTAUx="history.littleendian.direct.monthly.TAUx.dynescm-2."+suffix
histFnTAUy="history.littleendian.direct.monthly.TAUy.dynescm-2."+suffix
histFnQ="history.littleendian.direct.monthly.Q.Wm-2."+suffix

lat = new((/nyp1/),"float")
lon = new((/nxp1/),"float")
lat2d = new((/nyp1,nxp1/),"float")
lon2d = new((/nyp1,nxp1/),"float")

;;;;;;;;;;;;;;;;;
;linear lon
;nonlinear lat, lat is from model log file
;;;;;;;;;;;;;;;;;
lon(:)=fspan(123,289,nxp1)
lat(:)=(/-30.0000, -28.9966, -27.9945, -26.9939, -25.9950, -24.9, -24.0028, -23.0100, -22.0197, -21.0324, -20.0484, -19.0681, -18.0921, -17.1211, -16.1558, -15.1, -14.2462, -13.3045, -12.3735, -11.4554, -10.5525, -9.6678, -8.8048, -7.9678, -7.1615, -6.3, -5.6625, -4.9805, -4.3495, -3.7721, -3.2486, -2.7772, -2.3538, -1.9729, -1.6285, -1.3, -1.0243, -0.7527, -0.4946, -0.2452, 0.0000, 0.2452, 0.4946, 0.7527, 1.0243, 1.3, 1.6285, 1.9729, 2.3538, 2.7772, 3.2486, 3.7721, 4.3495, 4.9805, 5.6625, 6.3, 7.1615, 7.9678, 8.8048, 9.6678, 10.5525, 11.4554, 12.3735, 13.3045, 14.2462, 15.1, 16.1558, 17.1211, 18.0921, 19.0681, 20.0484, 21.0324, 22.0197, 23.0100, 24.0028, 24.9, 25.9950, 26.9940, 27.9945, 28.9966, 30.0000/)


do i=0,nxp1-1
lat2d(:,i)=lat
end do
do i=0,nyp1-1
lon2d(i,:)=lon
end do

;----------------------
;calculate NINO index
;----------------------
NINOindex="NINO3"
if(NINOindex.eq."NINO34") then
  latS  = -5.0
  latN  =  5.0
  lonL  = 190.0
  lonR  = 240.0
else if(NINOindex.eq."NINO3") then
  latS  = -5.0
  latN  =  5.0
  lonL  = 210.0
  lonR  = 270.0
else
print("wrong NINOindex")
exit
end if
end if

index_totalwind=new((nt),"float",999)
filename1 ="./NINO3_DegC.62yr.forced_totalwind.bin"
index_totalwind=fbindirread(filename1,0,(/nt/),"float")

index_totalwind240=new((nt),"float",999)
filename1 ="./NINO3_DegC.62yr.forced_totalwind.bin.rerunFromCNTL240"
index_totalwind240=fbindirread(filename1,0,(/nt/),"float")

index_CNTL=fbindirread("/Users/liang/pac.Tao/OCN2/tmp/run_forced_CNTL/NINO3_DegC.62yr.forced_CNTL.bin",0,(/nt/),"float")

;------------------------------------------------------------
print("read obs nino index")
ncol=10
nRowSkip=1
data=readAsciiTable("ersst4.nino.mth.81-10.ascii.txt",ncol,"float",nRowSkip)

if(NINOindex.eq."NINO34") then
index_obs=data(0:nt-1,8)
else if(NINOindex.eq."NINO3") then
index_obs=data(0:nt-1,4)
end if
end if

;as the model has no seasonal cycle in it, we need to remove the seasonal cycle in obs as well.
print("remove obs and totalwind seasonal cycle")
index_obs=rmAnnCycle1D(index_obs)+avg(index_obs)
index_totalwind=rmAnnCycle1D(index_totalwind)+avg(index_totalwind)

;------------------------------------------------------------
;print("read TAO's NINO3 obs data for verification")
;index_TAO=new((/nt/),"float")
;index_TAO(0:711)=fbindirread("gisst-nino3-ssta-1950.01-2011.12.dat",0,(/712/),"float")
;print(index_TAO)
;exit

yyyymm   = yyyymm_time(1950, 2011, "integer")
yyyyfrac = yyyymm_to_yyyyfrac(yyyymm,0)   ; create fractional years for 

indices1 = new((/2,nt/),"float")
indices1(0,:)=index_CNTL
indices1(1,0:12*47-1)=index_totalwind(0:12*47-1)

indices2 = new((/2,nt/),"float")
indices2(0,12*41:12*62-1)=index_totalwind240(12*41:12*62-1)
indices2(1,0:12*47-1)=index_totalwind(0:12*47-1)

indices3 = new((/2,nt/),"float")
indices3(0,:)=index_obs
indices3(1,:)=index_totalwind

indices4 = new((/3,nt/),"float")
indices4(0,:)=index_obs
indices4(1,:)=index_CNTL
indices4(2,:)=index_totalwind

;------------------------------------------------------------
wks = gsn_open_wks("x11","sst")   ;Open a workstation.
;wks = gsn_open_wks("pdf","TuneModle_forced_totalwind_"+NINOindex) ; Open a workstation.
gsn_define_colormap(wks,"NCV_blue_red")
;gsn_define_colormap(wks,"grads_rainbow")
plot = new(3,graphic)
plot1 = new(1,graphic)

res                   = True                     ; plot mods desired
res@xyDashPatterns  = (/1,0/)                         ; Make curves all solid
res@xyLineThicknesses = (/4.0,2.0/)          ; make second line thicker
res@xyLineColors      = (/"blue","red"/)          ; change line color
res@gsnDraw   = False  ; do not draw the plot (plots will be paneled)
res@gsnFrame  = False  ; do not advance the frame (plots will be paneled)

res@trXMinF   = 1950   ; starting point along X axis
res@trXMaxF   = 2011   ; ending point along X-axis
res@trYMinF   = 22   ; starting point along X axis
res@trYMaxF   = 31   ; ending point along X-axis
res@vpWidthF  = 0.66   ; stretch the plot to be wider (in NDC units)
res@vpHeightF = 0.25   ; and not as tall 
res@vpXF      = 0.15   ; set the start point along the X-axis in NDC units

;res@tmXBMode  ="Explicit"
;res@tmXBValues=ispan(1950,2011,5)
;res@tmXBLabels=""+(/1950,55,1960,65,1970,75,1980,85,1990,95,2000,05,2010/)
;res@tmXBLabels(11)="05"
;res@tmYLMode        = "Explicit" ;change y axis mark intervals to 0.5
;res@tmYLValues= fspan(-3,4,15)
;res@tmYLLabels      = "" + res@tmYLValues

res@tmXMajorGrid                = True          ; implement x grid 
res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
res@tmYMajorGrid                = True          ; implement y grid
res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
res@gsnXRefLine                 = (/0,1997/)

res@pmLegendDisplayMode    = "Always"         ; turn on legend
res@pmLegendSide           = "Top"            ; Change location of 
res@pmLegendParallelPosF   = .35               ; move units right
res@pmLegendOrthogonalPosF = -0.6             ; move units down
res@pmLegendWidthF         = 0.06              ; Change width and
res@pmLegendHeightF        = 0.07              ; height of legend.
res@lgPerimOn              = False            ; turn off box around
res@lgLabelFontHeightF     = .02              ; label font height
res@tiMainString      = "NINO34 index"                    ; add title
;plot(1)  = gsn_csm_xy (wks,yyyyfrac,indices_ano,res) ; create plot


res@gsnRightString="cor="+escorc(indices1(0,0:47*12-1),indices1(1,0:47*12-1))
res@xyExplicitLegendLabels = (/"AnoWind(No Seasonal Cycle)","TotalwindTo1997"/); create explicit labels
res@pmLegendParallelPosF   = .35               ; move units right
res@pmLegendOrthogonalPosF = -0.7             ; move units down
plot(0)  = gsn_csm_xy (wks,yyyyfrac,indices1,res) ; create plot

res@gsnRightString=""
res@gsnXRefLine                 = (/1991,1997/)
res@xyExplicitLegendLabels = (/"TotalwindRerunFromAno1991","TotalwindTo1997"/); create explicit labels
res@pmLegendParallelPosF   = .35               ; move units right
res@pmLegendOrthogonalPosF = -0.4             ; move units down
delete(res@tiMainString)
plot(1)  = gsn_csm_xy (wks,yyyyfrac,indices2,res) ; create plot

res@gsnRightString="cor="+escorc(indices3(0,:),indices3(1,:))
res@xyExplicitLegendLabels = (/"Observation","Totalwind"/); create explicit labels
res@pmLegendParallelPosF   = .2               ; move units right
res@pmLegendOrthogonalPosF = -0.5             ; move units down
res@gsnXRefLine                 = (/0,0/)
plot(2)  = gsn_csm_xy (wks,yyyyfrac,indices3,res) ; create plot

resP             = True
resP@gsnMaximize = True
;resP@gsnDraw   = False  ; do not draw the plot (plots will be paneled)
;resP@gsnFrame  = False  ; do not advance the frame (plots will be paneled)
gsn_panel(wks,plot,(/3,1/),resP)

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
;filename1 ="./NINO3_DegC.62yr.forced_totalwind.bin"
;system("rm -f " + filename1)
;fbindirwrite(filename1,indices(1,:))
;system("ls -l " + filename1)
end
