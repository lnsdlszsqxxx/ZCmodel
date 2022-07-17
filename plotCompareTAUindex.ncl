;plot input and output TAU index for verification
;liang, 2022.7.2

begin

yrStart=1
yrEnd=62

nx=34
ny=30
nyear=yrEnd-yrStart+1
nt=(yrEnd-yrStart+1)*12

U = new((/nt,ny,nx/),"float")
V = new((/nt,ny,nx/),"float")
H = new((/nt,ny,nx/),"float") ;H1 H2 in ZC
T = new((/nt,ny,nx/),"float") ;T is SST
TAU = new((/nt,2,ny,nx/),"float")
TAUx = new((/nt,ny,nx/),"float")
TAUy = new((/nt,ny,nx/),"float")

suffix=yrStart+"-"+yrEnd
histFnU="ZCHIST_SSTA_34x30_Deg.bin"
histFnV="ZCHIST_SSTA_34x30_Deg.bin"
histFnH="ZCHIST_SSTA_34x30_Deg.bin"
;histFnT="ZCHIST_SSTA_34x30_Deg.bin"
histFnT="ZCFCD_SSTA_34x30_Deg.bin"
histFnTAUZC="ZCFCD_TAUA_34x30_dyne.bin" ;ZC output
;histFnTAUGC="tauxy_ano_GC.bin" ;wind stress used for forcing GC model
histFnTAUGC="tauxy_ano_Tao.bin" ;wind stress used for forcing GC model
histFnTAUNCEP="tauxy_ano_NCEP.bin" ;wind stress used for forcing ZC model
histFnTAUxNCEP="~/calculateS/NCEPReanalysis/uflx.sfc.mon.mean.nc"
histFnTAUyNCEP="~/calculateS/NCEPReanalysis/vflx.sfc.mon.mean.nc";original wind stress I downloaded for ZC forcing


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
;TAU = fbindirread(histFnTAUZC,0,(/nt,2,ny,nx/),"float")
T=TAU(:,0,:,:)
T = where(T.eq.0, T@_FillValue, T)

setfileoption("bin","ReadByteOrder","BigEndian")
TAUGC = fbindirread(histFnTAUGC,0,(/nt,2,ny,nx/),"float")
TGC=TAUGC(:,0,:,:)

setfileoption("bin","ReadByteOrder","BigEndian")
TAUNCEP = fbindirread(histFnTAUNCEP,0,(/nt,2,ny,nx/),"float")
TNCEP=TAUNCEP(:,0,:,:)

f = addfile(histFnTAUxNCEP, "r")
TAUxNCEP = f->uflx ;894*94*192, N/m^2,1948.1-2022.6
;TAUxNCEP=dim_rmvmean_n_Wrap(TAUxNCEP,0) ;use rmMonAnnCycTLL instead of this
TAUxNCEP=TAUxNCEP*(-10) ;dy/cm2
TAUxN=TAUxNCEP(24:767,:,:) ;1950.1-2011.12, anom
TAUxN=rmMonAnnCycTLL(TAUxN)
;printVarSummary(TAUxN)

T!0="month"
T!1="lat"
T!2="lon"
T&lat=lat
T&lon=lon
T@lat2d=lat2d
T@lon2d=lon2d
T@long_name = "ZC output TAUx"
T@units     = "dynes/cm2"
TGC!0="month"
TGC!1="lat"
TGC!2="lon"
TGC&lat=lat
TGC&lon=lon
TGC@lat2d=lat2d
TGC@lon2d=lon2d
TGC@long_name = "GC forcing TAUx"
TGC@units     = "dynes/cm2"
TNCEP!0="month"
TNCEP!1="lat"
TNCEP!2="lon"
TNCEP&lat=lat
TNCEP&lon=lon
TNCEP@lat2d=lat2d
TNCEP@lon2d=lon2d
TNCEP@long_name = "NCEP forcing TAUx"
TNCEP@units     = "dynes/cm2"

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
index_TAUZC=dim_avg_n(T(:,{latS:latN},{lonL:lonR}),(/1,2/)) 
index_TAUGC=dim_avg_n(TGC(:,{latS:latN},{lonL:lonR}),(/1,2/)) 
index_TAUNCEP=dim_avg_n(TNCEP(:,{latS:latN},{lonL:lonR}),(/1,2/)) 
index_TAUxN=dim_avg_n(TAUxN(:,{latS:latN},{lonL:lonR}),(/1,2/)) 

;print(index_CNTL-index_TAUin)

;;------------------------------------------------------------
;print("read obs nino index")
;ncol=10
;nRowSkip=1
;data=readAsciiTable("ersst4.nino.mth.81-10.ascii.txt",ncol,"float",nRowSkip)
;if(NINOindex.eq."NINO34") then
;index_obs=data(0:nt-1,8)
;else if(NINOindex.eq."NINO3") then
;index_obs=data(0:nt-1,4)
;end if
;end if
;;;as the model has no seasonal cycle in it, we need to remove the seasonal cycle in obs as well.
;print("NOTE: remove obs seasonal cycle")
;;index_obs=rmAnnCycle1D(index_obs)+avg(index_obs)
;;------------------------------------------------------------


yyyymm   = yyyymm_time(1, nyear, "integer")
yyyyfrac = yyyymm_to_yyyyfrac(yyyymm,0)   ; create fractional years for 

indices = new((/4,nt/),"float")
indices(0,:)=index_TAUZC
indices(1,:)=index_TAUGC
indices(2,:)=index_TAUNCEP
indices(3,:)=index_TAUxN

indices!0="var"
indices!1="time"
indices&time=yyyyfrac

indices_ano=dim_rmvmean_n_Wrap(indices,(/1/))

;------------------------------------------------------------
print("plot obs&model nino index")
wks = gsn_open_wks("x11","sst")   ;Open a workstation.
;wks = gsn_open_wks("pdf","TAUA_index_62yr_GC_ZC") ;Open a workstation.
gsn_define_colormap(wks,"NCV_blue_red")
;gsn_define_colormap(wks,"grads_rainbow")
plot = new(3,graphic)

res                   = True                     ; plot mods desired
res@xyDashPatterns  = (/0,0/)                         ; Make curves all solid
res@xyLineThicknesses = (/2.0,2.0/)          ; make second line thicker
res@xyLineColors      = (/"blue","red"/)          ; change line color
res@gsnDraw   = False  ; do not draw the plot (plots will be paneled)
res@gsnFrame  = False  ; do not advance the frame (plots will be paneled)

res@vpWidthF  = 0.66   ; stretch the plot to be wider (in NDC units)
res@vpHeightF = 0.2   ; and not as tall 
res@vpXF      = 0.15   ; set the start point along the X-axis in NDC units

res@tmXBMode  ="Explicit"
res@tmXBValues=ispan(1,nyear,10)
res@tmXBLabels=""+ispan(1950,2010,10)
;res@tmXBLabels=""+ispan(10,nyear,10)
res@tmXBMinorValues=ispan(6,nyear,10)
;res@tmYLMode        = "Explicit" ;change y axis mark intervals to 0.5
;res@tmYLValues= fspan(-3,4,15)
;res@tmYLLabels      = "" + res@tmYLValues

res@tmXMajorGrid                = True          ; implement x grid 
res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
res@tmYMajorGrid                = True          ; implement y grid
res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
res@tmYMajorGridLineDashPattern = 2             ; select short dash lines


res@pmLegendDisplayMode    = "Always"         ; turn on legend
res@pmLegendSide           = "Top"            ; Change location of 
res@pmLegendParallelPosF   = .15               ; move units right
res@pmLegendOrthogonalPosF = -0.75             ; move units down
res@pmLegendWidthF         = 0.06              ; Change width and
res@pmLegendHeightF        = 0.05              ; height of legend.
res@lgPerimOn              = False            ; turn off box around
res@lgLabelFontHeightF     = .015              ; label font height
res@xyExplicitLegendLabels = (/"GC forcing","ZC forcing"/); create explicit     labels


res@tiMainString="Wind stress index in NINO3 region"
res@gsnRightString="dynes/cm2"
res@gsnLeftString="index anomaly"
res@trYMinF   = -0.5   ; starting point along X axis
res@trYMaxF   = 0.5    ; ending point along X-axis
res@trXMinF   = 1    ; starting point along X axis
res@trXMaxF   = 62   ; ending point along X-axis
plot(0)=gsn_csm_xy(wks,yyyyfrac,indices(1:2,:),res) ; create plot

resP             = True
resP@gsnMaximize = True
gsn_panel(wks,plot(0),(/1,1/),resP)

;**********************************
; now maximize the plot
;**********************************   
; Maximize these plots for PS or PDF output. If output is
; to an X11 window or an NCGM file, then this routine doesn't
; do anything.
;psres = True
;maximize_output(wks,psres)  ; calls draw and frame for you 

;;---------------------
;;save data
;;---------------------
;print("saving...")
;;---Remove file just in case
;filename1 ="./NINO3_coupled_90yr.bin"
;system("rm -f " + filename1)
;fbindirwrite(filename1,indices(1,:))
;system("ls -l " + filename1)

end
