
copy_VarMeta(T1,T2)
copy_VarMeta(T1,T3)
copy_VarMeta(T1,T4)

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
index_CNTL1=dim_avg_n(T1(:,{latS:latN},{lonL:lonR}),(/1,2/)) 
index_CNTL2=dim_avg_n(T2(:,{latS:latN},{lonL:lonR}),(/1,2/)) 
index_CNTL3=dim_avg_n(T3(:,{latS:latN},{lonL:lonR}),(/1,2/)) 
index_CNTL4=dim_avg_n(T4(:,{latS:latN},{lonL:lonR}),(/1,2/)) 


;;;------------------------------------------------------------
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
;index_obs=rmAnnCycle1D(index_obs)
;;------------------------------------------------------------


yyyymm   = yyyymm_time(1, nyear, "integer")
yyyyfrac = yyyymm_to_yyyyfrac(yyyymm,0)   ; create fractional years for 

indices = new((/4,nt/),"float")
indices(0,:)=index_CNTL1
indices(1,:)=index_CNTL2
indices(2,:)=index_CNTL3
indices(3,:)=index_CNTL4

indices!0="var"
indices!1="time"
indices&time=yyyyfrac

;indices_ano=dim_rmvmean_n_Wrap(indices,(/1/))

;------------------------------------------------------------
print("plot obs&model nino index")
;wks = gsn_open_wks("x11","sst")   ;Open a workstation.
wks = gsn_open_wks("pdf","NINO3_CPD_Tsub+1") ;Open a workstation.
;wks = gsn_open_wks("pdf","NINO3_CPD_Tsub-1") ;Open a workstation.
gsn_define_colormap(wks,"NCV_blue_red")
;gsn_define_colormap(wks,"grads_rainbow")
plot = new(4,graphic)

res                   = True                     ; plot mods desired
res@xyDashPattern  = 0                         ; Make curves all solid
res@xyLineThicknesses = (/1.0,2.0/)          ; make second line thicker
res@xyLineColors      = (/"blue","red"/)          ; change line color
res@gsnDraw   = False  ; do not draw the plot (plots will be paneled)
res@gsnFrame  = False  ; do not advance the frame (plots will be paneled)

res@vpWidthF  = 0.66   ; stretch the plot to be wider (in NDC units)
res@vpHeightF = 0.2   ; and not as tall 
res@vpXF      = 0.15   ; set the start point along the X-axis in NDC units

res@tmXBMode  ="Explicit"
res@tmXBValues=ispan(20,nyear,20)
;res@tmXBLabels=""+ispan(1950,2010,10)
res@tmXBLabels=""+ispan(20,nyear,20)
res@tmXBMinorValues=ispan(10,nyear,20)
;res@tmYLMode        = "Explicit" ;change y axis mark intervals to 0.5
;res@tmYLValues= fspan(-3,4,15)
;res@tmYLLabels      = "" + res@tmYLValues

res@tmXMajorGrid                = True          ; implement x grid 
res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
res@tmXMinorGrid                = True          ; implement x grid 
res@tmXMinorGridThicknessF      = 1.0           ; 2.0 is default
res@tmXMinorGridLineDashPattern = 2             ; select short dash lines
res@tmYMajorGrid                = True          ; implement y grid
res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
res@tmYMajorGridLineDashPattern = 2             ; select short dash lines

;res@tiMainString=NINOindex+" index anomaly in ZC coupled run"
res@gsnRightString="b1=1.25  b2=3.0";control
res@gsnLeftString="(a)"
res@tiYAxisString="~S~o~N~C"
res@trYMinF   = -2   ; starting point along X axis
res@trYMaxF   = 4    ; ending point along X-axis
res@trXMinF   = 1    ; starting point along X axis
res@trXMaxF   = nyear   ; ending point along X-axis
res@tiXAxisString="Year"
plot(0)=gsn_csm_xy(wks,yyyyfrac,indices(0,:),res) ; FCD


res@gsnLeftString="(b)"
res@gsnRightString="b1=1.0  b2=3.0" ;Tsub+1
;res@gsnRightString="b1=1.5  b2=3.0" ;Tsub-1
plot(1)=gsn_csm_xy(wks,yyyyfrac,indices(1,:),res) ; FCD
res@gsnLeftString="(c)"
res@gsnRightString="b1=1.25  b2=4.0";Tsub+1
;res@gsnRightString="b1=1.25  b2=2.0";Tsub-1
plot(2)=gsn_csm_xy(wks,yyyyfrac,indices(2,:),res) ; FCD
res@gsnLeftString="(d)"
res@gsnRightString="b1=1.0  b2=4.0";Tsub+1
;res@gsnRightString="b1=1.5  b2=2.0";Tsub-1
plot(3)=gsn_csm_xy(wks,yyyyfrac,indices(3,:),res) ; FCD

resP             = True
resP@gsnMaximize = True
resP@gsnPanelMainString  = "~Z90~       "+NINOindex+" index anomaly in ZC coupled run"
;resP@gsnPanelFigureStrings = (/"(a)","(b)"/)
;resP@gsnPanelFigureStringsJust="TopLeft"
gsn_panel(wks,plot(0:3),(/4,1/),resP)

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
