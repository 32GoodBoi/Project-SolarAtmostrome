pro example_aiasyn_hv,day=day,hour=hour,min=mins,$
  min_snr=min_snr,sat_lvl=sat_lvl,serr_per=serr_per,$
  tempdir=tempdir,respdir=respdir, use_diw=use_diw,em=em

  ;  ;   About X3, Impulsive start of X-flare
;    if (n_elements(day) ne 1) then day='24-Oct-2014'
;    if (n_elements(hour) ne 1) then hour='21'
;    if (n_elements(mins) ne 1) then mins='26'

  ; About B7, few  ARs, day of C-flares
    if (n_elements(day) ne 1) then day='25-Dec-2015'
    if (n_elements(hour) ne 1) then hour='12'
    if (n_elements(mins) ne 1) then mins='00'

  ;  ; About A5, quiet day
  ;  if (n_elements(day) ne 1) then day='01-Jul-2010'
  ;  if (n_elements(hour) ne 1) then hour='12'
  ;  if (n_elements(mins) ne 1) then mins='00'

  date=day+' '+hour+':'+mins+':00'


  if (n_elements(min_snr) lt 1) then min_snr=3.0

  if (n_elements(sat_lvl) lt 1) then sat_lvl=1.5e4

  if (n_elements(serr_per) lt 1) then serr_per=0.0


  if (n_elements(temp_dir) ne 1) then temp_dir=curdir()+'/temp/'


  if (n_elements(resp_dir) ne 1) then resp_dir=curdir()+'/resp/'

 
  ystr=strmid(anytim(date,/ccsds),0,4)
  mstr=strmid(anytim(date,/ccsds),5,2)
  dstr=strmid(anytim(date,/ccsds),8,2)
  durl='http://jsoc.stanford.edu/data/aia/synoptic/'+ystr+'/'+mstr+'/'+dstr+'/H'+hour+'00/'
  wavenum=['94','131','171','193','211','335']
  ffs='AIA'+ystr+mstr+dstr+'_'+hour+mins+'_0'+['094','131','171','193','211','335']+'.fits'
  furls=durl+ffs

  ftest=file_test(temp_dir+ffs)
  if (product(ftest) ne 1) then wgetc=ssw_wget_mirror(furls,temp_dir,/spawn)


  read_sdo,temp_dir+ffs,ind,data,/uncomp_delete

  wv_srt=sort(ind.wavelnth)
  ind=ind[wv_srt]
  data=data[*,*,wv_srt]

  nx=n_elements(data[*,0,0])
  ny=n_elements(data[0,*,0])
  nf=n_elements(data[0,0,*])

  xc=ind[0].xcen
  yc=ind[0].ycen
  dx=ind[0].cdelt1
  dy=ind[0].cdelt2

  id=where(data ge sat_lvl,nid)
  if (nid gt 1) then data[id]=0.0


  id=where(data le 0,nid)
  if (nid gt 1) then data[id]=0.0

  edg0=10
  data[0:edg0-1,*,*]=0.0
  data[*,0:edg0-1,*]=0.0
  data[nx-edg0:nx-1,*,*]=0.0
  data[*,nx-edg0:nx-1,*]=0.0


  edata=fltarr(nx,ny,nf)

  npix=4096.^2/(nx*ny)
  edata=fltarr(nx,ny,nf)
  gains=[18.3,17.6,17.7,18.3,18.3,17.6]
  dn2ph=gains*[94,131,171,193,211,335]/3397.
  rdnse=1.15*sqrt(npix)/npix
  drknse=0.17
  qntnse=0.288819*sqrt(npix)/npix
  ; error in DN/px
  for i=0, nf-1 do begin
    etemp=sqrt(rdnse^2.+drknse^2.+qntnse^2.+(dn2ph[i]*abs(data[*,*,i]))/(npix*dn2ph[i]^2))
    esys=serr_per*data[*,*,i]/100.
    edata[*,*,i]=sqrt(etemp^2. + esys^2.)
  endfor
  id=where(data/edata le min_snr,nid)
  if (nid gt 1) then data[id]=0.0


  durs=ind.exptime
  for i=0, nf-1 do data[*,*,i]=data[*,*,i]/durs[i]
  for i=0, nf-1 do edata[*,*,i]=edata[*,*,i]/durs[i]

  ; What temperature binning do you want for the DEM?
  ; temps variable are the bin edges
  ; These are the bin edges
  ;  temps=[0.5,1,1.5,2,3,4,6,8,11,14,19,25,32]*1e6
  ;  temps=[0.5,1,2,4,6,8,11,14,19]*1d6
  ;  logtemps=alog10(temps)


  logtemps=5.7+findgen(17)*0.1
  temps=10d^logtemps


  mlogt=get_edges(logtemps,/mean)
  nt=n_elements(mlogt)
  respfile='resp/aia_resp_'+mstr+ystr+'.dat'
  if (file_test(respfile) eq 0) then begin
    tresp=aia_get_response(/temperature,/dn,/chianti,/noblend,/evenorm,timedepend_date=ystr+'/'+mstr+'/01')
    save,file=respfile,tresp
  endif
  restore,file=respfile

  idc=[0,1,2,3,4,6]

  tr_logt=tresp.logte
  gdt=where(tr_logt ge min(logtemps) and tr_logt le max(logtemps),ngd)
  tr_logt=tr_logt[gdt]
  TRmatrix=tresp.all[*,idc]
  TRmatrix=TRmatrix[gdt,*]
  if keyword_set(use_diw) then begin

    dem_norm0=dblarr(nx,ny,nt)
    dem_norm_temp=exp(-(findgen(nt)+1-nt*0.5)^2/12)
    dem_norm_temp=smooth(dem_norm_temp,3)

    for xx=0,nx-1 do begin
      for yy=0,ny-1 do begin
        dem_norm0[xx,yy,*]=dem_norm_temp
      endfor

      dn2dem_pos_nb, data,edata,TRmatrix,tr_logt,temps,dem,edem,elogt,chisq,dn_reg,/timed,dem_norm0=dem_norm0
    endfor
  endif else begin
    
    dn2dem_pos_nb, data,edata,TRmatrix,tr_logt,temps,dem,edem,elogt,chisq,dn_reg,/timed
  endelse


  dem_maps=replicate(make_map(fltarr(nx,ny),xc=xc,yc=yc,dx=dx,dy=dy,date=date),nt)

 
  if keyword_set(em) then begin
    dt=fltarr(nt)
    for i=0,nt-1 do dt[i]=10d^(logtemps[i+1])-10d^(logtemps[i])
    em=dblarr(nx,ny,nt)
    for xx=0,nx-1 do begin
      for yy=0,ny-1 do begin
        em[xx,yy,*]=dem[xx,yy,*]*dt
      endfor
    endfor
    dem_maps.data=em
    idd='EM '
  endif else begin
    dem_maps.data=dem
    idd='DEM '
  endelse

  temp_ids=strarr(nt)
  for i=0,nt-1 do temp_ids[i]=string(logtemps[i],format='(f3.1)')+' - '+string(logtemps[i+1],format='(f3.1)')
  dem_maps.id=idd+'logT: '+temp_ids


  !p.multi=[0,4,4]
  loadct,5,/silent
  
  if keyword_set(em) then drang=[1d25,1d30] else drang=[1d19,1d23]
  
  for i=0,nt-1 do plot_map,dem_maps[i],title=dem_maps[i].id,/log,dmin=drang[0],dmax=drang[1],chars=2.0
  xyouts,10,10,date,chars=1.2,/device

  stop

  stop
end
