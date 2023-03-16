pro lws_stats_sanitized

;;;; Define times and frequencies at which data is presented in the database ;;;;
  PERIOD = 5                  ; Interval at which data is reported (generally spin-period = 5s)
  VOLTAGE_THRESHOLD = 1000    ; mV/m (V58/V12 threshold, where data beyond this value is considered NaN)
  TIME_ERROR = 2.5            ; Permissible deviation between absolute time and recorded data time
  
  ELF_START_freq = 1
  ELF_END_freq = 160
  ELF_INTERVAL_freq = 1.25
  VLF_START_freq = 64
  VLF_END_freq = 17800
  VLF_N_freq = 16
  
  elf_fixed_frequencies = findgen(fix((ELF_END_freq - ELF_START_freq)/ELF_INTERVAL_freq) + 1) * ELF_INTERVAL_freq + ELF_START_freq
  VLF_INTERVAL_LOG = (alog(VLF_END_freq) - alog(VLF_START_freq))/(VLF_N_freq - 1)
  vlf_fixed_frequencies = [32, exp(findgen(VLF_N_freq) * VLF_INTERVAL_LOG + alog(VLF_START_freq))]
  vlf_error = (vlf_fixed_frequencies/2)*((vlf_fixed_frequencies[1] - vlf_fixed_frequencies[0])/vlf_fixed_frequencies[1])
;;;;;;;;

;;;; Despin fields to spacecraft frame ;;;;
  ucla_mag_despin
  fa_fields_despin, V58, V12
  get_data, 'ORBIT', data=orbit  
;;;;;;;;

;;; Define times at which data will be included in the database ;;;
  start_time = min(orbit.x)
  end_time = max(orbit.x)
  all_times = (indgen(fix((end_time - start_time)/PERIOD)) * PERIOD) + start_time
;;;;;;

;;; Telemetry Data ;;;
  get_data, 'LAT', data=lat
  get_data, 'LNG', data=lng
  get_data, 'ALT', data=alt
  get_data, 'FLAT', data=flat
  get_data, 'FLNG', data=flng
  get_data, 'ILAT', data=ilat
  get_data, 'ILNG', data=ilng
  get_data, 'MLT', data=mlt
  get_data, 'fa_vel', data=fa_vel
  get_data, 'fa_pos', data=fa_pos
  get_data, 'BFOOT', data=bfoot
  get_data, 'B_model', data=bmodel
  
  lat_times = fit_to_times(lat, all_times, TIME_ERROR)
  lng_times = fit_to_times(lng, all_times, TIME_ERROR)
  alt_times = fit_to_times(alt, all_times, TIME_ERROR)
  flat_times = fit_to_times(flat, all_times, TIME_ERROR)
  flng_times = fit_to_times(flng, all_times, TIME_ERROR)
  ilat_times = fit_to_times(ilat, all_times, TIME_ERROR)
  ilng_times = fit_to_times(ilng, all_times, TIME_ERROR)
  mlt_times = fit_to_times(mlt, all_times, TIME_ERROR)
  fa_vel_times = fit_to_times(fa_vel, all_times, TIME_ERROR)
  fa_pos_times = fit_to_times(fa_pos, all_times, TIME_ERROR)
  bfoot_times = fit_to_times(bfoot, all_times, TIME_ERROR)
  bmodel_times = fit_to_times(bmodel, all_times, TIME_ERROR)
  
  _ = store_tplot_cdf('lat', data=lat_times, catdesc='Latitude (degrees)')
  _ = store_tplot_cdf('lng', data=lng_times, catdesc = 'Longitude (degrees)')
  _ = store_tplot_cdf('flat', data=flat_times, catdesc = 'Foot Latitude (degrees)')
  _ = store_tplot_cdf('flng', data=flng_times, catdesc = 'Foot Longitude (degrees)')
  _ = store_tplot_cdf('ilat', data=ilat_times, catdesc = 'Invariant Latitude (degrees)')
  _ = store_tplot_cdf('ilng', data=ilng_times, catdesc = 'Invariant Longitude (degrees)')
  _ = store_tplot_cdf('mlt', data=mlt_times, catdesc = 'Magnetic Local Time (hours)')
  _ = store_tplot_cdf('alt', data=alt_times, catdesc = 'Altitude (km)')
  _ = store_tplot_cdf('fa_vel', data=fa_vel_times, catdesc = 'Spacecraft Velocity (3-component) (km/s) in GEI')
  _ = store_tplot_cdf('fa_pos', data=fa_pos_times, catdesc = 'Spacecraft Position (3-component) (km) in GEI')
  _ = store_tplot_cdf('b_foot', data=bfoot_times, catdesc = 'Footpoint magnetic field (3-component) in Cartesian coordinates (nT)')
  _ = store_tplot_cdf('b_model', data=bmodel_times, catdesc = 'Magnetic field from model (3-component) in Cartesian Coordinates (nT)')
;;;;;;

;;; Spectral transform energy density ;;;
  magdata=get_fa_fields('MagDC')                                ; Get the fgm spec for fast survey data  
  _ = get_fa_fields('Mag3ac_S',/all,/calibrate,/store)          ; Get three-component AC magnetic field from search-coil magnetometer
  get_data,'Mag3ac_S',data=mag3ac
  get_data,'E_ALONG_V',data=eav
  get_data, 'E_NEAR_B',data=enb
  
  ; Remove points where electric field data is reported to be higher than 1 V/m in either component (heuristic to avoid instrumentation error)
  enb.y[where(abs(enb.y) ge VOLTAGE_THRESHOLD)] = !VALUES.F_NAN
  eav.y[where(abs(eav.y) ge VOLTAGE_THRESHOLD)] = !VALUES.F_NAN 
  
  eav_sp = get_spectral_transform(eav, elf_fixed_frequencies, freq_error=1e-2, min_max_freq=2.5, frequency_resolution=0.25)
  
  magz_spt_av = sp_window_averaged_times(get_spectral_transform({x: magdata.time, y: magdata.comp3}, elf_fixed_frequencies, freq_error=1e-2, min_max_freq=2.5, frequency_resolution=0.25), $
                                     all_times, PERIOD, time_error=TIME_ERROR)
  eav_spt_av = sp_window_averaged_times(eav_sp, all_times, PERIOD, time_error=TIME_ERROR)
  mag3ac_spt_av = sp_window_averaged_times(get_spectral_transform(mag3ac, elf_fixed_frequencies, freq_error=1e-2, min_max_freq=2.5, frequency_resolution=0.25), $
                                     all_times, PERIOD, time_error=TIME_ERROR)
                                     
  _ = store_tplot_cdf('magz-sp', data=magz_spt_av, catdesc = 'nT!U2!N/Hz, z-axis s/c spin axis spectral energy density from fluxgate magnetometer')
  _ = store_tplot_cdf('eav-sp', data=eav_spt_av, catdesc = '(mV/m)!U2!N/Hz, electric field spectral energy density along spacecraft trajectory')
  _ = store_tplot_cdf('mag3ac-sp', data=mag3ac_spt_av, catdesc = 'nT!U2!N/Hz, instrument (Mag3ac) z-axis search coil spectral energy density')
;;;;;;

;;; DSP Quantities ;;;              
  fa_fields_dsp,'DspADC_V5-V8', t1=orbit.x[0], t2=orbit.x[-1]
  fa_fields_dsp,'DspMag3ac', t1=orbit.x[0], t2=orbit.x[-1]
  
  freq=(1+findgen(512))*32768/1024.  ;; TODO: CHECK THIS FREQUENCY!!!!!!!!!!!!!!!!!!!!!!!
  
  get_data,'DSP_V5-V8',data=datv
  get_data,'DSP_Mag3a',data=datm
  
  dspv={x:datv.x,y:datv.y,v:freq,spec:1}  
  dspm3ac={x:datm.x,y:datm.y,v:freq,spec:1}            
  
  dspv_av = sp_window_averaged_times(fit_to_frequencies(dspv, vlf_fixed_frequencies, vlf_error), all_times, PERIOD, time_error=TIME_ERROR)
  dspv_av_cleaned = {x: dspv_av.x, y: clean_where_exceeds(dspv_av, {x: V58.time, y:V58.comp1}, VOLTAGE_THRESHOLD, PERIOD, /absolute), v: dspv_av.v}
  dspm3ac_av = sp_window_averaged_times(fit_to_frequencies(dspm3ac, vlf_fixed_frequencies, vlf_error), all_times, PERIOD, time_error=TIME_ERROR)
  
  _ = store_tplot_cdf('dspV58-sp', data=dspv_av_cleaned, catdesc = 'log((V/m)!U2!N/Hz), spin plane electric field spectral density')
  _ = store_tplot_cdf('dspmag3ac-sp', data=dspm3ac_av, catdesc = 'log(nT!U2!N/Hz), search coil z-axis spectral energy density')
  
  options,'dspV58-sp','spec',1
  options,'dspV58-sp','ylog',1
  options,'dspV58-sp','zlog',0
  ;options,'dspV58-sp','zrange',1
  ;options,'dspv-s','yrange',1
  ;options,'dspv-s','yrange',[1,10.0]
  options, 'dspV58-sp', 'x_no_interp', 0
  options, 'dspV58-sp', 'x_no_interp', 0
  options,'dspV58-sp','zrange',[-11,-1]
  options,'dspV58-sp','ztitle','log((V/m)!U2!N/Hz)'
  options,'dspV58-sp','panel_size',3
  options,'dspV58-sp','ytitle','DSPV. Freq (Hz)'
  ;tplot, 'dspV58-sp'
  
  eav_time_indices = where(eav_sp.x ge 844891025 and eav_sp.x le 844891200)
  eav_times = {x: eav_sp.x[eav_time_indices], y: alog10(eav_sp.y[eav_time_indices, *]*1e-6), v: eav_sp.v, spec:1}
  store_data, 'eav-sp', data=eav_times;{x: eav_sp.x, y: alog10(eav_sp.y*1e-6), v:eav_sp.v, spec:1}
  options,'eav-sp','spec',1
  options,'eav-sp','ylog',1
  options,'eav-sp','zlog',0
  options,'eav-sp','zrange',1
  options,'eav-sp','yrange',1
  options,'eav-sp','yrange',[1,16000.0]
  options, 'eav-sp', 'x_no_interp', 0
  options, 'eav-sp', 'x_no_interp', 0
  options,'eav-sp','zrange',[-11,-1]
  options,'eav-sp','ztitle','log((V/m)!U2!N/Hz)'
  options,'eav-sp','panel_size',3
  options,'eav-sp','ytitle','EAV. Freq (Hz)'
  tplot, 'eav-sp'
  
  good_eav = where(finite(eav_sp.y[*, -1]), counteav)
  eav_good = eav_sp.y[good_eav, *]
  eav_means_good = exp(mean(alog(eav_good), dimension=1, /nan))
  
  pdf = histogram(alog(eav_good[*, 24]*1e-6), nbins=100, locations=binvals, /NAN)
  myplot = barplot(binvals, 1.0*pdf/max(pdf), title='EAV Values at ' + str(eav_sp.v[24]) + ' Hz', xtitle='Value log((V/m)^2/Hz)', ytitle='Frequency', xrange=[-30,2])
  
  dspv_time_indices = where(dspv.x ge 844891025 and dspv.x le 844891200)
  dspv_times = {x: dspv.x[dspv_time_indices], y: dspv.y[dspv_time_indices, *], v: dspv.v}
  store_data, 'dspv-s', data=dspv_times;dspv

  options,'dspv-s','spec',1
  options,'dspv-s','ylog',1
  options,'dspv-s','zlog',0
  options,'dspv-s','zrange',1
  options,'dspv-s','yrange',1
  options,'dspv-s','yrange',[1,16000.0]
  options, 'dspv-s', 'x_no_interp', 0
  options, 'dspv-s', 'x_no_interp', 0
  options,'dspv-s','zrange',[-11,-1]
  options,'dspv-s','ztitle','log((V/m)!U2!N/Hz)'
  options,'dspv-s','panel_size',3
  options,'dspv-s','ytitle','DSPV. Freq (Hz)'
  tplot, 'dspv-s'
  
  p = plot(eav_times.v, 10^mean(eav_times.y, dimension=1, /nan), title = 'EAV SP', xtitle = 'Frequency (Hz)', ytitle = '(V/m)^2/Hz', /xlog, /ylog);, /overplot, color='red')
  p = plot(dspv_times.v, 10^mean(dspv_times.y, dimension=1, /nan), title = 'Electric Field Spectral Energy Density', xtitle = 'Frequency (Hz)', ytitle = '(V/m)^2/Hz', /xlog, /overplot, /ylog)
  
  ;freqs = exp(findgen(10) * (alog(16000) - alog(500))/(10 - 1) + alog(500))
  
  ;for i=0,9 do begin
 ind = 31;value_locate(dspv.v, ;freqs[i])
 pdf = histogram(dspv.y[*, ind], nbins=50, locations=binvals, /NAN)
 ;myplot = barplot(binvals, 1.0*pdf/max(pdf), title='DSPV Values at ' + str(dspv.v[ind]) + ' Hz', xtitle='Value log((V/m)^2/Hz)', ytitle='Frequency', xrange=[-30, 2])
 ;endfor
  
  ;p = plot(eav_sp.v, mean(eav_good*1e-6, dimension=1, /nan), title = 'EAV SP', xtitle = 'Frequency (Hz)', ytitle = '(V/m)^2/Hz', /xlog, /ylog)
  ;p = plot(dspv.v, mean(10^dspv.y, dimension=1, /nan), title = 'Electric Field Spectral Energy Density', xtitle = 'Frequency (Hz)', ytitle = '(V/m)^2/Hz', /xlog, /overplot, /ylog)
  p = plot(eav_sp.v, eav_means_good*1e-6, title = 'EAV SP', xtitle = 'Frequency (Hz)', ytitle = '(V/m)^2/Hz', /xlog, /ylog, /overplot, color='red')
  p = plot(dspv.v, 10^mean(dspv.y, dimension=1, /nan), title = 'Electric Field Spectral Energy Density', xtitle = 'Frequency (Hz)', ytitle = '(V/m)^2/Hz', /xlog, /overplot, /ylog, color='red')
  
  ;magz_sp = get_spectral_transform({x: magdata.time, y: magdata.comp3}, elf_fixed_frequencies, freq_error=1e-2, min_max_freq=2.5, frequency_resolution=0.25)
  ;magz_sp.y = alog10(magz_sp.y)
  ;mag3ac_sp = get_spectral_transform(mag3ac, elf_fixed_frequencies, freq_error=1e-2, min_max_freq=2.5, frequency_resolution=0.25)
  ;mag3ac_sp.y = alog10(mag3ac_sp.y)
  
  store_data, 'dspm3ac', data=dspm3ac
  ;store_data, 'magz_spt', data=magz_sp
  ;store_data, 'mag3ac_spt', data=mag3ac_sp

  options,'dspm3ac','spec',1
  options,'dspm3ac','ylog',1
  options,'dspm3ac','zlog',0
  ;options,'dspm3ac','yrange',1
  ;options,'dspm3ac','yrange',[1,1000.0]
  options,'dspm3ac','zrange',[-14,-10]
  options,'dspm3ac','ztitle','log(nT!U2!N/Hz)'
  options,'dspm3ac','panel_size',3
  options,'dspm3ac','ytitle','B. Freq (Hz)'
  tplot, 'dspm3ac'
  
  store_data, 'dspm3aca', data=dspm3ac_av
  ;store_data, 'magz_spt', data=magz_sp
  ;store_data, 'mag3ac_spt', data=mag3ac_sp

  options,'dspm3aca','spec',1
  options,'dspm3aca','ylog',1
  options,'dspm3aca','zlog',0
  ;options,'dspm3aca','yrange',1
  ;options,'dspm3aca','yrange',[1,1000.0]
  options,'dspm3aca','zrange',[-14,-10]
  options,'dspm3aca','ztitle','log(nT!U2!N/Hz)'
  options,'dspm3aca','panel_size',3
  options,'dspm3aca','ytitle','B. Freq (Hz)'
  tplot, 'dspm3aca'
  
  ;options,'magz_spt','spec',1
  ;options,'magz_spt','ylog',1
  ;options,'magz_spt','zlog',0
  ;options,'magz_spt','yrange',1
  ;options,'magz_spt','yrange',[0.1,10.0]
  ;options,'magz_spt','zrange',[-10,10]
  ;options,'magz_spt','ztitle','log(nT!U2!N2/Hz)'
  ;options,'magz_spt','panel_size',3
  ;options,'magz_spt','ytitle','B. Freq (Hz)'
  
  ;options,'mag3ac_spt','spec',1
  ;options,'mag3ac_spt','ylog',1
  ;options,'mag3ac_spt','zlog',0
  ;options,'magz_spt','yrange',1
  ;options,'magz_spt','yrange',[0.1,10.0]
  ;options,'magz_spt','zrange',[-10,10]
  ;options,'mag3ac_spt','ztitle','log(nT!U2!N2/Hz)'
  ;options,'mag3ac_spt','panel_size',3
  ;options,'mag3ac_spt','ytitle','B. Freq (Hz)'
  
  ;p = plot(magz_sp.v, mean(magz_sp.y, dimension=1, /nan), title = 'MagZ', xtitle = 'Frequency (Hz)', ytitle = '(nT)^2/Hz', /xlog)
  ;p = plot(mag3ac_sp.v, mean(mag3ac_sp.y, dimension=1, /nan), title = 'Mag3AC', xtitle = 'Frequency (Hz)', ytitle = '(nT)^2/Hz', /xlog)
  ;p = plot(dspm3ac.v, mean(dspm3ac.y, dimension=1, /nan), title = 'DSPM3AC', xtitle = 'Frequency (Hz)', ytitle = '(nT)^2/Hz', /xlog)
  
;;;;;;

;;; Averaged 3-Component Magnetic field ;;;
  get_data, 'dB_fac_v', data = db_fac_v
  
  b_avg_x = window_averaged_times({x: db_fac_v.x, y:db_fac_v.y[*, 0]}, all_times, PERIOD, time_error=TIME_ERROR)
  b_avg_y = window_averaged_times({x: db_fac_v.x, y:db_fac_v.y[*, 1]}, all_times, PERIOD, time_error=TIME_ERROR)
  b_avg_z = window_averaged_times({x: db_fac_v.x, y:db_fac_v.y[*, 2]}, all_times, PERIOD, time_error=TIME_ERROR)
  b_avgs = {x: b_avg_x.x, y: [[b_avg_x.y], [b_avg_y.y], [b_avg_z.y]]}
  
  _ = store_tplot_cdf('b_avg', data=b_avgs, catdesc = '3-component (nT) delta B (rel. to IGRF) along trajectory, delta B cross track, delta B field-aligned')
;;;;;;

;;; Averaged field values ;;;
  get_data,'B_gei',data=b_gei
  
  tot_b = {x: b_gei.x, y: sqrt(b_gei.y(*,0)*b_gei.y(*,0)+b_gei.y(*,1)*b_gei.y(*,1)+b_gei.y(*,2)*b_gei.y(*,2))}
  tot_e = {x: eav.x, y: sqrt(eav.y*eav.y + enb.y*enb.y)}
  
  tot_b_avg = window_averaged_times(tot_b, all_times, PERIOD, time_error=TIME_ERROR)
  e_avg = window_averaged_times(tot_e, all_times, PERIOD, time_error=TIME_ERROR)
  
  e_sdev = window_sdev_times(tot_e, all_times, PERIOD, time_error=TIME_ERROR)
  tot_b_sdev = window_sdev_times(tot_b, all_times, PERIOD, time_error=TIME_ERROR)
  
  _ = store_tplot_cdf('tot_b_avg', data=tot_b_avg, catdesc = '(nT) Spin Period Averaged Total Magnetic Field')
  _ = store_tplot_cdf('e_avg', data=e_avg, catdesc = '(mV/m) spin averaged electric field along spacecraft trajectory')
  _ = store_tplot_cdf('e_sdev', data=e_sdev, catdesc = '(mV/m) standard deviation in total electric field per spin interval')
  _ = store_tplot_cdf('tot_b_sdev', data=tot_b_sdev, catdesc = '(nT) standard deviation in total magnetic field per spin interval')
;;;;;;

;;; Particle Data ;;;
  get_2dt_ts,'j_2d','fa_ees',energy=[20.0,30000.]
  get_data,'j_2d_fa_ees',data=j_2d_ees
  j_2d_ees_avg = window_averaged_times(j_2d_ees, all_times, PERIOD, time_error=TIME_ERROR)

  get_2dt_ts,'je_2d','fa_ees',energy=[20.0,30000.]
  get_data,'je_2d_fa_ees',data=je_2d_ees
  je_2d_ees_avg = window_averaged_times(je_2d_ees, all_times, PERIOD, time_error=TIME_ERROR)
  
  get_2dt_ts_antispacecraftvel_integration,'j_2d','fa_ies', latitude=lat, energy=[20.0,1000.], /ANTISPACECRAFTVEL_INTEGRATION
  get_data,'j_2d_fa_ies',data=j_2d_ies
  j_2d_ies_avg = window_averaged_times(j_2d_ies, all_times, PERIOD, time_error=TIME_ERROR)
  
  get_2dt_ts_antispacecraftvel_integration,'je_2d','fa_ies', latitude=lat, energy=[20.0,1000.], /ANTISPACECRAFTVEL_INTEGRATION
  get_data,'je_2d_fa_ies',data=je_2d_ies
  je_2d_ies_avg = window_averaged_times(je_2d_ies, all_times, PERIOD, time_error=TIME_ERROR)
  
  get_2dt_ts_antispacecraftvel_integration,'j_2d','fa_ies', latitude=lat, energy=[20.0,1000.], /OUTFLOW
  get_data,'j_2d_fa_ies',data=j_2d_ies_out
  j_2d_ies_out_avg = window_averaged_times(j_2d_ies_out, all_times, PERIOD, time_error=TIME_ERROR)

  get_2dt_ts_antispacecraftvel_integration,'je_2d','fa_ies', latitude=lat, energy=[20.0,1000.], /OUTFLOW
  get_data,'je_2d_fa_ies',data=je_2d_ies_out
  je_2d_ies_out_avg = window_averaged_times(je_2d_ies_out, all_times, PERIOD, time_error=TIME_ERROR)
  
  get_2dt_ts_antispacecraftvel_integration,'j_2d','fa_ies', latitude=lat, energy=[20.0,1000.], /INFLOW
  get_data,'j_2d_fa_ies',data=j_2d_ies_in
  j_2d_ies_in_avg = window_averaged_times(j_2d_ies_in, all_times, PERIOD, time_error=TIME_ERROR)

  get_2dt_ts_antispacecraftvel_integration,'je_2d','fa_ies', latitude=lat, energy=[20.0,1000.], /INFLOW
  get_data,'je_2d_fa_ies',data=je_2d_ies_in
  je_2d_ies_in_avg = window_averaged_times(je_2d_ies_in, all_times, PERIOD, time_error=TIME_ERROR)
  
  _ = store_tplot_cdf('j_2d_fa_ees', data=j_2d_ees_avg, catdesc = 'Electron number flux integrated over pitch angles 0-360 (1/cm!U2!Ns) - 20eV to 30keV ')
  _ = store_tplot_cdf('je_2d_fa_ees', data=je_2d_ees_avg, catdesc = 'Electron energy flux integrated over pitch angles 0-360 (mW/m!U2!N) - 20eV to 30keV ')
  _ = store_tplot_cdf('j_2d_fa_ies', data=j_2d_ies_avg, catdesc = 'Ion number flux integrated opposite spacecraft velocity doubled (1/cm!U2!Ns) - 20eV to 1keV ')
  _ = store_tplot_cdf('je_2d_fa_ies', data=je_2d_ies_avg, catdesc = 'Ion energy flux integrated opposite spacecraft velocity doubled (mW/m!U2!N) - 20eV to 1keV')
  _ = store_tplot_cdf('j_2d_fa_ies_outflow', data=j_2d_ies_out_avg, catdesc = 'Ion number flux integrated in outflow quadrant opposite spacecraft velocity doubled (1/cm!U2!Ns) - 20eV to 1keV')
  _ = store_tplot_cdf('je_2d_fa_ies_outflow', data=je_2d_ies_out_avg, catdesc = 'Ion energy flux integrated in outflow quadrant opposite spacecraft velocity doubled (mW/m!U2!N) - 20eV to 1keV')
  _ = store_tplot_cdf('j_2d_fa_ies_inflow', data=j_2d_ies_in_avg, catdesc = 'Ion number flux integrated in inflow quadrant opposite spacecraft velocity doubled (1/cm!U2!Ns) - 20eV to 1keV')
  _ = store_tplot_cdf('je_2d_fa_ies_inflow', data=je_2d_ies_in_avg, catdesc = 'Ion energy flux integrated in inflow quadrant opposite spacecraft velocity doubled (mW/m!U2!N) - 20eV to 1keV')
;;;;;

;;; Store Data to CDF file ;;;  
  filename = '/mydisks/disk1/database_files/sanitizedRuns1-9/fast_orbit_' + strcompress(string(orbit.y[0]), /REMOVE_ALL)

  tplot2cdf, filename=filename, tvars=['fa_vel', 'fa_pos', 'lat', 'lng', 'alt', 'flat', 'flng', 'ilat', 'ilng', 'mlt', 'b_foot', 'b_model', 'magz-sp', $
                                       'mag3ac-sp', 'eav-sp', 'dspV58-sp', 'dspmag3ac-sp', 'e_avg', 'e_sdev', 'b_avg', 'tot_b_avg', 'tot_b_sdev', $ 
                                       'j_2d_fa_ees', 'je_2d_fa_ees', 'j_2d_fa_ies', 'je_2d_fa_ies', 'j_2d_fa_ies_outflow', 'je_2d_fa_ies_outflow', $
                                       'j_2d_fa_ies_inflow', 'je_2d_fa_ies_inflow'], /default_cdf_structure
;;;;;;

return 
end
