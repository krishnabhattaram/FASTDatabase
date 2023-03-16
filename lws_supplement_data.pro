pro lws_supplement_data
  
  ;;;; Despin fields to spacecraft frame ;;;;
  ucla_mag_despin
  fa_fields_despin, V58, V12
  get_data, 'ORBIT', data=orbit
  ;;;;;;;;
  
  ;;; Open CDF File ;;;
  filename = '/mydisks/disk1/database_files/sanitizedRuns1-9/fast_orbit_' + strcompress(string(orbit.y[0]), /REMOVE_ALL) + '.cdf'
  cdf2tplot, filename
  ;;;;;;;;
  
  ;;; Define times at which data will be included in the database ;;;
  TIME_ERROR = 2.5            ; Permissible deviation between absolute time and recorded data time
  ;;;;;;
  
  ;;; Spectral transform energy density ;;;
  magdata=get_fa_fields('MagDC')                                ; Get the fgm spec for fast survey data
  _ = get_fa_fields('Mag3ac_S',/all,/calibrate,/store)          ; Get three-component AC magnetic field from search-coil magnetometer
  get_data,'Mag3ac_S',data=mag3ac
  get_data,'E_ALONG_V',data=eav
  
  dtbz = {x: magdata.time, y: magdata.time - shift(magdata.time, 1)}
  dtb3 = {x: mag3ac.x, y: mag3ac.x - shift(mag3ac.x, 1)}
  dte = {x: eav.x, y: eav.x - shift(eav.x, 1)}
  
  
  get_data, 'magz_sp', data=magz_sp
  get_data, 'mag3ac_sp', data=mag3ac_sp
  get_data, 'eav_sp',data=eav_sp
  
  dtbz_times = fit_to_times(dtbz, magz_sp.x, TIME_ERROR)
  dtb3_times = fit_to_times(dtb3, mag3ac_sp.x, TIME_ERROR)
  dte_times = fit_to_times(dte, eav_sp.x, TIME_ERROR)

  _ = store_tplot_cdf('magz_sp_tres', data=dtbz_times, catdesc = 'Time resolution of fluxgate magnetometer z-axis s/c spin axis magnetic field (s)')
  _ = store_tplot_cdf('mag3ac_sp_tres', data=dtb3_times, catdesc = 'Time resolution of z-axis search coil magnetic field (s)')
  _ = store_tplot_cdf('eav_sp_tres', data=dte_times, catdesc='Time resolution of V5-V8 electric field (s)')
  ;;;;;;
    
  
  filename = '/mydisks/disk1/database_files/sanitizedRunsSupplemental2-2/fast_orbit_' + strcompress(string(orbit.y[0]), /REMOVE_ALL)

  tplot2cdf, filename=filename, tvars=['fa_vel', 'fa_pos', 'lat', 'lng', 'alt', 'flat', 'flng', 'ilat', 'ilng', 'mlt', 'b_foot', 'b_model', 'magz_sp', 'magz_sp_tres', $
    'mag3ac_sp', 'mag3ac_sp_tres', 'eav_sp', 'eav_sp_tres', 'dspV58_sp', 'dspmag3ac_sp', 'e_avg', 'e_sdev', 'b_avg', 'tot_b_avg', 'tot_b_sdev', $
    'j_2d_fa_ees', 'je_2d_fa_ees', 'j_2d_fa_ies', 'je_2d_fa_ies', 'j_2d_fa_ies_outflow', 'je_2d_fa_ies_outflow', $
    'j_2d_fa_ies_inflow', 'je_2d_fa_ies_inflow'], /default_cdf_structure
    
return
end