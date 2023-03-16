;+
; :Author: Krishnakumar Bhattaram
; 
; Contains helper functions which are used to generate the FAST database (CDF format)
; from SDT. This project is part of the Living with a Star program, and is conducted
; under the supervision of Christopher Chaston (SSL) and Ennio Sanchez (Stanford)
;-

function value_locate_error,times,x,error

  closest_indices = []
  foreach time, times do begin
    closest_indices = [closest_indices, min(where(x - time eq min(x - time, /ABSOLUTE)))]
  endforeach

  closest_times = x(closest_indices)

  dt = abs(closest_times - times)
  bad_indices = where(dt gt error, count)

  if (count gt 0) then closest_indices[bad_indices] = -1
  return, closest_indices

end

function float_generate,y,closest_indices
  bad_indices = where(closest_indices lt 0, count)
  rdata = y[closest_indices, *]

  if (count gt 0) then rdata[bad_indices, *] = !VALUES.F_NAN

  return, rdata
end

function fit_to_times, data, times, error
  time_indices = value_locate_error(times, data.x, error)
  return, {x: float_generate(data.x, time_indices), y: float_generate(data.y, time_indices)}
end

function fit_to_frequencies, data, frequencies, error
  freq_indices = value_locate_error(frequencies, data.v, error)
  return, {x: data.x, y: transpose(float_generate(transpose(data.y), freq_indices)), v: float_generate(data.v, freq_indices)}
end

;+
; :Description:
;    Describe the procedure.
;
; :Params:
;    data:
;    freqs
;
; :Keywords:
;    min_max_freq
;    frequency_resolution
;    freq_error
;    overlap
;    name
;
; :Author: krishnabhattaram
;-
function get_spectral_transform, data, freqs, min_max_freq=min_freq, frequency_resolution=min_delta_freq, freq_error=freq_error, overlap=overlap, name=name

  if not keyword_set(min_delta_freq) then min_delta_freq = 1 ; 1 Hz as minimum frequency delta (required frequency resolution)
  if not keyword_set(min_freq) then min_freq = 5 ; Minimum temporal resolution (converted to hertz by Nyquist)
  if not keyword_set(freq_error) then freq_errror = min_delta_freq * 1e-2 ; Permissible error in frequency bin
  if not keyword_set(overlap) then overlap = 0.5 ; FFT interval overlaps
  if not keyword_set(name) then name = 'dy_sp' ; Name of data quantity

  dt = data.x - shift(data.x, 1)
  startind = 0
  endind = 0

  total_data_times = []
  total_data_values = []

  while startind le n_elements(dt) - 1 and startind ge 0 do begin
    interval_dt = dt[startind + 1]
    endind = min(where(dt[startind + 1:*] ne interval_dt)) + startind
    if endind le startind then endind = n_elements(dt) - 1

    interval_max_freq = 0.5/interval_dt ; Nyquist limited highest frequency for this interval
    if (interval_max_freq lt min_freq) then begin
      startind = endind + 1
      continue
    endif

    n_pts_resolution = fix(2*interval_max_freq/min_delta_freq) ; Number of points required for FFT to have sufficient resolution
    if (endind - startind lt n_pts_resolution) then begin
      startind = endind + 1
      continue
    endif

    data_interval = {x: data.x[startind:endind], y: data.y[startind:endind]}
    clu_pow_spec_fixed,data_interval,av,n_pts_resolution,overlap

    get_data,'dy_sp',data=interval_spectral_transformed
    freq_vals_struct = fit_to_frequencies(interval_spectral_transformed, freqs, freq_error)
    freq_vals = freq_vals_struct.y

    total_data_times = [total_data_times, interval_spectral_transformed.x]
    total_data_values = [total_data_values, freq_vals]

    startind = endind + 1
  endwhile

  if keyword_set(name) then store_data, name, data={x: total_data_times, y: total_data_values, v: freqs}
  return, {x: total_data_times, y: total_data_values, v: freqs}
end

function clean_where_exceeds, data, measurement, threshold, period, absolute=absolute
  if keyword_set(absolute) then bad_times = measurement.x[where(abs(measurement.y) ge threshold)] $
  else bad_times = measurement.x[where(measurement.y ge threshold)]
  period_indices = value_locate_error(bad_times, data.x, period/2)

  ret_y = data.y
  ret_y[period_indices(where(period_indices ge 0)), *] = !VALUES.F_NAN

  return, ret_y
end

function window_average, data, indices, period, d

  avgs = []
  
  foreach index, indices do begin
    
    if (index lt 0) then begin
      if d eq 1 then avg = make_array(n_elements(data.y[index, *]), value=!VALUES.F_NAN) else avg = !VALUES.F_NAN
      avgs = [[avgs], [avg]]
      continue
    endif 
      
    if index eq 0 then begin
      dt = data.x[index+1] - data.x[index]
    endif else begin
      dt = data.x[index] - data.x[index-1]
    endelse
    
    if (dt le 0) then begin
      if d eq 1 then avg = make_array(n_elements(data.y[index, *]), value=!VALUES.F_NAN) else avg = !VALUES.F_NAN
      avgs = [[avgs], [avg]]
      continue
    endif

    points_per_period = fix(period/dt) + 1
    start_i = max([0, index - fix(points_per_period/2)])
    end_i = min([index + fix(points_per_period/2), n_elements(data.x) - 1])
    
    if (data.x[end_i] - data.x[start_i] ge 1.5*period) then begin
      if d eq 1 then avg = make_array(n_elements(data.y[index, *]), value=!VALUES.F_NAN) else avg = !VALUES.F_NAN
      avgs = [[avgs], [avg]]
      continue
    endif

    avg = mean(data.y[start_i:end_i, *], DIMENSION=d)
    avgs = [[avgs], [avg]]
    
  endforeach

  return, transpose(avgs)
end

function window_sdev, data, indices, period, d

  sdevs = []

  foreach index, indices do begin

    if (index lt 0) then begin
      if d eq 1 then sdev = make_array(n_elements(data.y[index, *]), value=!VALUES.F_NAN) else sdev = !VALUES.F_NAN
      sdevs = [[sdevs], [sdev]]
      continue
    endif

    if index eq 0 then begin
      dt = data.x[index+1] - data.x[index]
    endif else begin
      dt = data.x[index] - data.x[index-1]
    endelse

    if (dt le 0) then begin
      if d eq 1 then sdev = make_array(n_elements(data.y[index, *]), value=!VALUES.F_NAN) else sdev = !VALUES.F_NAN
      sdevs = [[sdevs], [sdev]]
      continue
    endif

    points_per_period = fix(period/dt) + 1
    start_i = max([0, index - fix(points_per_period/2)])
    end_i = min([index + fix(points_per_period/2), n_elements(data.x) - 1])

    if (data.x[end_i] - data.x[start_i] ge 1.5*period) then begin
      if d eq 1 then sdev = make_array(n_elements(data.y[index, *]), value=!VALUES.F_NAN) else sdev = !VALUES.F_NAN
      sdevs = [[sdevs], [sdev]]
      continue
    endif

    sdev = stddev(data.y[start_i:end_i, *], DIMENSION=d)
    sdevs = [[sdevs], [sdev]]

  endforeach

  return, transpose(sdevs)
end

function window_averaged_times, data, times, period, time_error=time_error, dim=dim
  if not keyword_set(time_error) then time_error = period * 5e-1
  if not keyword_set(dim) then dim = 0

  time_indices = value_locate_error(times, data.x, time_error)
  retdata = {x: float_generate(data.x, time_indices), y: window_average(data, time_indices, period, dim)}
  valid_time_indices = where(finite(retdata.x))

  return, {x: retdata.x[valid_time_indices], y: retdata.y[valid_time_indices, *]}
end

function window_sdev_times, data, times, period, time_error=time_error, dim=dim
  if not keyword_set(time_error) then time_error = period * 5e-1
  if not keyword_set(dim) then dim = 0

  time_indices = value_locate_error(times, data.x, time_error)
  retdata = {x: float_generate(data.x, time_indices), y: window_sdev(data, time_indices, period, dim)}
  valid_time_indices = where(finite(retdata.x))

  return, {x: retdata.x[valid_time_indices], y: retdata.y[valid_time_indices, *]}
end

function sp_window_averaged_times, data, times, period, time_error=time_error
  window_averaged = window_averaged_times({x: data.x, y: data.y}, times, period, time_error=time_error, dim=1)
  return, {x: window_averaged.x, y: window_averaged.y, v: data.v}
end

function get_dlimit, data, name, catdesc=catdesc
  if not keyword_set(catdesc) then catdesc = 'Undefined'

  newline = string(10B)
  tab = string(9B)

  rules_of_use = " In order to confirm the reliability of the data, you are requested to contact the PI, Dr. Christopher Chaston, before using any data in oral/poster presentations."
  rules_of_use += " In the publication, you should contact to the PI before the submission for the confirmation of the principle of authorship/acknowledgement." + newline
  rules_of_use += " PI: " + newline + tab + "Dr. Christopher Chaston (Scientist)" + newline + tab + "Space Science Laboratories" + newline + tab + "Email: cchaston@berkeley.edu" + newline + newline
  rules_of_use += " PI: " + newline + tab + "Krishnakumar Bhattaram" + newline + tab + "University of California, Berkeley" + newline + tab + "Email: krishnabhattaram@berkeley.edu" + newline + newline

  gatt_general = {DESCRIPTOR: 'Undefined', LOGICAL_SOURCE: 'Undefined', LOGICAL_SOURCE_DESCRIPTION: 'Undefined', $
    GENERATION_DATE: str(SYSTIME()), START_TIME: str(systime(elapsed=min(data.x))),  END_TIME: str(systime(elapsed=max(data.x))), TITLE: name, $
    ACKNOWLEDGEMENT:'FAST data were sourced from SDT', ADID_REF: 'Undefined', DATA_TYPE: catdesc, $
    DISCIPLINE:  'Space Physics > Magnetospheric Science', GENERATED_BY: 'SSL', $
    INSTRUMENT_TYPE: 'Radio and Plasma Waves (space)', MISSION_GROUP: 'FAST', MODS: 'Undefined', $
    PI_AFFILIATION: 'SSL', PI_NAME: 'K. Bhattaram and C. Chaston', RULES_OF_USE: rules_of_use, SOURCE_NAME: 'FAST', $
    TEXT: 'Undefined', TIME_RESOLUTION: '5s'}

  vatt_general = {CATDESC: catdesc, DICT_KEY: 'Undefined', FIELDNAM: 'Undefined', FILLVAL: !VALUES.F_NAN, LABLAXIS: name, $
    FORMAT: '%3d',  UNITS: 'Undefined', VALIDMAX: 'Undefined', VALIDMIN: 'Undefined', VAR_NOTE: 'Undefined', $
    VAR_TYPE: 'undefined', DISPLAY_TYPE: 'undefined'}

  dlimit_general = {gatt: gatt_general, vatt: vatt_general}
  return, dlimit_general
end

function store_tplot_cdf, name, data=data, catdesc=catdesc
  if not keyword_set(data) then print, ('No data is being stored for: ' + name)
  
  del_data, name
  
  store_data, name, data=data
  options, name, 'cdf', get_dlimit(data, name, catdesc=catdesc), /default
  
end

