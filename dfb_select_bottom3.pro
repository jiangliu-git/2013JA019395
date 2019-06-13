function dfb_select_bottom3, tplot_variable, window_size = window_size, pre_size = pre_size, sm_points = sm_points, db = db, increase = increase, strictness =strictness, c_secondary = c_secondary, btotal = btotal, normalize = normalize, Blobe_tv = Blobe_tv
;;;; same as select_bottom2, but tuned especially for bbf_flux
;; The faster version, using derivatives. control the bottom level of selecting dfbs: applying criterion to windows.
;; tplot_variable: the tplot variable to work on, should be a fgl or fgs with three components.
;; window_size: in minutes
;; sm_points: points to smooth the data
;; db: criteria for the derivative (for the smoothed data)
;; increase: the increase of next beginning of the window, in seconds
;; pre_size: in minutes, the previous time to the window, if not set, equals to window size.
;; strictness: 1 to 4
;; c_secondary: the criteria in seconds for a secondary event
;; btotal: use btotal instead of Bz to apply criteria
;; normalize: normalize the field to a value (scaler), e. g., lobe field, the db criterion should be changed too.
;; Blobe_tv: a tplot variable containing only lobe field, if set, do a point by point normalize the field to lobe field, the db criterion should be changed too.
;; return variable is 'event_list' (array): [peak dB/dt time (string), peak B time (string), type], sorted by time, is the time of peak Bz.
;; Jiang Liu, 3/21/2011

event_list = ['no event','no event', 'no type'] ; [peak db/dt time, peak bz time, type] type m for major, type s for secondary
;; if tplot_variable does not exist then return
tplot_names, tplot_variable, names = tv_exist
if ~tv_exist(tplot_variable) then return, event_list

if not keyword_set(sm_points) then sm_points = 3
if not keyword_set(db) then begin
	if ~keyword_set(normalize) and ~keyword_set(Blobe_tv) then db = 0.5 else db = 0.5/40.
endif
if not keyword_set(normalize) then normalize = 1.
if not keyword_set(window_size) then window_size = 0.5
if not keyword_set(pre_size) then pre_size = window_size 
if not keyword_set(increase) then increase = 4.
if not keyword_set(strictness) then strictness = 3
if not keyword_set(c_secondary) then c_secondary = 30.

get_data, tplot_variable, data = b_data
; first, 3-point smooth the data
if keyword_set(btotal) then begin
	strength = sqrt(total(b_data.y^2, 2))
	store_data, tplot_variable+'_ttl', data = {x:b_data.x, y:strength}
	ddt_variable = tplot_variable+'_ttl'
endif else begin
	split_vec, tplot_variable
	ddt_variable = tplot_variable+'_z'
endelse
tsmooth2, ddt_variable, sm_points
deriv_data, ddt_variable+'_sm'
get_data, ddt_variable+'_sm_ddt', data = ddt_data
good_ddt = where(ddt_data.y/normalize gt db)

if good_ddt(0) ne -1 then begin ; there are big enough derivatives
  if n_elements(good_ddt) gt 1 then begin
    good_shadow = [good_ddt(0)-1, good_ddt(0:n_elements(good_ddt)-2)]
    find_i = good_ddt - good_shadow
    i_crack = where(find_i gt 1)
  endif else i_crack = -1
  if i_crack(0) ne -1 then i_crack = [i_crack, n_elements(good_ddt)-1] ; must consider the last term of good
  for i = 0, n_elements(i_crack)-1 do begin ; for of cracks
	  if i_crack(i) eq -1 then begin
  		i_tgs = 0
  		i_tge = n_elements(good_ddt)-1
  	endif else begin
  		if i eq 0 then i_tgs = 0 else i_tgs = i_crack(i-1) ; tgs: time, good, start
  		if i eq n_elements(i_crack)-1 then i_tge = i_crack(i) else i_tge = i_crack(i)-1 ; tge: time, good, end
  	endelse
  	i_ts = good_ddt(i_tgs)
  	i_te = good_ddt(i_tge)
  	examine_trange = [ddt_data.x(i_ts)-window_size*60., ddt_data.x(i_te)+window_size*60]
;  	;;;;;;; strictness = 1: begin the sliding window within examine_trange ;;;;;;;;;;;;;;;;;;;;;;;;
;  	ts = examine_trange(0)
;    while ts+window_size*60 lt examine_trange(1) do begin
;      window_trange = [ts, ts+window_size*60]
;      pre_trange = [ts-pre_size*60, ts]
;      good_w = where((b_data.x ge window_trange(0)) and (b_data.x le window_trange(1)))
;      good_p = where((b_data.x ge pre_trange(0)) and (b_data.x le pre_trange(1)))
;      yes = 0
;      if (good_p(0) ne -1) and (good_w(0) ne -1) then begin
;        window_b_data = b_data.y(good_w, *)
;        pre_b_data = b_data.y(good_p, *)
;        yes = dfb_criteria(pre_b_data, window_b_data)
;      endif
;      if yes eq 1 then begin
;        max_nouse = max(ddt_data.y(i_ts:i_te), i_ddt)
;        time_temp = ddt_data.x(i_ts:i_te)
;        event_list = [event_list, time_temp(i_ddt)]
;        break
;      endif else ts = ts+increase
;    endwhile
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; strictness >= 2: use the peak of the derivative ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	time_ddt = ddt_data.x(good_ddt(i_tgs):good_ddt(i_tge))
	;; option 1: use the beginning of large db/dt as t_ddt
	t_ddt = time_ddt(0)
	;; option 2: use the max of dB/dt as t_ddt
	;max_nouse = max(ddt_data.y(good_ddt(i_tgs):good_ddt(i_tge)), i_ddt) ; the index (in ddt_data) of the max derivative point of each qualified derivative range 
	;t_ddt = time_ddt(i_ddt) ; time of max dB/dt
	;;;;;;;;;;;; determine the window
	window_trange = [t_ddt-0.0*window_size*60., t_ddt+1.0*window_size*60.] ; adjust window trange here
  good_w = where((b_data.x ge window_trange(0)) and (b_data.x le window_trange(1)))
  if good_w(0) eq -1 then continue ; no data in the window, continue
  window_b_data = b_data.y(good_w, *)
  time_b = b_data.x(good_w)
  if keyword_set(btotal) then begin
	max_nouse = max(sqrt(total(window_b_data^2, 2)), i_maxb)
  endif else begin
	max_nouse = max(window_b_data(*,2), i_maxb)
  endelse
  t_maxb = time_b(i_maxb) ; time of max Bz
  
	;;;;;;;;;;;; if it's a secondary event I don't bother test it with the criteria ;;;;;;;
	n_caught = n_elements(event_list(0,*))
	if n_caught gt 1 then begin
;		if t_ddt-time_double(event_list(1,n_caught-1)) lt c_secondary then begin ; time difference from t_ddt of the current event to the max Bz of the previous event
		if t_maxb-time_double(event_list(1,n_caught-1)) lt c_secondary then begin ; time difference of max Bz of this event to previous event
			event_list = [[event_list], [time_string(t_ddt), time_string(t_maxb), 's']]
			continue
		endif
	endif
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; strictness = 2: window is fixed, pre varies    ;;;;;;;;;;;;;;;;;;;;;;;;;;;
	if strictness eq 2 then begin
		pre_trange_all = [window_trange(0)-4*pre_size*60., window_trange(0)]
		pre_ts = pre_trange_all(0)
		while pre_ts+pre_size*60. le pre_trange_all(1) do begin
	  		pre_trange = [pre_ts, pre_ts+pre_size*60.]
			good_p = where((b_data.x ge pre_trange(0)) and (b_data.x le pre_trange(1)))
			if (good_p(0) ne -1) then begin
  				pre_b_data = b_data.y(good_p, *)
  				if dfb_criteria(pre_b_data, window_b_data) then begin
					event_list = [[event_list], [time_string(t_ddt),  time_string(t_maxb), 'm']]
					break ; if catch, no need to continue
				endif
			endif
			pre_ts = pre_ts+pre_size*60.
		endwhile
	endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; strictness =3: both window and pre are fixed ranges    ;;;;;;;;;;;;;;;;;;;;;;;;;;;
	if strictness eq 3 then begin
  		pre_trange = [window_trange(0)-pre_size*60., window_trange(0)]
  		good_p = where((b_data.x ge pre_trange(0)) and (b_data.x le pre_trange(1)))
  		if (good_p(0) ne -1) then begin
  			window_b_data = b_data.y(good_w, *)
  			pre_b_data = b_data.y(good_p, *)
			;;; ignore other dfb criteria
  			if dfb_criteria(pre_b_data, window_b_data, dbz = -100., m_dbz = -100., mm_dbz = -100.) then event_list = [[event_list], [time_string(t_ddt),  time_string(t_maxb), 'm']]
  		endif
	endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  endfor ; for of cracks (ranges of large dB/dt)
  ; no good_ddt, no events
endif else print, 'DFB_SELECT_BOTTOM2: '+tplot_variable+' do not have big enough derivative!'

if n_elements(event_list(0,*)) le 1 then print, 'DFB_SELECT_BOTTOM2: No event found.' else event_list = event_list(*,1:*)
return, event_list
end
