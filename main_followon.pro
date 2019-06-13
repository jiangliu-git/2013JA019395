pro main_followon
;;; get different ratios, but with event-by-event binning

thm_init

;; the event
computer = 'I:'
;computer = '/home/jliu'

@declare

listname = 'dfb_list_lead_tail_fgs.txt'

;;; all quantities
quantities_check = 'number'

events_all = load_list(listname, folder = listfolder)
;events_all = events_all(*,0:110)

t0 = time_double(events_all(0,*))

;; 0809 events only
i_0809 = where(t0 lt time_double('2010 1 1'))
events_all = events_all(*, i_0809)

;; revise the events to be t0+3min
t0 = time_double(events_all(0,*))
tmax = time_double(events_all(1,*))
events_all(0,*) = time_string(t0+30.)
events_all(1,*) = time_string(t0+180.)

;;;;; specify different choices
;;; choose how to select DFBs (normalize or not)
dfb_select_type = 'bz'
;dfb_select_type = 'bz_nor2quiet_lobe'

;;;;; specify run of locations
;locations = 'all'
locations = ['tail', 'midt', 'mide', 'earth']
mark1 = -9
mark2 = -12
mark3 = -16
mark4 = -20
;title_str = ['X<-15RE', '-9RE>X>-15RE', 'X>-9RE']

;; the datatype and event
minutes_load = 5
sm_pts_fgs = 3 ;; used for DFCS paper
sm_pts_fgl = 10
window_size = 0.5
pre_size = 0.5
minutes_select = 2.
c_secondary = 30.

;;; choose fgs or fgl
data_dfb = 'fgs'
bfolder = fgs_folder
sm_pts = sm_pts_fgs

case dfb_select_type of
'bz': cri_ddts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
'bz_nor2quiet_lobe': cri_ddts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]/40.
endcase

;;;;;;;; get the event location and seperate
if n_elements(locations) gt 1 then begin
	loc_suf = '_'+locations
	pos_list = event_status_instant(events_all, 'pos', vtrange = time_string(events_all(0:1,*)), datafolder=pos_folder)
	x = pos_list(2,*)/6371.
	
	i_earth = where(x gt mark1)
	i_mide = where((x lt mark1) and (x gt mark2))
	i_midt = where((x lt mark2) and (x gt mark3))
	i_tail = where((x lt mark3) and (x gt mark4))
	store_data, 'earth', data=i_earth
	store_data, 'mide', data=i_mide
	store_data, 'midt', data=i_midt
	store_data, 'tail', data=i_tail
endif else store_data, 'all', data=lindgen(n_elements(events_all(0,*)))


number_loc = dblarr(n_elements(locations), n_elements(cri_ddts))

for l = 0, n_elements(locations)-1 do begin
	;;; get the quiet time lobe field for DFB selection
	get_data, locations(l), data = i_use
	events = events_all(*, i_use)
	if strcmp(dfb_select_type, 'bz_nor2quiet_lobe') then h = scale_height(events, b_folder = fgs_folder, p_folder = Pall_folder, time_length = 1., pre_time = 2.5, Blobe = Blobe)
	;;; repeat for criteria
	for k = 0, n_elements(cri_ddts)-1 do begin
		cri_ddt = cri_ddts(k)
		n_events_used = 0
		empty_arr = dblarr(n_elements(events(0,*)))
		empty_arr(*) = !values.f_nan
		n_dfbs_in_bbf = empty_arr
		;;;; begin survey the events
		for i = 0, n_elements(events(0,*))-1 do begin
			del_data, 'th*'
			event = events(*,i)
			if strcmp(dfb_select_type, 'bz_nor2quiet_lobe') then nor_v = Blobe(i) else nor_v = 1
			time_start = time_double(event(0))
			time_end = time_double(event(1))
			sc = event(3)
			;;;; select DFB
			trange_select = [time_start-minutes_select*60.,time_end+minutes_select*60.] ;;; for DFB
			del_data, 'th'+sc+'_'+data_dfb+'_gsm_tclip'
			load_bin_data, trange = trange_select, probe = sc, /tclip, datatype = data_dfb, datafolder = bfolder 
			if tv_exist('th'+sc+'_'+data_dfb+'_gsm_tclip') then begin
				secs_bbf_this = time_end-time_start ;; must be inside this if
				dfblist = dfb_select_bottom3('th'+sc+'_'+data_dfb+'_gsm_tclip', db = cri_ddt, window_size = window_size, sm_points = sm_pts, pre_size=pre_size, c_secondary = c_secondary, normalize = nor_v)
				if strcmp(dfblist(0), 'no event') then begin
					n_dfb_in = 0
				endif else begin
					n_dfb_in = n_elements(dfblist(0,*))
				endelse
				n_dfbs_in_bbf(i) = n_dfb_in
			endif ;; if of b data existence
		endfor ;; for of i, events
		;;;;;;; save the median value for events.
		number_loc(l,k) = median(n_dfbs_in_bbf, /even)
	endfor ;; for of k, criteria
endfor ;; for of l, locations

print, number_loc

stop
end
