pro main_ratio_indiv
;;; get different ratios, but with event-by-event binning

thm_init

;; the event
computer = 'I:'
;computer = '/home/jliu'

@declare

listname = 'bbf_ps_list_original2.txt'

;;; all quantities
;quantities_check = 'number time bflux Bttl Bz Vttl Vx Ey'
;quantities_check = 'Bttl Bz Vttl Vx Ey'
;quantities_check = 'number Vx Bz Ey' ;;; for paper use
;quantities_check = 'Vx Bz Ey' ;;; for paper use
;quantities_check = 'ni Pth' ;;; for paper use
;quantities_check = 'time bflux'
;quantities_check = 'number'
;quantities_check = 'Bz Bttl'
quantities_check = 'Vx'

;;; whether to use vexb for removeing electric field offset or 0 to remove
dtrd_suf = ''
;dtrd_suf = '_vexb'
if ~(strmatch(quantities_check, '*Ey*') or strmatch(quantities_check, '*bflux*')) then dtrd_suf= ''
if strcmp(dtrd_suf, '_vexb') then edtrd_folder = vexb_dsl_folder else edtrd_folder = ''

;;; choose whether refine the list regarding vperp
;vperp_ave_max = -1. ;; no requirement
vperp_ave_max = 50.
;vperp_std_max = -1. ;; no requirement
vperp_std_max = 30.

events_all = load_list(listname, folder = listfolder)
;events_all = events_all(*,100:110)
loc_suf = ''
title_str = 'All'

;;;;; specify different choices
;;; choose how to select DFBs (normalize or not)
dfb_select_type = 'bz'
;dfb_select_type = 'bz_nor2quiet_lobe'
;;; choose which flux to use for magnetic flux
flux_use = 'efs_dsl'
;flux_use = 'vixb'

;;;;; specify run of locations
locations = 'all'
;locations = ['tail', 'mid', 'earth']
;mark1 = -9
;mark2 = -15
;title_str = ['X<-15RE', '-9RE>X>-15RE', 'X>-9RE']

;;;;; limit to ratio so that results are not contaminated
ratio_limit = 1000.

;; the datatype and event
minutes_load = 5
sm_pts_fgs = 3 ;; used for DFCS paper
;sm_pts_fgs = 1 ;; used for more DFBs
sm_pts_fgl = 10
window_size = 0.5
pre_size = 0.5
minutes_select = 2.
c_secondary = 30.
;; the duration of each DFB
;dfb_length = 0.5 ;; in minutes
dfb_length = 2./3. ;; in minutes

;; filename
filename_pre = 'dfb_bbf_indiv_ratio_'
;; size of the figure
xsize = 4.5
ysize = 3

case dfb_select_type of
'bz': dfb_select_suf = ''
'bz_nor2quiet_lobe': dfb_select_suf = '_dfbnor2lobe'
endcase
case dfb_select_type of
'bz': cri_ddts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
'bz_nor2quiet_lobe': cri_ddts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]/40.
endcase

;;; for test
;cri_ddts = [0.4, 0.5]

;;; choose fgs or fgl
data_dfb = 'fgs'
bfolder = fgs_folder
sm_pts = sm_pts_fgs

;;;;;;;; limiting the events regarding quiet time requirement
lim_suf = ''
if vperp_ave_max gt 0 then begin
	lim_suf = '_limit'
	vperp_ave_list = event_status_instant(events_all, 'viperp', pre_time = 2.5, time_length=1, datafolder=viperp_folder)
	vperp_ave = vperp_ave_list(5,*)
	i_use = where(vperp_ave lt vperp_ave_max)
	events_all = events_all(*, i_use)
endif
if vperp_std_max gt 0 then begin
	lim_suf = '_limit'
	vperp_std_list = event_status_instant(events_all, 'viperp', pre_time = 2.5, time_length=1, datafolder=viperp_folder, /std)
	vperp_std = vperp_std_list(5,*)
	i_use = where(vperp_std lt vperp_std_max)
	events_all = events_all(*, i_use)
endif
print, 'Event used:'
print, n_elements(events_all(0,*))

;;;;;;;; get the event location and seperate
if n_elements(locations) gt 1 then begin
	loc_suf = '_'+locations
	pos_list = event_status_instant(events_all, 'pos', vtrange = time_string(events_all(0:1,*)), datafolder=pos_folder)
	x = pos_list(2,*)/6371.
	
	i_earth = where(x gt mark1)
	i_mid = where((x lt mark1) and (x gt mark2))
	i_tail = where(x lt mark2)
	store_data, 'earth', data=i_earth
	store_data, 'mid', data=i_mid
	store_data, 'tail', data=i_tail
endif else store_data, 'all', data=lindgen(n_elements(events_all(0,*)))


for l = 0, n_elements(locations)-1 do begin
	;;; get the quiet time lobe field for DFB selection
	get_data, locations(l), data = i_use
	events = events_all(*, i_use)
	if strcmp(dfb_select_type, 'bz_nor2quiet_lobe') then h = scale_height(events, b_folder = fgs_folder, p_folder = Pall_folder, time_length = 1., pre_time = 2.5, Blobe = Blobe)
	;;;;;;;;;;;; bbf time durations
	secs_bbf = time_double(events(1,*))-time_double(events(0,*))
	;;;;;;;;;;; bbf average values
	if strmatch(quantities_check, '*V*') then begin
		v_bbf_list = event_status_instant(events, 'vi', vtrange = bbf_ranges, datafolder=vi_folder)
		vx_bbf = v_bbf_list(2,*)
		vttl_bbf = v_bbf_list(5,*)
	endif
	if strmatch(quantities_check, '*B*') then begin
		b_bbf_list = event_status_instant(events, 'fgs', vtrange = bbf_ranges, datafolder=fgs_folder)
		bz_bbf = b_bbf_list(4,*)
		bttl_bbf = b_bbf_list(5,*)
	endif
	if strmatch(quantities_check, '*Ey*') then begin
		e_bbf_list = event_status_instant_efi(events, 'efs_dsl', vtrange = bbf_ranges, rpre_time = 2.5, rtime_length=1., e_folder=efs_folder, /reverse_bc, vexb_dsl_folder = edtrd_folder)
		ey_bbf = e_bbf_list(3,*)
		rtranges = [time_double(events(0,*))-180., time_double(events(0,*))-160.]
	endif
	if strmatch(quantities_check, '*ni*') then begin
		ni_bbf_list = event_status_instant(events, 'ni', vtrange = bbf_ranges, datafolder=ni_folder)
		ni_bbf = ni_bbf_list(2,*)
	endif
	if strmatch(quantities_check, '*Pth*') then begin
		Pth_bbf_list = event_status_instant(events, 'Pth', vtrange = bbf_ranges, datafolder=Pth_folder)
		Pth_bbf = Pth_bbf_list(4,*)
	endif
	;;;;;;;;;; prepare arrays ready
	empty_arr = dblarr(n_elements(cri_ddts), n_elements(events(0,*)))
	empty_arr(*) = !values.f_nan
	n_dfbs_in_bbf = empty_arr
	if strmatch(quantities_check, '*time*') then ratio_t = empty_arr
	if strmatch(quantities_check, '*Vttl*') then ratio_v = empty_arr
	if strmatch(quantities_check, '*Vx*') then ratio_vx = empty_arr
	if strmatch(quantities_check, '*Bttl*') then ratio_b = empty_arr
	if strmatch(quantities_check, '*Bz*') then ratio_bz = empty_arr
	if strmatch(quantities_check, '*Ey*') then ratio_ey = empty_arr
	if strmatch(quantities_check, '*ni*') then ratio_ni = empty_arr
	if strmatch(quantities_check, '*Pth*') then ratio_Pth = empty_arr
	if strmatch(quantities_check, '*bflux*') then ratio_f = empty_arr

	for k = 0, n_elements(cri_ddts)-1 do begin
		cri_ddt = cri_ddts(k)
		n_events_used = 0
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
					dfb_tranges = [time_start-180., time_start-90]
					n_dfb_in = 0
				endif else begin
					dfb_tranges = [time_double(dfblist(0,*)), time_double(dfblist(0,*))+60.*dfb_length]
					n_dfb_in = n_elements(dfblist(0,*))
				endelse
				n_dfbs_in_bbf(k,i) = n_dfb_in
	
				;;;;; calculate the flux of DFB and BBF
				if strmatch(quantities_check, '*time*') or strmatch(quantities_check, '*bflux*') then begin
					case flux_use of
					'vixb': bbf_flux_cal, [time_start, time_end], dfb_tranges, probe = sc, b_type = data_dfb, b_folder = bfolder, v_type = 'vi', v_folder = vi_folder, flux_bbf = flux_bbf_this, flux_dfb = flux_dfb_this, /x_only, secs_dfb = secs_dfb_this, percent_t_dfb = percent_t_dfb, dfb_tranges_actual = dfb_tranges_actual
					'efs_dsl': bbf_flux_cal, [time_start, time_end], dfb_tranges, probe = sc, e_folder = efs_folder, rt_pre = 2.5, rt_length = 1., flux_bbf = flux_bbf_this, flux_dfb = flux_dfb_this, secs_dfb = secs_dfb_this, percent_t_dfb = percent_t_dfb, dfb_tranges_actual = dfb_tranges_actual, vexb_dsl_folder = edtrd_folder
					endcase
					if strmatch(quantities_check, '*time*') then ratio_t(k,i) = secs_dfb_this/secs_bbf_this
					if strmatch(quantities_check, '*bflux*') then ratio_f(k,i) = flux_dfb_this/flux_bbf_this
				endif else begin
					;;;; need to get the actual dfb tranges anyway
					dfb_tranges_actual = ranges_compress(time_double(dfb_tranges), limit = [time_start, time_end])
				endelse

				;;;;; calculate the average v, b, e of dfbs
				if strmatch(quantities_check, '*V*') then begin
					v_dfb = value_tranges(dfb_tranges_actual, 'vi', probe = sc, datafolder=vi_folder, /less)
					vx_dfb = v_dfb(0)
					vttl_dfb = v_dfb(3)
					if strmatch(quantities_check, '*Vttl*') then ratio_v(k,i) = vttl_dfb/vttl_bbf(i)
					if strmatch(quantities_check, '*Vx*') then ratio_vx(k,i) = vx_dfb/vx_bbf(i)
				endif
				if strmatch(quantities_check, '*B*') then begin
					b_dfb = value_tranges(dfb_tranges_actual, 'fgs', probe = sc, datafolder=fgs_folder, /less)
					bz_dfb = b_dfb(2)
					bttl_dfb = b_dfb(3)
					if strmatch(quantities_check, '*Bz*') then ratio_bz(k,i) = bz_dfb/bz_bbf(i)
					if strmatch(quantities_check, '*Bttl*') then ratio_b(k,i) = bttl_dfb/bttl_bbf(i)
				endif
				if strmatch(quantities_check, '*Ey*') then begin
					e_dfb = value_tranges(dfb_tranges_actual, 'efs_dsl', probe = sc, rtrange = rtranges(*,i), datafolder=efs_folder, /reverse_bc, /less, vexb_dsl_folder = edtrd_folder)
					ey_dfb = e_dfb(1)
					ratio_ey(k,i) = ey_dfb/ey_bbf(i)
				endif
				if strmatch(quantities_check, '*ni*') then begin
					ni_dfb = value_tranges(dfb_tranges_actual, 'ni', probe = sc, datafolder=ni_folder, /less)
					ratio_ni(k,i) = ni_dfb/ni_bbf(i)
				endif
				if strmatch(quantities_check, '*Pth*') then begin
					Pth_dfb_all = value_tranges(dfb_tranges_actual, 'Pth', probe = sc, datafolder=Pth_folder, /less)
					Pth_dfb = Pth_dfb_all(2)
					ratio_Pth(k,i) = Pth_dfb/Pth_bbf(i)
				endif
			endif ;; if of b data existence
		endfor ;; for of i, events
	endfor ;; for of k, the trys for ddt

	;;;;;;; set nan values: negative ratios and ratios exceeding limit will be set NaN
	if strmatch(quantities_check, '*bflux*') then begin
		i_bad_flux = where(ratio_f lt 0, j_bad_flux)
		if j_bad_flux gt 0 then ratio_f(i_bad_flux) = !values.f_nan
		i_big_flux = where(ratio_f gt ratio_limit, j_big_flux)
		if j_big_flux gt 0 then begin
			ratio_f(i_big_flux) = !values.f_nan
			ind_asnd = sort(ratio_f)
			ratio_f(ind_asnd(0:j_big_flux-1)) = !values.f_nan
		endif
	endif
	if strmatch(quantities_check, '*Vx*') then begin
		i_bad_vx = where(ratio_vx lt 0, j_bad_vx)
		if j_bad_vx gt 0 then ratio_vx(i_bad_vx) = !values.f_nan
		i_big_vx = where(ratio_vx gt ratio_limit, j_big_vx)
		if j_big_vx gt 0 then begin
			ratio_vx(i_big_vx) = !values.f_nan
			ind_asnd = sort(ratio_vx)
			ratio_vx(ind_asnd(0:j_big_vx-1)) = !values.f_nan
		endif
	endif
	if strmatch(quantities_check, '*Ey*') then begin
		i_bad_ey = where(ratio_ey lt 0, j_bad_ey)
		if j_bad_ey gt 0 then ratio_ey(i_bad_ey) = !values.f_nan
		i_big_ey = where(ratio_ey gt ratio_limit, j_big_ey)
		if j_big_ey gt 0 then begin
			ratio_ey(i_big_ey) = !values.f_nan
			ind_asnd = sort(ratio_ey)
			ratio_ey(ind_asnd(0:j_big_ey-1)) = !values.f_nan
		endif
	endif
	if strmatch(quantities_check, '*Bz*') then begin
		i_bad_bz = where(ratio_bz lt 0, j_bad_bz)
		if j_bad_bz gt 0 then ratio_bz(i_bad_bz) = !values.f_nan
		i_big_bz = where(ratio_bz gt ratio_limit, j_big_bz)
		if j_big_bz gt 0 then begin
			ratio_bz(i_big_bz) = !values.f_nan
			ind_asnd = sort(ratio_bz)
			ratio_bz(ind_asnd(0:j_big_bz-1)) = !values.f_nan
		endif
	endif
	
	;;;;;;; plot the data
	case dfb_select_type of
	'bz': xrange = [0., 1.]
	'bz_nor2quiet_lobe': xrange = [0., 1.]/40.
	endcase
	case dfb_select_type of
	'bz': xunit = 'nT/s'
	'bz_nor2quiet_lobe': xunit = '/s'
	endcase

	xtitle = 'cri_ddt ['+xunit+']'
	
	;;;;;;; plot the ratios
	if strmatch(quantities_check, '*number*') then begin
		dataout_simple, save_folder+'/number'+loc_suf+lim_suf, [[cri_ddts], [double(n_dfbs_in_bbf)]]
		popen, pic_folder+'/'+filename_pre+'number'+loc_suf+dfb_select_suf
		print_options,xsize=xsize,ysize=ysize
		stat_plot_disc, cri_ddts, n_dfbs_in_bbf, qtt_range = xrange, qtt_2_range = [0., 20.], qtt_title = xtitle, qtt_2_title = '# of DFBs in each BBF', title = 'DFBs in BBFs';+title_str(l)+' DFB:'+dfb_select_suf
		xyouts, 0.9, 0.9*!y.crange(1), '(a)'
		pclose
	endif
	if strmatch(quantities_check, '*time*') then begin
		popen, pic_folder+'/'+filename_pre+'t'+loc_suf+dfb_select_suf
		print_options,xsize=xsize,ysize=ysize
		stat_plot_disc, cri_ddts, ratio_t, qtt_range = xrange, qtt_2_range = [0., 1.], qtt_title = xtitle, qtt_2_title = 'ratio of time';, title = 'DFB/FB Ratio '+title_str(l)+' DFB:'+dfb_select_suf
		xyouts, 0.9, 0.9*!y.crange(1), '(e)'
		pclose
	endif
	if strmatch(quantities_check, '*bflux*') then begin
		popen, pic_folder+'/'+filename_pre+'flux'+loc_suf+dfb_select_suf
		print_options,xsize=xsize,ysize=ysize
		stat_plot_disc, cri_ddts, ratio_f, qtt_range = xrange, qtt_2_range = [0., 1.], qtt_title = xtitle, qtt_2_title = 'ratio of flux';, title = 'DFB/FB Ratio '+title_str(l)+' DFB:'+dfb_select_suf
		xyouts, 0.9, 0.9*!y.crange(1), '(f)'
		pclose
	endif
	if strmatch(quantities_check, '*Vttl*') then begin
		popen, pic_folder+'/'+filename_pre+'vttl'+loc_suf+dfb_select_suf
		print_options,xsize=xsize,ysize=ysize
		stat_plot_disc, cri_ddts, ratio_v, qtt_range = xrange, qtt_2_range = [0., 3.], qtt_title = xtitle, qtt_2_title = 'ratio of <Vi_ttl>';, title = 'DFB/FB Ratio '+title_str(l)+' DFB:'+dfb_select_suf;, qtt_2_range = 0
		xyouts, 0.9, 0.9*!y.crange(1), '(c)'
		pclose
	endif
	if strmatch(quantities_check, '*Vx*') then begin
		dataout_simple, save_folder+'/ratio_vx'+loc_suf+lim_suf, [[cri_ddts], [ratio_vx]]
		popen, pic_folder+'/'+filename_pre+'vx'+loc_suf+dfb_select_suf
		print_options,xsize=xsize,ysize=ysize
		stat_plot_disc, cri_ddts, ratio_vx, qtt_range = xrange, qtt_2_range = [0., 6.], qtt_title = xtitle, qtt_2_title = 'ratio of <Vi_x>';, title = 'DFB/FB Ratio '+title_str(l)+' DFB:'+dfb_select_suf;, qtt_2_range = 0
		xyouts, 0.9, 0.9*!y.crange(1), '(c)'
		pclose
	endif
	if strmatch(quantities_check, '*Bttl*') then begin
		popen, pic_folder+'/'+filename_pre+'bttl'+loc_suf+dfb_select_suf
		print_options,xsize=xsize,ysize=ysize
		stat_plot_disc, cri_ddts, ratio_b, qtt_range = xrange, qtt_2_range = [0., 3.], qtt_title = xtitle, qtt_2_title = 'ratio of <B_ttl>', title = 'DFB/BBF Ratio';+title_str(l)+' DFB:'+dfb_select_suf;, qtt_2_range = 0
		xyouts, 0.9, 0.9*!y.crange(1), '(b)'
		pclose
	endif
	if strmatch(quantities_check, '*Bz*') then begin
		dataout_simple, save_folder+'/ratio_bz'+loc_suf+lim_suf, [[cri_ddts], [ratio_bz]]
		popen, pic_folder+'/'+filename_pre+'bz'+loc_suf+dfb_select_suf
		print_options,xsize=xsize,ysize=ysize
		stat_plot_disc, cri_ddts, ratio_bz, qtt_range = xrange, qtt_2_range = [0., 6.], qtt_title = xtitle, qtt_2_title = 'ratio of <Bz>', title = 'DFB/FB Ratio';+title_str(l)+' DFB:'+dfb_select_suf;, qtt_2_range = 0
		xyouts, 0.9, 0.9*!y.crange(1), '(b)'
		pclose
	endif
	if strmatch(quantities_check, '*Ey*') then begin
		dataout_simple, save_folder+'/ratio_ey'+dtrd_suf+loc_suf+lim_suf, [[cri_ddts], [ratio_ey]]
		popen, pic_folder+'/'+filename_pre+'ey'+dtrd_suf+loc_suf+dfb_select_suf
		print_options,xsize=xsize,ysize=ysize
		stat_plot_disc, cri_ddts, ratio_ey, qtt_range = xrange, qtt_2_range = [0., 10.], qtt_title = xtitle, qtt_2_title = 'ratio of <Ey_dsl>';, title = 'DFB/FB Ratio '+title_str(l)+' DFB:'+dfb_select_suf;, qtt_2_range = 0
		xyouts, 0.9, 0.9*!y.crange(1), '(d)'
		pclose
	endif
	if strmatch(quantities_check, '*ni*') then begin
		dataout_simple, save_folder+'/ratio_ni'+loc_suf+lim_suf, [[cri_ddts], [ratio_ni]]
		popen, pic_folder+'/'+filename_pre+'ni'+loc_suf+dfb_select_suf
		print_options,xsize=xsize,ysize=ysize
		stat_plot_disc, cri_ddts, ratio_ni, qtt_range = xrange, qtt_2_range = [0., 2.], qtt_title = xtitle, qtt_2_title = 'ratio of <ni>';, title = 'DFB/FB Ratio '+title_str(l)+' DFB:'+dfb_select_suf;, qtt_2_range = 0
		xyouts, 0.9, 0.9*!y.crange(1), '(d)'
		pclose
	endif
	if strmatch(quantities_check, '*Pth*') then begin
		dataout_simple, save_folder+'/ratio_Pth'+loc_suf+lim_suf, [[cri_ddts], [ratio_Pth]]
		popen, pic_folder+'/'+filename_pre+'Pth'+loc_suf+dfb_select_suf
		print_options,xsize=xsize,ysize=ysize
		stat_plot_disc, cri_ddts, ratio_Pth, qtt_range = xrange, qtt_2_range = [0., 2.], qtt_title = xtitle, qtt_2_title = 'ratio of <Pth>';, title = 'DFB/FB Ratio '+title_str(l)+' DFB:'+dfb_select_suf;, qtt_2_range = 0
		xyouts, 0.9, 0.9*!y.crange(1), '(d)'
		pclose
	endif
endfor ;; for of l, locations

stop
end
