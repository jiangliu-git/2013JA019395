pro main_ratio_acum
;;; see how much flux transport is contributed by DFBs in BBFs
;;; select DFBs on the fly

thm_init

;; the event
computer = 'I:'
;computer = '/home/jliu'

@declare

listname = 'bbf_ps_list_original2'

;;;;;; select what quantities to look at
quantities_check = 'time bflux'
;quantities_check = 'Bttl Vttl Ey'
;quantities_check = 'Bz Vx Ey'

;;; whether to use vexb for removeing electric field offset or 0 to remove
dtrd_suf = ''
;dtrd_suf = '_vexb'
if ~(strmatch(quantities_check, '*Ey*') or strmatch(quantities_check, '*bflux*')) then dtrd_suf= ''
if strcmp(dtrd_suf, '_vexb') then edtrd_folder = vexb_dsl_folder else edtrd_folder = ''

;;;;; get the absolute value or not
abs_suf = ''
;abs_suf = '_abs'

events_all = load_list(listname+'.txt', folder = listfolder)
;;;; for test
;events_all = events_all(*,100:110)

;;;;; specify different choices
;;; if plot events
plot_events = 'no'
;plot_events = 'yes'

;;; if compare with Flow bursts
compare_fb = 'no'
;compare_fb = 'yes'

;;; choose how to select DFBs (normalize or not)
;dfb_select_type = 'bz'
dfb_select_type = 'bz_nor2quiet_lobe'

;;; choose which flux to use for magnetic flux
flux_use = 'efs_dsl'
;flux_use = 'vixb'

;;; choose whether refine the list regarding vperp
vperp_ave_max = -1. ;; no requirement
;vperp_ave_max = 50.
vperp_std_max = -1. ;; no requirement
;vperp_std_max = 30.

;;; choose and how to select FBs (use which data)
if strcmp(compare_fb, 'yes') then begin
	fb_select_type = 'vi'
	;fb_select_type = 'vixb'
	;fb_select_type = 'efs_dsl'
	;;;; the suffix
	fb_select_suf = '_'+fb_select_type
endif else begin
	fb_select_type = ''
	fb_select_suf = ''
endelse

;;;;; specify run of locations
;locations = 'all'
locations = ['tail', 'mid', 'earth']
mark1 = -9.
mark2 = -15.
title_str = ['X<-15RE', '-9RE>X>-15RE', 'X>-9RE']

;;;;;;;;;; run
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
dfb_length = 2./3. ;; in minutes
;dfb_length = 1. ;; in minutes

;; file name
filename_pre = 'dfb_bbf_acum_ratio_'

;; size of the figure
xsize = 4.5
ysize = 3

case dfb_select_type of
'bz': dfb_select_suf = ''
'bz_nor2quiet_lobe': dfb_select_suf = '_dfbnor2lobe'
endcase
case dfb_select_type of
;'bz': cri_ddts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
;'bz_nor2quiet_lobe': cri_ddts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]/40.
'bz': cri_ddts = [0.5]
'bz_nor2quiet_lobe': cri_ddts = [0.5]/40.
endcase

;;; choose fgs or fgl
data_dfb = 'fgs'
bfolder = fgs_folder
sm_pts = sm_pts_fgs

ratio_secs = dblarr(n_elements(cri_ddts))
ratio_flux = dblarr(n_elements(cri_ddts))
ratio_bz = dblarr(n_elements(cri_ddts))
ratio_bttl = dblarr(n_elements(cri_ddts))
ratio_vx = dblarr(n_elements(cri_ddts))
ratio_vttl = dblarr(n_elements(cri_ddts))
ratio_ey = dblarr(n_elements(cri_ddts))
ratio_inter_dfb = dblarr(n_elements(cri_ddts))
ratio_inter_fb = dblarr(n_elements(cri_ddts))

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
pos_list = event_status_instant(events_all, 'pos', vtrange = time_string(events_all(0:1,*)), datafolder=pos_folder)
x = pos_list(2,*)/6371.

if n_elements(locations) gt 1 then begin
	loc_suf = '_'+locations
	i_earth = where(x gt mark1)
	i_mid = where((x lt mark1) and (x gt mark2))
	i_tail = where(x lt mark2)
	store_data, 'earth', data=i_earth
	store_data, 'mid', data=i_mid
	store_data, 'tail', data=i_tail
endif else begin
	store_data, 'all', data=lindgen(n_elements(events_all(0,*)))
	loc_suf = ''
	title_str = ''
endelse

;temp = !p.multi
;!p.multi = [0,1,3,0,0]
qtt_bbf_3 = [0., 0., 0., 0.]
for l = 0, n_elements(locations)-1 do begin
;for l = 0, 0 do begin
	;;; get the quiet time lobe field for DFB selection
	get_data, locations(l), data = i_use
	events = events_all(*, i_use)
	if strcmp(dfb_select_type, 'bz_nor2quiet_lobe') then h = scale_height(events, b_folder = fgs_folder, p_folder = Pall_folder, time_length = 1., pre_time = 2.5, Blobe = Blobe)
	;;;;;;;;;;; bbf average values
	if strmatch(quantities_check, '*V*') then v_bbf = qtt_bbf_3
	if strmatch(quantities_check, '*B*') then b_bbf = qtt_bbf_3
	if strmatch(quantities_check, '*E*') then begin
		e_bbf = qtt_bbf_3
		rtranges = [time_double(events(0,*))-180., time_double(events(0,*))-160.]
	endif
	for i = 0, n_elements(events(0,*))-1 do begin
		sc = events(3,i)
		if strmatch(quantities_check, '*V*') then begin
			no_use = value_tranges(events(0:1,i), 'vi', probe = sc, datafolder=vi_folder, /less, combined_data = v_this)
			if strcmp(abs_suf, '_abs') then v_this = abs(v_this)
			v_bbf = [[v_bbf], [transpose(v_this)]]
		endif
		if strmatch(quantities_check, '*B*') then begin
			no_use = value_tranges(events(0:1,i), 'fgs', probe = sc, datafolder=fgs_folder, /less, combined_data = b_this)
			if strcmp(abs_suf, '_abs') then b_this = abs(b_this)
			b_bbf = [[b_bbf], [transpose(b_this)]]
		endif
		if strmatch(quantities_check, '*Ey*') then begin
			no_use = value_tranges(events(0:1,i), 'efs_dsl', probe = sc, rtrange = rtranges(*,i), datafolder=efs_folder, /reverse_bc, /less, combined_data = e_this, vexb_dsl_folder = edtrd_folder)
			if strcmp(abs_suf, '_abs') then e_this = abs(e_this)
			e_bbf = [[e_bbf], [transpose(e_this)]]
		endif
	endfor
	;;;; get the average values
	if strmatch(quantities_check, '*V*') then begin
		v_bbf = v_bbf(*, 1:*)
		v_bbf_x = v_bbf(0,*)
		v_bbf_ttl = v_bbf(3,*)
		v_bbf_x_ave = mean(v_bbf_x, /nan)
		v_bbf_ttl_ave = mean(v_bbf_ttl, /nan)
	endif
	if strmatch(quantities_check, '*B*') then begin
		b_bbf = b_bbf(*, 1:*)
		b_bbf_z = b_bbf(2,*)
		b_bbf_ttl = b_bbf(3,*)
		b_bbf_z_ave = mean(b_bbf_z, /nan)
		b_bbf_ttl_ave = mean(b_bbf_ttl, /nan)
	endif
	if strmatch(quantities_check, '*E*') then begin
		e_bbf = e_bbf(*, 1:*)
		e_bbf_y = e_bbf(1,*)
		e_bbf_ttl = e_bbf(3,*)
		e_bbf_y_ave = mean(e_bbf_y, /nan)
		e_bbf_ttl_ave = mean(e_bbf_ttl, /nan)
	endif
	;;;;;;;;;; begin surveying
	for k = 0, n_elements(cri_ddts)-1 do begin
		cri_ddt = cri_ddts(k)
		secs_bbf = 0.
		secs_dfb = 0.
		secs_fb = 0.
		secs_overlap = 0.
		flux_bbf = 0.
		flux_dfb = 0.
		n_events_used = 0
		;;;;;;;;;;; dfb average values
		if strmatch(quantities_check, '*V*') then v_dfb = qtt_bbf_3
		if strmatch(quantities_check, '*B*') then b_dfb = qtt_bbf_3
		if strmatch(quantities_check, '*E*') then e_dfb = qtt_bbf_3
		;;;; begin survey the events
		for i = 0, n_elements(events(0,*))-1 do begin
			del_data, 'th*'
			event = events(*,i)
			if strcmp(dfb_select_type, 'bz_nor2quiet_lobe') then nor_v = Blobe(i) else nor_v = 1
			time_start = time_double(event(0))
			time_end = time_double(event(1))
			sc = event(3)
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
	
				;;;;; calculate the flux of DFB and BBF
				if strmatch(quantities_check, '*time*') or strmatch(quantities_check, '*bflux*') then begin
					case flux_use of
					'vixb': bbf_flux_cal, [time_start, time_end], dfb_tranges, probe = sc, b_type = data_dfb, b_folder = bfolder, v_type = 'vi', v_folder = vi_folder, flux_bbf = flux_bbf_this, flux_dfb = flux_dfb_this, /x_only, secs_dfb = secs_dfb_this, percent_t_dfb = percent_t_dfb, dfb_tranges_actual = dfb_tranges_actual
					'efs_dsl': bbf_flux_cal, [time_start, time_end], dfb_tranges, probe = sc, e_folder = efs_folder, rt_pre = 2.5, rt_length = 1., flux_bbf = flux_bbf_this, flux_dfb = flux_dfb_this, secs_dfb = secs_dfb_this, percent_t_dfb = percent_t_dfb, dfb_tranges_actual = dfb_tranges_actual, vexb_dsl_folder = edtrd_folder
					endcase
	
					if finite(flux_bbf_this) and finite(flux_dfb_this) then begin
						if strcmp(abs_suf, '_abs') then begin
							flux_bbf_this = abs(flux_bbf_this)
							flux_dfb_this = abs(flux_dfb_this)
						endif
						secs_bbf = secs_bbf+secs_bbf_this
						secs_dfb = secs_dfb+secs_dfb_this
						flux_bbf = flux_bbf+flux_bbf_this
						flux_dfb = flux_dfb+flux_dfb_this
						n_events_used = n_events_used+1
					endif
				endif else begin
					;;;; need to get the actual dfb tranges anyway
					dfb_tranges_actual = ranges_compress(time_double(dfb_tranges), limit = [time_start, time_end])
				endelse
	
				;;;;;; get the values of different quantities
				if strmatch(quantities_check, '*V*') then begin
					no_use = value_tranges(dfb_tranges_actual, 'vi', probe = sc, datafolder=vi_folder, /less, combined_data = v_this)
					if strcmp(abs_suf, '_abs') then v_this = abs(v_this)
					v_dfb = [[v_dfb], [transpose(v_this)]]
				endif
				if strmatch(quantities_check, '*B*') then begin
					no_use = value_tranges(dfb_tranges_actual, 'fgs', probe = sc, datafolder=fgs_folder, /less, combined_data = b_this)
					if strcmp(abs_suf, '_abs') then b_this = abs(b_this)
					b_dfb = [[b_dfb], [transpose(b_this)]]
				endif
				if strmatch(quantities_check, '*Ey*') then begin
					no_use = value_tranges(dfb_tranges_actual, 'efs_dsl', probe = sc, rtrange = rtranges(*,i), datafolder=efs_folder, /reverse_bc, /less, combined_data = e_this, vexb_dsl_folder = edtrd_folder)
					if strcmp(abs_suf, '_abs') then e_this = abs(e_this)
					e_dfb = [[e_dfb], [transpose(e_this)]]
				endif

				;;;;; compare DFB and FB
				if strcmp(compare_fb, 'yes') then begin
					;;;;; get FBs
					case fb_select_type of
					'vi': fb_tranges = tr_flow_burst([time_start, time_end], probe = sc, v_folder = vi_folder, secs_fb = secs_fb_this, percent_t_fb = percent_t_fb_this)
					'vixb': fb_tranges = tr_flow_burst([time_start, time_end], probe = sc, v_folder = vi_folder, b_folder = fgs_folder, secs_fb = secs_fb_this, percent_t_fb = percent_t_fb_this)
					'efs_dsl': fb_tranges = tr_flow_burst([time_start, time_end], probe = sc, v_folder = vi_folder, e_folder = efs_folder, rt_pre = 2.5, rt_length = 1., secs_fb = secs_fb_this, percent_t_fb = percent_t_fb_this)
					endcase
					;;; seconds
					if finite(secs_fb_this) then secs_fb = secs_fb+secs_fb_this
					if finite(dfb_tranges_actual(0)) and finite(fb_tranges(0)) then begin
						trange_overlap = ranges_inter(dfb_tranges_actual, fb_tranges)
						if finite(trange_overlap(0)) then begin
							secs_overlap_this = total(trange_overlap(1,*)-trange_overlap(0,*))
							secs_overlap = secs_overlap+secs_overlap_this
						endif
					endif
				endif
		
				;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; plot the event to examine ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
				if strcmp(plot_events, 'yes') then begin
					trange_load = [time_start-minutes_load*60., time_end+minutes_load*60.] ;;; for DFB
					load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'fgs', datafolder = fgs_folder 
					load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'fgl', datafolder = fgl_folder 
					;; satellite positions
					load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'pos', datafolder = pos_folder
					;; pressure
					load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'Pall', datafolder = Pall_folder
					;; beta
					load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'beta', datafolder = beta_folder
					;; vi
					load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'vi', datafolder = vi_folder
					;; vixb
					load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'vixb', datafolder = vixb_folder
					;; ni
					load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'ni', datafolder = ni_folder
					;; velocity from efi
					;load_efi_data, trange = trange_load, probe = sc, rtrange = [time_start-180., time_start-30.], /tclip, e_folder = efs_folder, b_folder = fgs_folder
					;; efs, dsl
					load_efi_data, trange = trange_load, probe = sc, rtrange = [time_start-180., time_start-120.], /tclip, e_folder = efs_folder, /dsl, vexb_dsl_folder = edtrd_folder
					;; load plasma quantities
					;thm_load_esansst2, trange = , probe = sc
		
					;; smooth and get the B derivatives
					split_vec, 'th'+sc+'_fgs_gsm_tclip'
					split_vec, 'th'+sc+'_fgl_gsm_tclip'
					tsmooth2, 'th'+sc+'_fgs_gsm_tclip_z', sm_pts_fgs
					tsmooth2, 'th'+sc+'_fgl_gsm_tclip_z', sm_pts_fgl
					deriv_data, 'th'+sc+'_fgs_gsm_tclip_z'
					deriv_data, 'th'+sc+'_fgs_gsm_tclip_z_sm'
					deriv_data, 'th'+sc+'_fgl_gsm_tclip_z_sm'
					options, 'th'+sc+'_fgs_gsm_tclip', colors = [2,4,6], ytitle = 'B FGS', ysubtitle = '[nT]', labels = ['Vx', 'Vy', 'Vz']
					options, 'th'+sc+'_fgl_gsm_tclip', colors = [2,4,6], ytitle = 'B FGL', ysubtitle = '[nT]', labels = ['Vx', 'Vy', 'Vz']
					options, 'th'+sc+'_fgs_gsm_tclip_z_ddt', colors = 6, ytitle = 'dBz/dt FGS', ysubtitle = '[nT/s]'
					options, 'th'+sc+'_fgs_gsm_tclip_z_sm_ddt', colors = 6, ytitle = 'dBz/dt FGS SM', ysubtitle = '[nT/s]'
					options, 'th'+sc+'_fgl_gsm_tclip_z_sm_ddt', colors = 6, ytitle = 'dBz/dt FGL SM', ysubtitle = '[nT/s]'
					options, 'th'+sc+'_fgs_gsm_tclip_z_sm', colors = 6, ytitle = 'Bz FGS SM', ysubtitle = '[nT/s]'
					options, 'th'+sc+'_fgl_gsm_tclip_z_sm', colors = 6, ytitle = 'Bz FGL SM', ysubtitle = '[nT/s]'
					options, 'th'+sc+'_ptix_velocity_gsm_tclip', colors = [2,4,6], ytitle = 'Vi', ysubtitle = '[km/s]', labels = ['Vx', 'Vy', 'Vz']
		
					;; make another magnetic field
					copy_data, 'th'+sc+'_fgs_gsm_tclip', 'th'+sc+'_fgs_gsm_tclip1'
		
					;; plasma parameters
					options, 'th'+sc+'_Pall_tclip', colors=[2,4,0], labels=['P!db','P!dth','P!dall'], ylog=1
					options, 'th'+sc+'_beta_tclip', ylog=1
					options, 'th'+sc+'_ptix_density_tclip', ylog=1
					get_data, 'th'+sc+'_state_pos_tclip', data=data
					store_data, 'th'+sc+'_pos_gsm_re', data={x:data.x, y:data.y/6371.}
					split_vec, 'th'+sc+'_pos_gsm_re'
					options, 'th'+sc+'_pos_gsm_re_x', ytitle='XGSM'
					options, 'th'+sc+'_pos_gsm_re_y', ytitle='YGSM'
					options, 'th'+sc+'_pos_gsm_re_z', ytitle='ZGSM'
		
					;; transport
					split_vec, 'th'+sc+'_vixb_gsm_tclip'
					options, 'th'+sc+'_vixb_gsm_tclip_y', ytitle = '(VixB)_Y', ysubtitle = '[mV/m]'
					if tv_exist('th'+sc+'_efs_dsl_tclip') then begin
						get_data, 'th'+sc+'_efs_dsl_tclip', data = efs_dsl
						Ey = efs_dsl.y(*,1)
						if strcmp(sc, 'b') or strcmp(sc, 'c') then Ey=-Ey
						store_data, 'th'+sc+'_efs_dsl_tclip_y', data = {x:efs_dsl.x, y:Ey}
						options, 'th'+sc+'_efs_dsl_tclip_y', ytitle = 'Ey (DSL)', ysubtitle = '[mV/m]'
					endif
		
					;;;electrons: flux spectra, use burst
					if tv_exist('th'+sc+'_pseb_en_eflux') then begin
						store_data,'th'+sc+'_ptex_en_eflux',data='th'+sc+'_pseb_en_eflux th'+sc+'_peer_en_eflux'
						ylim,'th'+sc+'_ptex_en_eflux',5,1e6,1
						zlim,'th'+sc+'_ptex_en_eflux',1e2,5e8,1
					endif
		
					tplot, ['th'+sc+'_fgs_gsm_tclip', 'th'+sc+'_fgs_gsm_tclip1', 'th'+sc+'_fgs_gsm_tclip_z_ddt', 'th'+sc+'_fgs_gsm_tclip_z_sm_ddt', 'th'+sc+'_ptix_velocity_gsm_tclip', 'th'+sc+'_vixb_gsm_tclip_y', 'th'+sc+'_efs_dsl_tclip_y'], var_label = ['th'+sc+'_pos_gsm_re_z', 'th'+sc+'_pos_gsm_re_y', 'th'+sc+'_pos_gsm_re_x'], title = time_string(time_start, format = 6, precision = -1)+'_'+'th'+sc ;; to examine criteria
					timebar, time_start
					timebar, time_end
					timebar_mass, 0, varname=['th'+sc+'_fgs_gsm_tclip', 'th'+sc+'_fgs_gsm_tclip1', 'th'+sc+'_ptix_velocity_gsm_tclip', 'th'+sc+'_vixb_gsm_tclip_y', 'th'+sc+'_efs_dsl_tclip_y'], /databar
					timebar_mass, cri_ddt, varname='th'+sc+'_fgs_gsm_tclip_z*_ddt', /databar
		
					;;;; mark the DFBs
					if n_dfb_in ge 1 then begin
						for j = 0, n_dfb_in-1 do begin
							if strcmp(dfblist(2,j), 'm') then line = 2 else line = 1
							timebar_mass, dfblist(0,j), line = line, varname = ['th'+sc+'_fgs_gsm_tclip', 'th'+sc+'_fgs_gsm_tclip_z_sm_ddt']
						endfor ;; for of j
					endif
					if secs_dfb_this gt 0. then begin
						for j = 0, n_elements(dfb_tranges_actual(0,*))-1 do begin
							timebar_mass, dfb_tranges_actual(0,j), line = 0, varname = ['th'+sc+'_fgs_gsm_tclip1', 'th'+sc+'_vixb_gsm_tclip_y', 'th'+sc+'_efs_dsl_tclip_y']
							timebar_mass, dfb_tranges_actual(1,j), line = 1, varname = ['th'+sc+'_fgs_gsm_tclip1', 'th'+sc+'_vixb_gsm_tclip_y', 'th'+sc+'_efs_dsl_tclip_y']
						endfor ;; for of j
					endif
					;;;; write values
					xyouts, 0.15, 0.7, 'time ratio: '+strcompress(string(percent_t_dfb)), /normal
					xyouts, 0.15, 0.55, 'bflux ratio: '+strcompress(string(flux_dfb_this/flux_bbf_this)), /normal
					xyouts, 0.15, 0.45, 'BBF flux: '+strcompress(string(flux_bbf_this))+' mWb/m', /normal
		
		    		makepng, picfolder_events+'/bbf'+time_string(time_start, format = 6, precision = -1)+'_'+'th'+sc
				endif ;; if of plotting events
				;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			endif ;; if of having B data
		endfor ;; for of i, events
	
		if strmatch(quantities_check, '*time*') or strmatch(quantities_check, '*bflux*') then begin
			print, 'Events used:'+string(n_events_used)
			print, 'Ratio of time:'
			print, secs_dfb/secs_bbf
			print, 'Ratio of flux:'
			print, flux_dfb/flux_bbf
			ratio_secs(k) = secs_dfb/secs_bbf
			ratio_flux(k) = flux_dfb/flux_bbf
		endif
		;;;;; compare DFB and FB
		if strcmp(compare_fb, 'yes') then begin
			print, secs_fb/secs_bbf ;; the flow birst ratio, only one run is enough
			print, secs_overlap/secs_dfb
			print, secs_overlap/secs_fb
			ratio_inter_dfb(k) = secs_overlap/secs_dfb
			ratio_inter_fb(k) = secs_overlap/secs_fb
		endif
		;;;; get the average values
		if strmatch(quantities_check, '*V*') then begin
			v_dfb = v_dfb(*, 1:*)
			v_dfb_x = v_dfb(0,*)
			v_dfb_ttl = v_dfb(3,*)
			v_dfb_x_ave = mean(v_dfb_x, /nan)
			v_dfb_ttl_ave = mean(v_dfb_ttl, /nan)
			ratio_vx(k) = v_dfb_x_ave/v_bbf_x_ave
			ratio_vttl(k) = v_dfb_ttl_ave/v_bbf_ttl_ave
		endif
		if strmatch(quantities_check, '*B*') then begin
			b_dfb = b_dfb(*, 1:*)
			b_dfb_z = b_dfb(2,*)
			b_dfb_ttl = b_dfb(3,*)
			b_dfb_z_ave = mean(b_dfb_z, /nan)
			b_dfb_ttl_ave = mean(b_dfb_ttl, /nan)
			ratio_bz(k) = b_dfb_z_ave/b_bbf_z_ave
			ratio_bttl(k) = b_dfb_ttl_ave/b_bbf_ttl_ave
		endif
		if strmatch(quantities_check, '*E*') then begin
			e_dfb = e_dfb(*, 1:*)
			e_dfb_y = e_dfb(1,*)
			e_dfb_ttl = e_dfb(3,*)
			e_dfb_y_ave = mean(e_dfb_y, /nan)
			e_dfb_ttl_ave = mean(e_dfb_ttl, /nan)
			ratio_ey(k) = e_dfb_y_ave/e_bbf_y_ave
		endif
	endfor ;; for of k, the trys for ddt
	
	;;;;;;; plot the data
	case dfb_select_type of
	'bz': xrange = [0., 1.]
	'bz_nor2quiet_lobe': xrange = [0., 1.]/40.
	endcase
	case dfb_select_type of
	'bz': xunit = 'nT/s'
	'bz_nor2quiet_lobe': xunit = '/s'
	endcase
	
	;;;;;;; plot the dfb-fb ratio dependence on cri_ddt
	if strmatch(quantities_check, '*time*') and strmatch(quantities_check, '*bflux*') then begin
		;dataout_simple, save_folder+'/ratio_t_flux'+dtrd_suf+loc_suf+lim_suf, [[cri_ddts], [ratio_secs], [ratio_flux]]
		popen, pic_folder+'/'+filename_pre+'flux'+dtrd_suf+loc_suf(l)+dfb_select_suf+fb_select_suf+abs_suf+lim_suf
		print_options,xsize=xsize,ysize=ysize
		plot, cri_ddts, ratio_secs, xrange = xrange, yrange = [0., 1.], xtitle = 'DFB criterion of dB!dz!n/dt ['+xunit+']', ytitle = 'DFB/BBF Ratio', title = '', yticklen = 1, ygridstyle = 1;'DFB/FB Ratio '+title_str(l)+' DFB:'+dfb_select_type+' FB:'+fb_select_type
		oplot, cri_ddts, ratio_secs, psym = 4;, color = 2
		oplot, cri_ddts, ratio_flux, color = 6
		oplot, cri_ddts, ratio_flux, psym = 4, color = 6
		;;; add a 0.5 line for future use
		oplot, [0.5, 0.5], !y.crange, line = 1
		xyouts, 0.2, 0.3, 'duration', /data, align = 0.5
		xyouts, 0.6, 0.62, 'flux transport', /data, color = 6
		xyouts, 0.87, 0.87, '(e)', charsize = 1.2, /data
		if strcmp(compare_fb, 'yes') then begin
			oplot, cri_ddts, ratio_inter_dfb, color = 4
			oplot, cri_ddts, ratio_inter_dfb, psym = 4, color = 4
			oplot, cri_ddts, ratio_inter_fb, color = 2
			oplot, cri_ddts, ratio_inter_fb, psym = 4, color = 2
			xyouts, 0.6, 0.8, 'dfb-fb overlap/dfb', /data, color = 4
			xyouts, 0.6, 0.75, 'dfb-fb overlap/fb', /data, color = 2
			xyouts, 0.5, 0.7, 'fb-bbf ratio:'+strcompress(string(secs_fb/secs_bbf))
		endif
		pclose
	endif
	if strmatch(quantities_check, '*Ey*') then begin
		popen, pic_folder+'/'+filename_pre+'Ey'+loc_suf(l)+dfb_select_suf+fb_select_suf+abs_suf
		print_options,xsize=xsize,ysize=ysize
		plot, cri_ddts, ratio_ey, xrange = xrange, yrange = [0., 6.], xtitle = 'DFB criterion of dB!dz!n/dt ['+xunit+']', ytitle = '<DFB>/<BBF> Ratio', /nodata
		oplot, cri_ddts, ratio_ey, color = 4
		oplot, cri_ddts, ratio_ey, psym = 4, color = 4
		xyouts, 0.6, 0.8*!y.crange(1), 'Ey', /data, color = 4
		xyouts, 0.9, 0.9*!y.crange(1), '(b)'
		if strmatch(quantities_check, '*Bttl*') and strmatch(quantities_check, '*Vttl*') then begin
			oplot, cri_ddts, ratio_bttl, color = 6
			oplot, cri_ddts, ratio_bttl, psym = 4, color = 6
			xyouts, 0.6, 0.9*!y.crange(1), 'Bttl', /data, color = 6
			oplot, cri_ddts, ratio_vttl, color = 2
			oplot, cri_ddts, ratio_vttl, psym = 4, color = 2
			xyouts, 0.6, 0.85*!y.crange(1), 'Vttl', /data, color = 2
		endif
		if strmatch(quantities_check, '*Bz*') and strmatch(quantities_check, '*Vx*') then begin
			oplot, cri_ddts, ratio_bz, color = 6
			oplot, cri_ddts, ratio_bz, psym = 4, color = 6
			xyouts, 0.6, 0.9*!y.crange(1), 'Bz', /data, color = 6
			oplot, cri_ddts, ratio_vx, color = 2
			oplot, cri_ddts, ratio_vx, psym = 4, color = 2
			xyouts, 0.6, 0.85*!y.crange(1), 'Vx', /data, color = 2
		endif
		pclose
	endif
endfor ;; for of l, locations

stop
end
