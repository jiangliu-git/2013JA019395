pro main_ratioacvsx
;;; see how much flux transport is contributed by DFBs in BBFs
;;; select DFBs on the fly

thm_init

computer = 'I:'
;computer = '/home/jliu'

@declare

listname = 'bbf_ps_list_original2'

;;; for x dependence one should only use 0809 events.
;season_suf = ''
season_suf = '_0809'

;;; whether to use vexb for removeing electric field offset or 0 to remove
dtrd_suf = ''
;dtrd_suf = '_vexb'
if strcmp(dtrd_suf, '_vexb') then edtrd_folder = vexb_dsl_folder else edtrd_folder = ''

;;;;; get the absolute value or not
abs_suf = ''
;abs_suf = '_abs'


;;;;; specify different choices
;;; choose how to select DFBs (normalize or not)
dfb_select_type = 'bz'
;dfb_select_type = 'bz_nor2blobeq'
;dfb_select_type = 'bz_trend'
;dfb_select_type = 'bz_nor2bqz'

;;; choose which DFB duration to use
dfb_dur_type = 'fixed'
;dfb_dur_type = 'xdependent'

;;; choose which flux to use for magnetic flux
flux_use = 'efs_dsl'
;flux_use = 'vixb'

;;;;; specify run of locations
;locs = [-30., -6.] ;;; for test
;locs = [-30., -15., -9., -6.] ;;; for test
;; different settings
if strcmp(dfb_select_type, 'bz') and strcmp(dfb_dur_type, 'fixed') then begin
	locs = [-30., -23, -18.2, -16.,-13.4,-11,-9,-6]
endif
if strcmp(dfb_select_type, 'bz') and strcmp(dfb_dur_type, 'xdependent') then begin
	locs = [-30., -24, -20, -16,-13,-11,-8.7,-6]
endif
if strcmp(dfb_select_type, 'bz_nor2blobeq') and strcmp(dfb_dur_type, 'fixed') then begin
	locs = [-30., -23, -21, -18,-16,-12.5,-10.3,-8,-6]
endif
if strcmp(dfb_select_type, 'bz_nor2blobeq') and strcmp(dfb_dur_type, 'xdependent') then begin
	locs = [-30., -22, -17,-13,-11,-9,-6] ;; done, use for paper is used
endif
if strcmp(dfb_select_type, 'bz_trend') and strcmp(dfb_dur_type, 'xdependent') then begin
	locs = [-30., -24, -20., -16,-13,-11,-9,-6]
;	locs = [-30., -24, -19.9, -16,-13,-11,-9,-6]
endif
if strcmp(dfb_select_type, 'bz_nor2bqz') and strcmp(dfb_dur_type, 'xdependent') then begin
	locs = [-30., -22, -17,-13,-11,-9,-6]
endif

;;; set the minimum bin number
k_c = 5

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
;; the fixed duration of each DFB
dfb_dur_fix = 2./3. ;; in minutes

;; file name
filename_pre = 'dfb_bbf_acum_ratio_'

;;; choose fgs or fgl
data_dfb = 'fgs'
bfolder = fgs_folder
sm_pts = sm_pts_fgs

case dfb_select_type of
'bz': dfb_select_suf = ''
'bz_nor2blobeq': dfb_select_suf = '_dfbnor2blobeq'
'bz_trend': dfb_select_suf = '_dfbnor2bztrend'
'bz_nor2bqz': dfb_select_suf = '_dfbnor2bqz'
endcase
case dfb_select_type of
'bz': cri_ddt = 0.5
'bz_nor2blobeq': cri_ddt = 0.5/35.520242 ;; the determinant is the median Blobe for ALL L13a events with plasma data
'bz_trend': cri_ddt = 0.5/8.9961452 ;; the determinant is the median peak dBz for ALL L13a events with fgs
'bz_nor2bqz': cri_ddt = 0.5/8.8975143 ;; the determinant is the median Bqz for ALL L13a events with fgs
endcase

case dfb_dur_type of
'fixed': dur_suf = ''
'xdependent': dur_suf = '_xdur'
endcase

events_all = load_list(listname+'.txt', folder = listfolder)
;;;; for test
;events_all = events_all(*,100:130)

if strcmp(season_suf, '_0809') then begin
	i_0809 = where(time_double(events_all(0,*)) lt time_double('2010 1 1'))
	events_all = events_all(*, i_0809)
endif

ratio_secs = dblarr(n_elements(locs)-1)
ratio_flux = dblarr(n_elements(locs)-1)
n_xbin = dblarr(n_elements(locs)-1)

;;;;;;;; get the event location and seperate
pos_list = event_status_instant(events_all, 'pos', vtrange = time_string(events_all(0:1,*)), datafolder=pos_folder)
x = pos_list(2,*)/6371.

for l = 0, n_elements(locs)-2 do begin
	;;; get the events in this location
	i_use = where((x gt locs(l)) and (x le locs(l+1)), j_use)
	if j_use gt 0 then begin
		events = events_all(*, i_use)
		;;; get the quiet time lobe field for DFB selection
		if strcmp(dfb_select_type, 'bz_nor2blobeq') then h = scale_height(events, b_folder = fgs_folder, p_folder = Pall_folder, time_length = 1., pre_time = 2.5, Blobe = Bnor)
		if strcmp(dfb_select_type, 'bz_trend') then Bnor = bz_trend(events, method = 'x', datafolder = pos_folder)
		if strcmp(dfb_select_type, 'bz_nor2bqz') then begin
			bq_list = event_status_instant(events, 'fgs', pre_time=2.5, time_length=1., datafolder=fgs_folder)
			Bnor = bq_list(4,*)
		endif
		;;;;;;;;;; begin surveying
		secs_bbf = 0.
		secs_dfb = 0.
		flux_bbf = 0.
		flux_dfb = 0.
		n_events_used = 0
		;;;; begin survey the events
		for i = 0, n_elements(events(0,*))-1 do begin
			del_data, 'th*'
			event = events(*,i)
			case dfb_select_type of
			'bz': nor_v = 1.
			else: nor_v = Bnor(i)
			endcase
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
					if strcmp(dfb_dur_type, 'fixed') then begin
						dfb_dur = dfb_dur_fix
						dfb_tranges = [time_double(dfblist(0,*)), time_double(dfblist(0,*))+60.*dfb_dur]
					endif else begin
						dfblist_4c = [dfblist, replicate(sc, 1, n_elements(dfblist(0,*)))]
						if strcmp(dfb_dur_type, 'xdependent') then durs = dfb_duration(dfblist_4c, method = 'x', datafolder = pos_folder)
						dfb_tranges = [time_double(dfblist(0,*)), time_double(dfblist(0,*))+durs]
					endelse
					n_dfb_in = n_elements(dfblist(0,*))
				endelse
		
				;;;;; calculate the flux of DFB and BBF
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
			endif ;; if of having B data
		endfor ;; for of i, events
		
		n_xbin(l) = n_events_used
		print, 'Events used:'+string(n_events_used)
		if n_events_used gt 0 then begin
			print, 'Ratio of time:'
			print, secs_dfb/secs_bbf
			print, 'Ratio of flux:'
			print, flux_dfb/flux_bbf
			ratio_secs(l) = secs_dfb/secs_bbf
			ratio_flux(l) = flux_dfb/flux_bbf
		endif else begin
			ratio_secs(l) = !values.f_nan
			ratio_flux(l) = !values.f_nan
		endelse
	endif else begin ;; if of whether there is bbf in this location.
		ratio_secs(l) = !values.f_nan
		ratio_flux(l) = !values.f_nan
		n_xbin(l) = 0.
	endelse
endfor ;; for of l, locations

;;;;;;; plot the dfb-bbf ratio dependence on cri_ddt
loc_values = 0.5*(locs(1:*)+locs(0:n_elements(locs)-2))
print, loc_values
print, n_xbin
dataout_simple, save_folder+'/ratio_t_flux_vsx'+dtrd_suf+'_'+dfb_select_type+dur_suf+season_suf, transpose([[loc_values], [ratio_secs], [ratio_flux], [n_xbin]])
i_bad = where(n_xbin le k_c, j_bad)
if j_bad gt 0 then begin
	ratio_secs(i_bad) = !values.f_nan
	ratio_flux(i_bad) = !values.f_nan
endif
plot, loc_values, ratio_secs, xrange = [-6, -30], yrange = [0., 1.], xtitle = 'X [RE]', ytitle = 'DFB/BBF Ratio', title = '', yticklen = 1, ygridstyle = 1;'DFB/FB Ratio '+title_str(l)+' DFB:'+dfb_select_type+' FB:'+fb_select_type
oplot, loc_values, ratio_secs, psym = 4;, color = 2
oplot, loc_values, ratio_flux, color = 6
oplot, loc_values, ratio_flux, psym = 4, color = 6
;;; add a 0.5 line for future use
oplot, [0.5, 0.5], !y.crange, line = 1
xyouts, 0.2, 0.3, 'duration', /data, align = 0.5
xyouts, 0.6, 0.62, 'flux transport', /data, color = 6
makepng, pic_folder+'/ratio_v_x'+dtrd_suf+'_'+dfb_select_type+dur_suf+season_suf

stop
end
