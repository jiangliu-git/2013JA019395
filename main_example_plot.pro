pro main_example_plot
;;;;; present an example of bbf-embeded dfbs, along with the superposed Bz on determing DFB duration, used for Figure 1 of the paper.
;;; uses saved file from main_superpose.pro
thm_init
computer = 'I:'
;computer = '/home/jliu'
@declare

dtrd_suf = ''
;dtrd_suf = '_vexb'

if strcmp(dtrd_suf, '_vexb') then begin
	pref_e = ''
endif else begin 
	pref_e = delta_letter
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; draw the one-event example plot ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; choose the event
listname = 'bbf_ps_list_original2.txt'
events_all = load_list(listname, folder = listfolder)
;event = events_all(*, 100)
;; may be used:
;; see example folder in the publication pics of bbf_flux

;;; choose the event to use
;i_use = 524 ;; the "transport and go" event
i_use = 64 ;; the "transport and stay" event

;;; whether to use vexb for removeing electric field offset or 0 to remove
if strcmp(dtrd_suf, '_vexb') then edtrd_folder = vexb_dsl_folder else edtrd_folder = ''

;;; plot constants : no density and pressure
xsize = 6
ysize = 6.5
abc = ['(a)','(b)','(c)','(d)','(e)']
x_abc = 0.2
y_abc = [0.9, 0.73, 0.58, 0.4, 0.21]
;;; plot constants : with density and pressure
;xsize = 7
;ysize = 9
;abc = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)']
;x_abc = 0.2
;y_abc = [0.9, 0.73, 0.58, 0.4, 0.21, 0.1, 0.1]

;; some constants
minutes_load = 5
t_show = 0.4 ;; as a ratio of the bbf duration
sm_pts_fgs = 3 ;; used for DFCS paper
;sm_pts_fgs = 1 ;; used for more DFBs
sm_pts_fgl = 10
window_size = 0.5
pre_size = 0.5
cri_ddt = 0.5
minutes_select = 2.
c_secondary = 30.
flux_use = 'efs_dsl'
;; the duration of each DFB
dfb_length = 2./3. ;; in minutes, 40 secs
;dfb_length = 1. ;; in minutes

for i = i_use, i_use do begin
;for i = 0, n_elements(events_all(0,*))-1 do begin
	event = events_all(*,i)
	del_data, 'th*'
	time_start = time_double(event(0))
	time_end = time_double(event(1))
	sc = event(3)
	;;; get sc_num
	case sc of
	'a': sc_num = 5
	'b': sc_num = 1
	'c': sc_num = 2
	'd': sc_num = 3
	'e': sc_num = 4
	endcase
	trange_select = [time_start-minutes_select*60.,time_end+minutes_select*60.] ;;; for DFB
	load_bin_data, trange = trange_select, probe = sc, /tclip, datatype = 'fgs', datafolder = fgs_folder 
	if tv_exist('th'+sc+'_fgs_gsm_tclip') then begin
		secs_bbf_this = time_end-time_start ;; must be inside this if
		dfblist = dfb_select_bottom3('th'+sc+'_fgs_gsm_tclip', db = cri_ddt, window_size = window_size, sm_points = sm_pts, pre_size=pre_size, c_secondary = c_secondary)
		if strcmp(dfblist(0), 'no event') then begin
			dfb_tranges = [time_start-180., time_start-90]
			n_dfb_in = 0
		endif else begin
			dfb_tranges = [time_double(dfblist(0,*)), time_double(dfblist(0,*))+60.*dfb_length]
			n_dfb_in = n_elements(dfblist(0,*))
		endelse
		case flux_use of
		'vixb': bbf_flux_cal, [time_start, time_end], dfb_tranges, probe = sc, b_type = 'fgs', b_folder = fgs_folder, v_type = 'vi', v_folder = vi_folder, flux_bbf = flux_bbf_this, flux_dfb = flux_dfb_this, /x_only, secs_dfb = secs_dfb_this, percent_t_dfb = percent_t_dfb, dfb_tranges_actual = dfb_tranges_actual
		'efs_dsl': bbf_flux_cal, [time_start, time_end], dfb_tranges, probe = sc, e_folder = efs_folder, rt_pre = 2.5, rt_length = 1., flux_bbf = flux_bbf_this, flux_dfb = flux_dfb_this, secs_dfb = secs_dfb_this, percent_t_dfb = percent_t_dfb, dfb_tranges_actual = dfb_tranges_actual
		endcase

		flux_dfb_this = flux_dfb_this *6371000./1e9
		flux_bbf_this = flux_bbf_this *6371000./1e9
	
		trange_load = [time_start-minutes_load*60., time_end+minutes_load*60.] ;;; for DFB
		load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'pos', datafolder = pos_folder 
		load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'fgs', datafolder = fgs_folder 
		;; vi
		load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'vi', datafolder = vi_folder
		;; ni
		load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'ni', datafolder = ni_folder
		;; Pth
		load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'Pth', datafolder = Pth_folder
		split_vec, 'th'+sc+'_Pth_tclip'
		;; vixb
		load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'vixb', datafolder = vixb_folder
		;; efs, dsl
		load_efi_data, trange = trange_load, probe = sc, rtrange = [time_start-180., time_start-120.], /tclip, e_folder = efs_folder, /dsl, vexb_dsl_folder = edtrd_folder
	
		;; smooth and get the B derivatives
		split_vec, 'th'+sc+'_fgs_gsm_tclip'
		tsmooth2, 'th'+sc+'_fgs_gsm_tclip_z', sm_pts_fgs
		deriv_data, 'th'+sc+'_fgs_gsm_tclip_z'
		deriv_data, 'th'+sc+'_fgs_gsm_tclip_z_sm'
		options, 'th'+sc+'_fgs_gsm_tclip', colors = [2,4,6], ytitle = 'B!dGSM!n', ysubtitle = '!c[nT]', labels = ['B!dx!n', 'B!dy!n', 'B!dz!n'], labflag = 1
		options, 'th'+sc+'_fgs_gsm_tclip_z_sm_ddt', colors = 6, ytitle = 'dB!S!dz!R!Usm!n/dt', ysubtitle = '!c[nT/s]'
		options, 'th'+sc+'_ptix_velocity_gsm_tclip', colors = [2,4,6], ytitle = 'V!di!n', ysubtitle = '!c[km/s]', labels = ['V!dx!n', 'V!dy!n', 'V!dz!n'], labflag = 1

		;; density and pressure
		options, 'th'+sc+'_Pth_tclip_z', ytitle = 'P!dth!n', ysubtitle = '!c[nPa]'
		options, 'th'+sc+'_ptix_density_tclip', ytitle = 'n!di!n', ysubtitle = '!c[cm!U-3!n]'
	
		;; transport
		split_vec, 'th'+sc+'_vixb_gsm_tclip'
		options, 'th'+sc+'_vixb_gsm_tclip_y', ytitle = '(VixB)_Y', ysubtitle = '![mV/m]'
		if tv_exist('th'+sc+'_efs_dsl_tclip') then begin
			get_data, 'th'+sc+'_efs_dsl_tclip', data = efs_dsl
			t = efs_dsl.x
			Ey = efs_dsl.y(*,1)
			if strcmp(sc, 'b') or strcmp(sc, 'c') then Ey=-Ey
			dt = t(1:*)-t(0:n_elements(t)-2)
			flux_tran = [0,total(Ey(1:*)*dt, /cum)]*6371000./1e9 ;; convert to MWb/RE
			store_data, 'th'+sc+'_efs_dsl_tclip_y', data = {x:efs_dsl.x, y:Ey}
			store_data, 'th'+sc+'_intEy', data = {x:efs_dsl.x, y:flux_tran}
			options, 'th'+sc+'_efs_dsl_tclip_y', ytitle = pref_e+'E!dy,DSL'+"'"+'!n', ysubtitle = '!c[mV/m]'
			options, 'th'+sc+'_intEy', ytitle = phi_letter+'='+intg_sign+pref_e+'E!dy,DSL'+"'"+'!ndt', ysubtitle = '!c[MWb/R!dE!n]'
		endif

		options, '*', thick = l_thick
	
		bbf_dur = time_end-time_start
		;;;;;; only plot favorable events
		bbf_length_req = time_end-time_start gt 120.
		n_dfb_req = n_elements(dfblist(0,*)) ge 3
		combine_req = n_elements(dfb_tranges_actual(0,*)) lt n_elements(dfblist(0,*))
		combine_req2 = n_elements(dfb_tranges_actual(0,*)) gt 1
		flux_req1 = flux_dfb_this gt 0
		flux_req2 = flux_bbf_this gt 0
		flux_req3 = flux_dfb_this/flux_bbf_this lt 1.
		importance_req1 = flux_dfb_this/flux_bbf_this - percent_t_dfb gt 0.15
		importance_req2 = (flux_dfb_this/flux_bbf_this)/percent_t_dfb gt 1.5
		
		print, bbf_length_req
		print, n_dfb_req
		print, combine_req
		print, combine_req2
		print, flux_req1
		print, flux_req2
		print, flux_req3
		print, importance_req1
		print, importance_req2
	
		if bbf_length_req and n_dfb_req and combine_req and combine_req2 and flux_req1 and flux_req2 and flux_req3 and importance_req1 and importance_req2 then begin
			popen, pic_folder+'/example'+dtrd_suf
			print_options,xsize=xsize, ysize=ysize ;; use this for single plot
			tplot, ['th'+sc+'_ptix_velocity_gsm_tclip', 'th'+sc+'_fgs_gsm_tclip_z_sm_ddt', 'th'+sc+'_fgs_gsm_tclip', 'th'+sc+'_efs_dsl_tclip_y', 'th'+sc+'_intEy'], trange = [time_start-t_show*bbf_dur, time_end+t_show*bbf_dur] ;; to examine criteria
			;tplot, ['th'+sc+'_ptix_velocity_gsm_tclip', 'th'+sc+'_fgs_gsm_tclip_z_sm_ddt', 'th'+sc+'_fgs_gsm_tclip', 'th'+sc+'_ptix_density_tclip', 'th'+sc+'_Pth_tclip_z', 'th'+sc+'_efs_dsl_tclip_y', 'th'+sc+'_intEy'], title = 'An Example from P'+strcompress(string(sc_num),/remove), trange = [time_start-t_show*bbf_dur, time_end+t_show*bbf_dur] ;; to examine criteria
			timebar, time_start
			timebar, time_end
			timebar_mass, 0, varname=['th'+sc+'_fgs_gsm_tclip', 'th'+sc+'_ptix_velocity_gsm_tclip', 'th'+sc+'_efs_dsl_tclip_y'], /databar, line = 3
			timebar_mass, cri_ddt, varname='th'+sc+'_fgs_gsm_tclip_z_sm_ddt', /databar
	
			;;;; mark the DFBs
			if n_dfb_in ge 1 then begin
				for j = 0, n_dfb_in-1 do begin
					if strcmp(dfblist(2,j), 'm') then line = 1 else line = 1
					timebar, dfblist(0,j), line = line;, varname = ['th'+sc+'_fgs_gsm_tclip', 'th'+sc+'_fgs_gsm_tclip_z_sm_ddt']
				endfor ;; for of j
			endif
			;;;; mark the end of DFBs
			if finite(dfb_tranges_actual(0)) then begin
				for j = 0, n_elements(dfb_tranges_actual(0,*))-1 do begin
					;timebar_mass, dfb_tranges_actual(0,j), line = 1, varname = ['th'+sc+'_efs_dsl_tclip_y', 'th'+sc+'_intEy']
					;timebar_mass, dfb_tranges_actual(1,j), line = 3, varname = ['th'+sc+'_ptix_velocity_gsm_tclip', 'th'+sc+'_fgs_gsm_tclip', 'th'+sc+'_efs_dsl_tclip_y', 'th'+sc+'_intEy']
					timebar, dfb_tranges_actual(1,j), line = 3;, varname = ['th'+sc+'_ptix_velocity_gsm_tclip', 'th'+sc+'_fgs_gsm_tclip', 'th'+sc+'_efs_dsl_tclip_y', 'th'+sc+'_intEy']
					;; for density and pressure
					;timebar_mass, dfb_tranges_actual(1,j), varname=['th'+sc+'_ptix_density_tclip', 'th'+sc+'_Pth_tclip_z'], line = 3
				endfor ;; for of j
			endif
			;;;; write abcs
			xyouts, replicate(x_abc, 5), y_abc, abc, /normal
			;;;; write values
			;xyouts, 0.15, 0.17, 'Embedded DFBs lasts for'+strcompress(string(fix(percent_t_dfb*100.)))+'% of the BBF duration,!Cbut contribute'+strcompress(string(fix(flux_dfb_this/flux_bbf_this*100.)))+'% of the BBF flux transport.', /normal
			pclose
			;makepng, pic_folder+'/examples/example'+strcompress(string(i),/remove)
		endif ;; if of meeting good event criteria
	endif ;; if of existence of data
endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; draw the superposed all DFBs ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; the range for average
;aft_range = [2., 3.]
;;aft_range = [0.95, 2.]
;
;del_data, '*'
;
;sign_suf = '' ;; when doing both use this
;;sign_suf = '_negative'
;
;names_all = ['dbz_superpose_', 'vx_superpose_', 'ey_dsl_superpose_'+dtrd_suf]
;names_neg = ['vx_superpose_', 'ey_dsl_superpose_'+dtrd_suf]+'_negative'
;;names = [names_all, names_neg]
;names = names_all
;names = 'dbz_superpose_'
;
;for i = 0, n_elements(names)-1 do begin
;  name = names(i)
;  datain, name, filename = save_folder+'/'+name+'.dat', dim = 1
;endfor
;
;time_clip, 'dbz_superpose_'+sign_suf, time_double('2000 1 1')+aft_range(0)*60., time_double('2000 1 1')+aft_range(1)*60.
;get_data,  'dbz_superpose_'+sign_suf+'_tclip', data = after
;after_v = mean(after.y, /nan)
;
;get_data, 'dbz_superpose_'+sign_suf, data = bz
;t_end = bz.x(n_elements(bz.x)-1)
;no_use = max(bz.y, i_max, /nan)
;bz_behind = bz.y(i_max:*)
;t_behind = bz.x(i_max:*)
;i_neg = where(bz_behind-after_v lt 0)
;i_first_neg = i_neg(0)
;i_last_pos = i_neg(0)-1
;t_cut = interpol(t_behind([i_last_pos, i_first_neg]), bz_behind([i_last_pos, i_first_neg]), after_v)
;print, t_cut-time_double('2000 1 1')
;stop
;
;store_data, 'horiz_line', data = {x:[t_cut-0.06*(t_end-t_cut), t_end], y:[after_v, after_v]}
;options, 'horiz_line', thick = 0.7
;
;;; generate the plots for MAIN_SUPERPOSE, run it first
;options, 'dbz_superpose_*', ytitle=delta_letter+'B!dz!n [nT]'
;options, 'vi_superpose_*', ytitle='V [km/s]'
;options, 'vx_superpose_*', ytitle='V!dx!n [km/s]'
;options, 'ey_dsl_superpose_'+dtrd_suf+'*', ytitle=pref_e+'E!dy,DSL'+"'"+'!n [mV/m]'
;options, '*superpose*', xtickname=['-2','0','2'], thick=l_thick
;
;store_data, 'dbz_superpose_line', data=['dbz_superpose_'+sign_suf, 'horiz_line']
;
;names(0) = 'dbz_superpose_line'
;
;;;;; for Bz only ;;;;;;;;;;;;;;;;;;;
;popen, pic_folder+'/superpose_all'+dtrd_suf+sign_suf
;print_options,xsize=6,ysize=2.1
;tplot, names, trange = ['1999 12 31 23 57', '2000 1 1 0 3'], title = 'Medians of Superposed Profiles'
;timebar, '2000 1 1 0 0', line = 1
;timebar, t_cut, line = 1
;timebar_mass, 0, varname=names, /databar, line=3
;xyouts, 0.42,0.15, 'Minutes to t!d0!n', /normal
;pclose
;
;;;;; for all events only ;;;;;;;;;;;;;;;;;;;;;
;;popen, pic_folder+'/superpose_all'+dtrd_suf+sign_suf
;;print_options,xsize=6,ysize=4.2
;;tplot, names, trange = ['1999 12 31 23 57', '2000 1 1 0 3'], title = 'Medians of Superposed Profiles'
;;timebar, '2000 1 1 0 0', line = 1
;;timebar, t_cut, line = 1
;;timebar_mass, 0, varname=names, /databar, line=3
;;xyouts, 0.42,0.08, 'Minutes to t!d0!n', /normal
;;;;; label panels
;;xs=0.2+dblarr(3)
;;ys=[0.88, 0.6, 0.33]
;;ss=['(a)', '(b)', '(c)']
;;xyouts, xs, ys, ss, charsize=1.2, /normal
;;pclose
;
;;;;; for all events and negative events ;;;;;;;;;;;;;;;;;;;;;
;;popen, pic_folder+'/superpose_all'+dtrd_suf
;;print_options,xsize=6,ysize=6.2
;;tplot, names, trange = ['1999 12 31 23 57', '2000 1 1 0 3'], title = 'Medians of Superposed Profiles'
;;timebar, '2000 1 1 0 0', line = 1
;;timebar, t_cut, line = 1
;;timebar_mass, 0, varname=names, /databar, line=3
;;xyouts, 0.43,0.053, 'Minutes to t!d0!n', /normal
;;;;; label panels
;;xs=0.2+dblarr(5)
;;ys=[0.91, 0.73, 0.55, 0.33, 0.17]
;;ss=['(a)', '(b)', '(c)', '(d)', '(e)']
;;xyouts, xs, ys, ss, charsize=1.2, /normal
;;;;; make the labels
;;x_bar = 0.85
;;;; all
;;ver=0.45
;;hor=0.02
;;pos_right=[x_bar, 0.7]
;;make_px, pos_right, ver=ver, hor=hor, charname='All', lchar=0.08, orientation=-90, charsize = 1.3
;;;; negative
;;ver=0.29
;;hor=0.02
;;pos_right=[x_bar, 0.275]
;;make_px, pos_right, ver=ver, hor=hor, charname='Negative Only', lchar=0.2, orientation=-90, charsize = 1.08
;;pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
stop
end
