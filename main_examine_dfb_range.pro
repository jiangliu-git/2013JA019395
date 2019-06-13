pro main_examine_dfb_range
;;;; use this interactive routine to examine DFB range
thm_init

;; the event
;computer = 'I:'
computer = '/home/jliu'

@declare

listname = 'dfb_list_lead_tail'

events = load_list(listname+'.txt', folder = listfolder)
;;;; only use part of them
events = events(*,102:150)

;;;; prior and after the dfb to show in minutes
minutes_load_pre = 2.
minutes_load_aft = 3.

;;;; smooth points
sm_pts_fgs = 3 ;; used for DFCS paper
sm_pts_fgl = 10

for i = 0, n_elements(events(0,*))-1 do begin 
	time_start = time_double(events(0,i))
	time_end = time_double(events(1,i))
	sc = events(3,i)
	trange_load = [time_start-minutes_load_pre*60., time_end+minutes_load_aft*60.] ;;; for DFB
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
	load_efi_data, trange = trange_load, probe = sc, rtrange = [time_start-180., time_start-120.], /tclip, e_folder = efs_folder, /dsl, /reverse
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
	split_vec, 'th'+sc+'_efs_dsl_tclip'
	options, 'th'+sc+'_efs_dsl_tclip_y', ytitle = 'Ey (DSL)', ysubtitle = '[mV/m]'
	options, 'th'+sc+'_intEy_dsl_tclip', ytitle = 'Flux: intEy (DSL)', ysubtitle = '[mWb/m]'
	
	;;;electrons: flux spectra, use burst
	if tv_exist('th'+sc+'_pseb_en_eflux') then begin
		store_data,'th'+sc+'_ptex_en_eflux',data='th'+sc+'_pseb_en_eflux th'+sc+'_peer_en_eflux'
		ylim,'th'+sc+'_ptex_en_eflux',5,1e6,1
		zlim,'th'+sc+'_ptex_en_eflux',1e2,5e8,1
	endif
	
	tplot, ['th'+sc+'_fgs_gsm_tclip', 'th'+sc+'_fgs_gsm_tclip_z_sm_ddt', 'th'+sc+'_ptix_velocity_gsm_tclip', 'th'+sc+'_vixb_gsm_tclip_y', 'th'+sc+'_efs_dsl_tclip_y', 'th'+sc+'_intEy_dsl_tclip'], var_label = ['th'+sc+'_pos_gsm_re_z', 'th'+sc+'_pos_gsm_re_y', 'th'+sc+'_pos_gsm_re_x'], title = time_string(time_start, format = 6, precision = -1)+'_'+'th'+sc ;; to examine criteria
	timebar, time_start
	timebar_mass, 0, varname=['th'+sc+'_fgs_gsm_tclip', 'th'+sc+'_ptix_velocity_gsm_tclip', 'th'+sc+'_vixb_gsm_tclip_y', 'th'+sc+'_efs_dsl_tclip_y'], /databar
	ctime, time_end
	timebar, time_end, line = 1
	
	makepng, pic_folder_events+'/bbf'+time_string(time_start, format = 6, precision = -1)+'_'+'th'+sc
endfor
end
