pro main_bbf_events
;; look at events
thm_init

;; the event
computer = 'I:'
;computer = '/home/jliu'

;cd, computer+'/My Dropbox/Work/DFB/statistics'

;; the folder to put pics
dataroot_folder = '/Work/data'
plasma_folder = '/plasma_new'

fgs_folder=computer + '/Work/data/fgs'
fgh_folder=computer + '/Work/data/fgh'
fgl_folder=computer + '/Work/data/fgl'
pos_folder=computer + '/Work/data/pos'
efs_folder=computer + '/Work/data/efs_dsl'
Pallfolder=computer + dataroot_folder + plasma_folder + '/Pall'
beta_folder=computer + dataroot_folder + plasma_folder + '/beta'
vperp_folder = computer + dataroot_folder + plasma_folder + '/vperp'
vi_folder = computer + dataroot_folder + plasma_folder + '/vi'
vixb_folder = computer + dataroot_folder + plasma_folder + '/vixb'
ni_folder = computer + dataroot_folder + plasma_folder + '/ni'

listfolder='../../lists'
picfolder = '../../pics_test/events'

listname = 'bbf_ps_list_original2.txt'
;listname = 'dfb_list_lobe_test.txt'
;listname = 'dfb_list_lead_tail_tailward_df_binxbout.txt'
;listname = 'dfb_list_lead_tail_sharp.txt'
;listname = 'dfb_list_standalone_tail.txt'
dfb_list = 'dfb_ignore2_list_original.txt'

events_all = load_list(listname, folder = listfolder)
events_dfb = load_list(dfb_list, folder = listfolder)
;; use max time
;t_dfb = time_double(events_dfb(1,*))
;; use onset time
t_dfb = time_double(events_dfb(0,*))

;; the datatype and event
minutes = 5
sm_pts_fgs = 3
sm_pts_fgl = 10
cri_ddt = 0.5

;; 
events = events_all(*,*)
;;; b_in as the peak value right behind the front
;b_in_list = event_status_instant(events, datatype, time_length =0.25, pre_time = -0.125, datafolder=datafolder, /max, comp_v = 2, time_interest = t_max_arr)
;;; b_out as the min value right ahead of the front
;b_out_list = event_status_instant(events, datatype, vtrange = [time_double(events(0,*))-15,t_max_arr(1,*)], datafolder=datafolder, /min, comp_v = 2, time_interest = t_min_arr)

for i = 0, n_elements(events(0,*))-1 do begin
	del_data, '*'
	event = events(*,i)
	time = time_double(event(0))
	;trange_load = [time-minutes*60.,time+minutes*60.] ;;; for DFB
	trange_load = [time-minutes*60., time_double(event(1))+minutes*60.] ;;; for DFB
	sc = event(3)
	load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'fgs', datafolder = fgs_folder 
	load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'fgl', datafolder = fgl_folder 
	;; satellite positions
	load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'pos', datafolder = pos_folder
	;; pressure
	load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'Pall', datafolder = Pallfolder
	;; beta
	load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'beta', datafolder = beta_folder
	;; vi
	load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'vi', datafolder = vi_folder
	;; vixb
	load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'vixb', datafolder = vixb_folder
	;; ni
	load_bin_data, trange = trange_load, probe = sc, /tclip, datatype = 'ni', datafolder = ni_folder
	;; velocity from efi
	load_efi_data, trange = trange_load, probe = sc, rtrange = [time-180., time-30.], /tclip, e_folder = efs_folder, b_folder = fgs_folder
	;; load plasma quantities
	;thm_load_esansst2, trange = , probe = sc
	catch, err
    if err eq 0 then begin
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

		;;;electrons: flux spectra, use burst
		if tv_exist('th'+sc+'_pseb_en_eflux') then begin
			store_data,'th'+sc+'_ptex_en_eflux',data='th'+sc+'_pseb_en_eflux th'+sc+'_peer_en_eflux'
			ylim,'th'+sc+'_ptex_en_eflux',5,1e6,1
			zlim,'th'+sc+'_ptex_en_eflux',1e2,5e8,1
		endif

		;; make the limits for velocities
		
		;tplot, ['th'+sc+'_'+datatype+'_gsm_all', 'th'+sc+'_ptix_vperp_gsm', 'th'+sc+'_ptex_vperp_gsm', 'th'+sc+'_efs_velocity_gsm', 'th'+sc+'_ptix_en_eflux', 'th'+sc+'_ptex_en_eflux', 'th'+sc+'_psif_an_eflux_phi', 'th'+sc+'_peir_an_eflux_phi'], var_label = ['th'+sc+'_pos_gsm_re_z', 'th'+sc+'_pos_gsm_re_y', 'th'+sc+'_pos_gsm_re_x']
		;tplot, ['th'+sc+'_'+datatype+'_gsm_all', 'th'+sc+'_ptix_velocity_gsm_tclip', 'th'+sc+'_ptix_density_tclip', 'th'+sc+'_beta_tclip', 'th'+sc+'_Pall_tclip', 'th'+sc+'_vixb_gsm_tclip_y'], var_label = ['th'+sc+'_pos_gsm_re_z', 'th'+sc+'_pos_gsm_re_y', 'th'+sc+'_pos_gsm_re_x']
		tplot, ['th'+sc+'_fgs_gsm_tclip', 'th'+sc+'_fgl_gsm_tclip', 'th'+sc+'_fgs_gsm_tclip_z_sm', 'th'+sc+'_fgl_gsm_tclip_z_sm', 'th'+sc+'_fgs_gsm_tclip_z_ddt', 'th'+sc+'_fgs_gsm_tclip_z_sm_ddt', 'th'+sc+'_fgl_gsm_tclip_z_sm_ddt', 'th'+sc+'_ptix_velocity_gsm_tclip', 'th'+sc+'_vixb_gsm_tclip_y'], var_label = ['th'+sc+'_pos_gsm_re_z', 'th'+sc+'_pos_gsm_re_y', 'th'+sc+'_pos_gsm_re_x']
		timebar, time
		timebar, time_double(event(1))
		;timebar, 0, varname='th'+sc+'_'+datatype+'_gsm_all', /databar
		timebar, 0, varname='th'+sc+'_ptix_velocity_gsm_tclip', /databar
		timebar, 0, varname='th'+sc+'_vixb_gsm_tclip_y', /databar
		timebar, cri_ddt, varname='th'+sc+'_fgs_gsm_tclip_z_ddt', /databar
		timebar, cri_ddt, varname='th'+sc+'_fgs_gsm_tclip_z_sm_ddt', /databar
		timebar, cri_ddt, varname='th'+sc+'_fgl_gsm_tclip_z_sm_ddt', /databar
		;timebar, 0, varname='th'+sc+'_ptix_vperp_gsm', /databar
		;timebar, 0, varname='th'+sc+'_ptex_vperp_gsm', /databar
		;timebar, 0, varname='th'+sc+'_efs_velocity_gsm', /databar
		;;; examine max and min value
		;timebar, t_max_arr(1,i), line=1
		;timebar, t_min_arr(1,i), line=1
    	;makepng, picfolder+'/dfb'+time_string(time, format = 6, precision = -1)+'_'+'th'+sc

		;;;; find the DFBs inside the bbf
		;pre_min = 1. ;; include also a pre time
		;i_dfb = where((t_dfb gt time - pre_min*60.) and (t_dfb lt time_double(event(1))), n_dfb_in)
		;if n_dfb_in ge 1 then begin
		;	for j = 0, n_dfb_in-1 do begin
		;		i_dfb_this = i_dfb(j)
		;		if strcmp(events_dfb(2,i_dfb_this), 'm') then line = 2 else line = 1
		;		timebar, t_dfb(i_dfb_this), line = line;, varname = 'th'+sc+'_'+datatype+'_gsm_all'
		;		;timebar, t_dfb(i_dfb_this), line = line, varname = 'th'+sc+'_vixb_gsm_tclip_y'
		;	endfor
		;endif
		;;;;;;;;;;;;

    	makepng, picfolder+'/bbf'+time_string(time, format = 6, precision = -1)+'_'+'th'+sc
    endif else begin
      dprint, !error_state.msg
    endelse
    catch, /cancel
endfor
end
