pro main_relations_dfb
;;; examine all kinds of relations (flow, field strength, etc)

thm_init

;;; the event
computer = 'I:'
;computer = '/home/jliu'

@declare

listname_original = 'dfb_list_lead_tail'

;;; select list
;suffix = ''
suffix = '_fgs'
;suffix = '_earthward_df_binxbout'

listname = listname_original+suffix+'.txt'

events_all = load_list(listname, folder = listfolder)
events = events_all
;events = events_all(*,386) ;; this event has a gap at t0-3min, t0-2min of FGS

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; run
;;;;;;;;;;;; get the event location and seperate
;pos_list = event_status_instant(events, 'pos', time_length=3, datafolder=pos_folder)
;x = pos_list(2,*)/6371.
;;y = pos_list(3,*)/6371.
;;z = pos_list(4,*)/6371.
;;;dataout_simple, save_folder+'/pos'+suffix, [x,y,z]
;;;dataout_simple, save_folder+'/dfb_x'+suffix, x
;;dataout_simple, save_folder+'/dfb_y'+suffix, y
;;stop

;;;;;;;; get event scale height
h = scale_height(events, b_folder = fgs_folder, p_folder = Pall_folder, Blobe = Bl)
print, median(bl, /even)
stop
;;dataout_simple, save_folder+'/dfb_h'+suffix, h
;;;;; load
;h = datain_simple(save_folder+'/dfb_h'+suffix+'.dat', dim=1, type='double')

;;;;;;;; get normal direction
;;n = df_normal(events, method=3, datatype='fgl', bfolder=fgl_folder, c_angle = 30)
;;dataout_simple, save_folder+'/n_binxbout'+suffix, n
;;;; load
;n = datain_simple(save_folder+'/n_binxbout'+suffix+'.dat', dim=3, type='double')
;;;; components
;ny = n(1,*)

;;;;;;;;;;;;;;;  ground signatures
;;;;;;;; get ae index
;kyoto_ae_peak_list = event_status_instant(events, 'kyoto_ae', pre_time = 0, time_length=3., datafolder=kyoto_ae_folder, /peak)
;kyoto_ae_peak = kyoto_ae_peak_list(2,*)

;;;;;;;; get ae index
;kyoto_ae_ave_list = event_status_instant(events, 'kyoto_ae', pre_time = 0, time_length=3., datafolder=kyoto_ae_folder)
;kyoto_ae_ave = kyoto_ae_ave_list(2,*)

;;;;;;;; get ae index
;pseudo_ae_peak_list = event_status_instant(events, 'pseudo_ae', pre_time = -3, time_length=3., datafolder=pseudo_ae_folder, /peak)
;pseudo_ae_peak = pseudo_ae_peak_list(2,*)

;;;;;;;;; get ae index
;pseudo_ae_ave_list = event_status_instant(events, 'pseudo_ae', pre_time = -3., time_length=3., datafolder=pseudo_ae_folder)
;pseudo_ae_ave = pseudo_ae_ave_list(2,*)
;;dataout_simple, save_folder+'/pseudo_ae_ave'+suffix, pseudo_ae_ave
;;stop

;;;;;;;; get al index
;pseudo_al_ave_list = event_status_instant(events, 'pseudo_al', pre_time = -3., time_length=3., datafolder=pseudo_al_folder)
;pseudo_al_ave = pseudo_al_ave_list(2,*)
;dataout_simple, save_folder+'/pseudo_al_ave'+suffix, pseudo_al_ave
;stop

;;;;;;;;;;;; solar wind parameters
;;;;;;;; get IMF
;;;; peak IMF
;omni_b_peak_list = event_status_instant(events, 'omni_b_gsm', pre_time = 10, time_length=15., datafolder=omni_b_folder, /peak)
;omni_bz_peak = omni_b_peak_list(4,*)
;
;;;; mean IMF
;omni_b_ave_list = event_status_instant(events, 'omni_b_gsm', pre_time = 10, time_length=15, datafolder=omni_b_folder)
;omni_bz_ave = omni_b_ave_list(4,*)
;
;;;;;;;;; get Pdyn
;;;; peak Pdyn
;omni_Pdyn_peak_list = event_status_instant(events, 'omni_Pdyn', pre_time = 10, time_length=15, datafolder=omni_Pdyn_folder, /peak)
;omni_Pdyn_peak = omni_Pdyn_peak_list(2,*)
;;
;;;;; mean Pdyn
;omni_Pdyn_ave_list = event_status_instant(events, 'omni_Pdyn', pre_time = 10, time_length=15, datafolder=omni_Pdyn_folder)
;omni_Pdyn_ave = omni_Pdyn_ave_list(2,*)

;;;;;;;; get electric field
;;; peak Pdyn
;omni_vxb_peak_list = event_status_instant(events, 'omni_vxb', pre_time = 10, time_length=5, datafolder=omni_vxb_folder, /peak)
;omni_vxb_y_peak = omni_vxb_peak_list(3,*)
;;
;;;;; mean vxb
;omni_vxb_ave_list = event_status_instant(events, 'omni_vxb', pre_time = 6, time_length=6, datafolder=omni_vxb_folder)
;omni_vxb_y_ave = omni_vxb_ave_list(3,*)

;;;;;;;;;;;; velocity
;;;; peak velocity
;; 40s
;v_peak_list = event_status_instant(events, 'vi', pre_time = -1./3., time_length=2./3., datafolder=vi_folder, /peak)
;dur_suf = '_40s'
;;; 2min
;;v_peak_list = event_status_instant(events, 'vi', pre_time = -1., time_length=2., datafolder=vi_folder, /peak)
;v_peak_x = v_peak_list(2,*)
;v_peak_y = v_peak_list(3,*)
;v_peak_z = v_peak_list(4,*)
;v_peak = v_peak_list(5,*)
;dataout_simple, save_folder+'/dfb_vx_peak'+dur_suf+suffix, v_peak_x
;stop

;;; average velocity
;v_ave_list = event_status_instant(events, 'vi', time_length=4, datafolder=vi_folder)
;v_ave_x = v_ave_list(2,*)
;v_ave_y = v_ave_list(3,*)
;v_ave_z = v_ave_list(4,*)
;v_ave = v_ave_list(5,*)

;;;; perp velocity
;;;;; peak velocity
;vperp_peak_list = event_status_instant(events, 'viperp', pre_time = -0.25, time_length=0.5, datafolder=viperp_folder, /peak)
;vperp_peak_x = vperp_peak_list(2,*)
;vperp_peak_y = vperp_peak_list(3,*)
;vperp_peak_z = vperp_peak_list(4,*)
;vperp_peak_ttl = vperp_peak_list(5,*)
;i_fast = where(finite(vperp_peak_x))
;output_txt, events(*,i_fast), filename = 'dfb_list_lead_tail_fs.txt'
;stop

;;; average perp velocity
;vperp_ave_list = event_status_instant(events, 'viperp', time_length=4, datafolder=viperp_folder)
;vperp_ave_x = vperp_ave_list(2,*)
;vperp_ave_y = vperp_ave_list(3,*)
;vperp_ave_z = vperp_ave_list(4,*)
;vperp_ave_ttl = vperp_ave_list(5,*)

;;;; xy velocity
;;;; peak velocity
;vxy_peak_list = event_status_instant(events, 'vixy_infer', pre_time = -0.25, time_length=0.5, datafolder=vixy_folder, /peak)
;vxy_peak_x = vxy_peak_list(2,*)
;vxy_peak_y = vxy_peak_list(3,*)
;vxy_peak_z = vxy_peak_list(4,*)
;vxy_peak = vxy_peak_list(5,*)

;; average xy velocity
;vxy_ave_list = event_status_instant(events, 'vixy_infer', pre_time = -0.25, time_length=0.5, datafolder=vixy_folder)
;vxy_ave_x = vxy_ave_list(2,*)
;vxy_ave_y = vxy_ave_list(3,*)
;vxy_ave_z = vxy_ave_list(4,*)
;vxy_ave = vxy_ave_list(5,*)

;;; inferred vperp from efi
;;;;;;;;;;;; peak
;;;vperp_efi_peak_list = event_status_instant_efi(events, 'vperp', pre_time = -0.25, time_length=0.5, rpre_time = 2.5, rtime_length=1., /peak, e_folder=efs_folder, b_folder = fgs_folder)
;;;vperp_efi_peak = vperp_efi_peak_list(2:5,*)
;;;dataout_simple, save_folder+'/vperp_efi_peak'+suffix, vperp_efi_peak
;;;;; load
;vperp_efi_peak = datain_simple(save_folder+'/vperp_efi_peak'+suffix+'.dat', dim=4, type='double')
;vperp_efi_peak_x = vperp_efi_peak(0,*)
;vperp_efi_peak_ttl = vperp_efi_peak(3,*)

;;;;;;;;;;; ave
;;vperp_efi_ave_list = event_status_instant_efi(events, 'vperp', pre_time = -0.25, time_length=0.5, rpre_time = 2.5, rtime_length=1., e_folder=efs_folder, b_folder = fgs_folder)
;;vperp_efi_ave = vperp_efi_ave_list(2:5,*)
;;dataout_simple, save_folder+'/vperp_efi_ave'+suffix, vperp_efi_ave
;;; load
;vperp_efi_ave = datain_simple(save_folder+'/vperp_efi_ave'+suffix+'.dat', dim=4, type='double')
;vperp_efi_ave_x = vperp_efi_ave(0,*)
;vperp_efi_ave_ttl = vperp_efi_ave(3,*)

;;; inferred vxy from efi
;;;;;;;;;;;; peak
;;vxy_efi_peak_list = event_status_instant_efi(events, 'vxy', pre_time = -0.25, time_length=0.5, rpre_time = 2.5, rtime_length=1., /peak, e_folder=efs_folder, b_folder = fgs_folder)
;;vxy_efi_peak = vxy_efi_peak_list(2:5,*)
;;dataout_simple, save_folder+'/vxy_efi_peak'+suffix, vxy_efi_peak
;;;; load
;vxy_efi_peak = datain_simple(save_folder+'/vxy_efi_peak'+suffix+'.dat', dim=4, type='double')
;vxy_efi_peak_x = vxy_efi_peak(0,*)

;;;;;;;;;;; ave
;vxy_efi_ave_list = event_status_instant_efi(events, 'vxy', pre_time = -0.25, time_length=0.5, rpre_time = 2.5, rtime_length=1., e_folder=efs_folder, b_folder = fgs_folder)
;vxy_efi_ave = vxy_efi_ave_list(2:5,*)
;dataout_simple, save_folder+'/vxy_efi_ave'+suffix, vxy_efi_ave
;;; load
;vxy_efi_ave = datain_simple(save_folder+'/vxy_efi_ave'+suffix+'.dat', dim=4, type='double')
;vxy_efi_ave_x = vxy_efi_ave(0,*)

;;;;; DF propagation speed: average velocity inside the DF layer
;;;; get t_in and t_out
;;b_type = 'fgl'
;;b_folder = fgl_folder
;;seconds_check = 15.
;;b_in_list = event_status_instant(events, b_type, time_length = seconds_check/60., pre_time = -0.5*seconds_check/60., datafolder=b_folder, /max, comp_v = 2, time_interest = t_max_arr)
;;b_out_list = event_status_instant(events, b_type, vtrange = [time_double(events(0,*))-seconds_check,t_max_arr(1,*)], datafolder=b_folder, /min, comp_v = 2, time_interest = t_min_arr)
;;t_in = t_max_arr(1,*)
;;t_out = t_min_arr(1,*)
;;;;; get the values and save data
;;vperp_df_list = event_status_instant_efi(events, 'vperp', vtrange = [t_out, t_in], rpre_time = 2.5, rtime_length=1., e_folder=efs_folder, b_folder = fgs_folder)
;;vperp_df = vperp_df_list(2:5,*)
;;dataout_simple, save_folder+'/df_vperp'+suffix, vperp_df
;;vxy_df_list = event_status_instant_efi(events, 'vxy', vtrange = [t_out, t_in], rpre_time = 2.5, rtime_length=1., e_folder=efs_folder, b_folder = fgs_folder)
;;vxy_df = vxy_df_list(2:5,*)
;;dataout_simple, save_folder+'/df_vxy'+suffix, vxy_df
;;;;; load
;vperp_df = datain_simple(save_folder+'/df_vperp'+suffix+'.dat', dim=4, type='double')
;vxy_df = datain_simple(save_folder+'/df_vxy'+suffix+'.dat', dim=4, type='double')
;vperp_df_x = vperp_df(0,*)
;vperp_df_ttl = vperp_df(3,*)
;vxy_df_x = vxy_df(0,*)
;vxy_df_ttl = vxy_df(3,*)


;;;; current strength
;;; save
;;df_i = df_current_dens(events, /linear, b_in = b_in, b_out = b_out, bfolder = fgl_folder)
;;dataout_simple, 'df_i'+suffix, df_i
;;;; load
;df_i = datain_simple(save_folder+'/df_i'+suffix+'.dat', dim=1, type='double')

;;;;;; magnetic field value
;bq_list = event_status_instant(events, 'fgs', pre_time=2.5, time_length=1., datafolder=fgs_folder)
;bqz = bq_list(4,*)
;print, median(bqz, /even)
;stop
;;;;;; peak
;;;; 40s
;bpeak_list = event_status_instant(events, 'fgs', pre_time=-0.25, time_length=0.5, datafolder=fgs_folder, /max)
;dur_suf = '_40s'
;;; 2min
;;bpeak_list = event_status_instant(events, 'fgs', pre_time=-1., time_length=2., datafolder=fgs_folder, /max)
;bz_peak = bpeak_list(4,*)-bqz 
;;dataout_simple, save_folder+'/dfb_dbz_peak'+dur_suf+suffix, bz_peak

;;;; flux transport
;intEy_dsl_list = event_status_instant_efi(events, 'intEy_dsl', pre_time = -0.5, time_length=4., rpre_time = 2.5, rtime_length=1., /inc, e_folder=efs_folder, /reverse_bc)
;intEy_dsl = intEy_dsl_list(2,*)

;;;;; all types of flux transport
;fluxtypes = ['vxbz', 'vzbx']
;;; 40 seconds
;flux = bflux_cal([time_double(events(0,*)), time_double(events(0,*))+40.], probe_list = events(3, *), fluxtypes = fluxtypes, v_folder = vi_folder, fgs_folder = fgs_folder, peak_ey = peak_ey) 
;dur_suf = '_40s'
;;; 2 min
;;flux = bflux_cal([time_double(events(0,*)), time_double(events(0,*))+120.], probe_list = events(3, *), fluxtypes = fluxtypes, v_folder = vi_folder, fgs_folder = fgs_folder, peak_ey = peak_ey) 
;vxbz_peak = peak_ey(0,*)
;vzbx_peak = peak_ey(1,*)
;flux_vxbz = flux(0,*)
;flux_vzbx = flux(1,*)
;vxbz_o_vzbx = abs(flux_vxbz/flux_vzbx)
;dataout_simple, save_folder+'/dfb_vxbz_peak'+dur_suf+suffix, vxbz_peak
;dataout_simple, save_folder+'/dfb_vzbx_peak'+dur_suf+suffix, vzbx_peak

;;;;; peak Ey DSL
;;;;;;;;; 40s
;Epeak_dsl_list = event_status_instant_efi(events, 'efs_dsl', pre_time = -1./3., time_length=2./3., rpre_time = 2.5, rtime_length=1., /peak, e_folder=efs_folder, /reverse_bc)
;dur_suf = '_40s'
;;;;; use vexb for offset removal
;;;Epeak_dsl_list = event_status_instant_efi(events, 'efs_dsl', pre_time = -0.25, time_length=0.5, rpre_time = 2.5, rtime_length=1., /peak, e_folder=efs_folder, /reverse_bc, vexb_dsl_folder = vexb_dsl_folder)
;;;;; 2min
;;;Epeak_dsl_list = event_status_instant_efi(events, 'efs_dsl', pre_time = -1., time_length=2., rpre_time = 2.5, rtime_length=1., /peak, e_folder=efs_folder, /reverse_bc)
;;;dur_suf = '2min'
;;;;; use vexb for offset removal
;;;Epeak_dsl_list = event_status_instant_efi(events, 'efs_dsl', pre_time = -1., time_length=2., rpre_time = 2.5, rtime_length=1., /peak, e_folder=efs_folder, /reverse_bc, vexb_dsl_folder = vexb_dsl_folder)
;Ey_peak_dsl = Epeak_dsl_list(3,*)
;;dataout_simple, save_folder+'/dfb_ey_peak'+dur_suf+suffix, Ey_peak_dsl

;;;; average Ey DSL
;;; 30s
;Eave_dsl_list = event_status_instant_efi(events, 'efs_dsl', pre_time = -0.25, time_length=0.5, rpre_time = 2.5, rtime_length=1., e_folder=efs_folder, /reverse_bc)
;Eave_dsl_list = event_status_instant_efi(events, 'efs_dsl', pre_time = -0.25, time_length=0.5, rpre_time = 2.5, rtime_length=1., e_folder=efs_folder, /reverse_bc, vexb_dsl_folder = vexb_dsl_folder)
;; 2min
;Eave_dsl_list = event_status_instant_efi(events, 'efs_dsl', pre_time = -1., time_length=2., rpre_time = 2.5, rtime_length=1., e_folder=efs_folder, /reverse_bc)
;Eave_dsl_list = event_status_instant_efi(events, 'efs_dsl', pre_time = -1., time_length=2., rpre_time = 2.5, rtime_length=1., e_folder=efs_folder, /reverse_bc, vexb_dsl_folder = vexb_dsl_folder)
;;; use x duration
;dur_suf = '_xdur'
;durations = dfb_duration(events, method = 'x', datafolder = pos_folder)
;vtrange = [time_double(events(0,*)), time_double(events(0,*))+durations]
;Eave_dsl_list = event_status_instant_efi(events, 'efs_dsl', vtrange = vtrange, rpre_time = 2.5, rtime_length=1., e_folder=efs_folder, /reverse_bc)
;Ey_ave_dsl = Eave_dsl_list(3,*)
;dataout_simple, save_folder+'/dfb_ey_ave'+dur_suf+suffix, Ey_ave_dsl
;dataout_simple, save_folder+'/dfb_duration'+suffix, durations
;stop

;;;; length in X
dur_suf = '_xdur'
durations = dfb_duration(events, method = 'x', datafolder = pos_folder)
vtrange = [time_double(events(0,*)), time_double(events(0,*))+durations]
vx_ave_list = event_status_instant(events, 'vi', vtrange = vtrange, datafolder=vi_folder)
vx_ave = vx_ave_list(2,*)
lengths = vx_ave*durations
dataout_simple, save_folder+'/vx_ave'+dur_suf+suffix, vx_ave
dataout_simple, save_folder+'/dfb_length'+dur_suf+suffix, lengths
stop

;; peak vexb
;vexb_peak_list = event_status_instant(events, 'vexb', pre_time=-0.25, time_length=0.5, datafolder=vexb_folder, /peak)
;vexb_peak_y = vexb_peak_list(3,*)

;;; average vexb
;vexb_ave_list = event_status_instant(events, 'vexb', pre_time=-0.25, time_length=0.5, datafolder=vexb_folder)
;vexb_ave_y = vexb_ave_list(3,*)

;;;;;;;;;;;;;;;;;;;; make traditional statistical plot

title = ''

;;;;;;;;;; qtt bin

;qtt_bin = x
;qtt_title = 'X [RE]'
;;qtt_range = [-12, -6]
;;qtt_range = [-25, -6]
;binsize = 1.5
;;bin_boundaries = [-25,-22,-20,-18,-16,-14,-12,-11,-10,-9,-8,-7,-6]
;;bin_boundaries = [-25,-22,-20,-17,-14,-11.5,-10.5,-9.5,-8.5,-7.5,-6]
;bin_boundaries = [-20., -16., -12., -9., -6.]

;qtt_bin = ny
;qtt_title = 'ny'
;vertical = 1
;binsize = 0.25

;qtt_bin = h
;qtt_title = 'h'
;vertical = 1
;binsize = 0.25

;qtt_bin = omni_bz_peak
;qtt_title = 'OMNI Bz peak [nT]'
;binsize = 1

;qtt_bin = omni_bz_ave
;qtt_title = 'OMNI Bz average [nT]'
;binsize = 2

;qtt_bin = omni_Pdyn_peak
;qtt_title = 'OMNI Pdyn peak [nPa]'
;binsize = 0.5

;qtt_bin = omni_Pdyn_ave
;qtt_title = 'OMNI Pdyn average [nPa]'
;binsize = 0.5

;qtt_bin = omni_vxb_y_peak
;qtt_title = 'OMNI Ey vxb peak [mV/m]'
;binsize = 1.

;qtt_bin = omni_vxb_y_ave
;qtt_title = 'OMNI Ey vxb average [mV/m]'
;binsize = 1.

;qtt_bin = omni_bz_peak*omni_Pdyn_peak*1e12
;qtt_title = 'OMNI Pdyn x Bz peak [nPa*nT]'
;binsize = 1.5
;qtt_range = [-20, 20]

;qtt_bin = omni_bz_ave*omni_Pdyn_ave*1e12
;qtt_title = 'OMNI Pdyn x Bz ave [nPa*nT]'
;binsize = 1.5
;qtt_range = [-20, 20]

;qtt_bin = kyoto_ae_peak
;qtt_title = 'kyoto AE peak [nT]'
;binsize = 100

;qtt_bin = kyoto_ae_ave
;qtt_title = 'kyoto AE average [nT]'
;binsize = 100

;qtt_bin = pseudo_ae_peak
;qtt_title = 'pseudo AE peak [nT]'
;binsize = 80
;k_c = 5

qtt_bin = pseudo_ae_ave
qtt_title = 'pseudo AE average [nT]'
;binsize = 100
bin_boundaries = [0, 100, 200, 400, 10000000.]
;bin_boundaries = [0, 100, 200, 350, 10000000.]
k_c = 5

;qtt_bin = v_peak 
;qtt_title = 'Vpeak [km/s]'

;qtt_bin = vperp_peak_x
;qtt_title = 'Vperpx_peak [km/s]'
;binsize = 80

;qtt_bin = df_i
;qtt_title = 'i [nA/m2]'

;;;;;;;;;;;;; qtt 2 bin

qtt_2_bin = x
qtt_2_title = 'X [RE]'
;qtt_2_range = [-13, -7]

;qtt_2_bin = pseudo_ae_ave
;qtt_2_title = 'pseudo AE average [nT]'

;qtt_2_bin = abs(h)
;qtt_2_title = '|h|'

;qtt_2_bin = df_i
;qtt_2_title = 'i [nA/m2]'
;qtt_2_range = [0,3e7]
;k_c = 5

;qtt_2_bin = bz_peak
;qtt_2_title = 'Bz peak [nT]'
;qtt_2_range = [5,15]

;k_c = 5
;qtt_2_bin = intEy_dsl
;qtt_2_title = 'Fluxtran [nWb/m]'

;qtt_2_bin = v_ave
;qtt_2_title = 'Vave [km/s]'
;
;qtt_2_bin = v_peak
;qtt_2_title = 'Vpeak [km/s]'

;qtt_2_bin = vxy_peak_x
;qtt_2_title = 'inferred vxy_x, peak [km/s]'
;pm = 1
;k_c = 5
;qtt_2_range = [-200,400]

;qtt_2_bin = vxy_efi_peak_x
;qtt_2_title = 'inferred vxy_x from EFI, peak [km/s]'
;pm = 1
;k_c = 5
;qtt_2_range = [-500,1000]

;qtt_2_bin = v_peak_x
;qtt_2_title = 'Vx_peak [km/s]'
;pm = 1
;k_c = 5
;;qtt_2_range = [-200,400]

;qtt_2_bin = v_peak_z
;title = 'Vz peak'
;ratio_pm = 2
;qtt_2_title = '#_negative/#_positive'

;;qtt_2_range = [-200,400]
;qtt_2_bin = vperp_peak_x
;qtt_2_title = 'Vperpx_peak [km/s]'
;pm = 1
;k_c = 5
;qtt_2_range = [-150,250]

;qtt_2_bin = vperp_peak_ttl
;qtt_2_title = 'Vperp_ttl_peak [km/s]'
;k_c = 5
;qtt_2_range = [0,300]

;qtt_2_bin = vperp_ave_ttl
;qtt_2_title = 'Vperp_ttl_ave [km/s]'
;k_c = 5
;qtt_2_range = [0,200]

;qtt_2_bin = vperp_efi_peak_x
;qtt_2_title = 'Vperp_efi_x_peak [km/s]'
;pm = 1
;k_c = 5
;qtt_2_range = [-300,400]

;qtt_2_bin = vperp_efi_peak_ttl
;qtt_2_title = 'Vperp_efi_ttl_peak [km/s]'
;k_c = 5
;qtt_2_range = [0,600]

;qtt_2_bin = vperp_efi_ave_ttl
;qtt_2_title = 'Vperp_efi_ttl_ave [km/s]'
;k_c = 5
;qtt_2_range = [0,400]

;qtt_2_bin = vperp_df_x
;qtt_2_title = 'V_DF_perp_x [km/s]'
;pm = 1
;k_c = 5
;qtt_2_range = [-150,250]

;qtt_2_bin = vperp_df_ttl
;qtt_2_title = 'V_DF_perp_ttl [km/s]'
;k_c = 5
;qtt_2_range = [0,300]

;qtt_2_bin = vxy_df_x
;qtt_2_title = 'V_DF_XY_x [km/s]'
;pm = 1
;k_c = 5
;qtt_2_range = [-200,400]

;;;;;; vxbz vs vzbx ;;;;;;;
;qtt_2_bin = vxbz_o_vzbx
;qtt_2_title = '|VxBz|/|VzBx|'
;qtt_2_range = [0.,100.]

;qtt_2_bin = flux_vxbz
;qtt_2_title = 'flux of VxBz for 2 min [mWb/m]'
;pm = 1
;qtt_2_range = [-100.,200.]

;qtt_2_bin = flux_vzbx
;qtt_2_title = 'flux of VzBx for 2 min [mWb/m]'
;pm = 1
;;qtt_2_range = [0.,100.]

;qtt_2_bin = vxbz_peak
;qtt_2_title = 'peak VxBz [mV/m]'
;pm = 1
;qtt_2_range = [-5.,10.]

;qtt_2_bin = vzbx_peak
;qtt_2_title = 'peak VzBx [mV/m]'
;pm = 1
;qtt_2_range = [-5.,10.]

;;; check when it's positive when it's negative

;qtt_2_bin = flux_vzbx
;title = 'flux of vzbx'
;;ratio_pm = 1
;;qtt_2_title = '#_positive/#_negative'
;ratio_pm = 2
;qtt_2_title = '#_negative/#_positive'

;qtt_2_bin = vzbx_peak
;title = 'peak VzBx'
;ratio_pm = 2
;qtt_2_title = '#_negative/#_positive'

;qtt_2_bin = Ey_peak_dsl
;title = 'Ey_dsl peak'
;qtt_2_title = 'Ey_peak_dsl [mV/m]'
;qtt_2_range = [-10,15]
;pm = 1
;k_c = 5
;;;;; when positive when negative
;;ratio_pm = 2
;;qtt_2_title = '#_negative/#_positive'
                
;qtt_2_bin = Ey_ave_dsl
;qtt_2_title = 'Ey_ave_dsl [mV/m]'
;qtt_2_range = [-5,10]
;pm = 1
;k_c = 5

;qtt_2_bin = vexb_peak_y
;qtt_2_title = 'vexb_y_peak [mV/m]'
;qtt_2_range = [-20,20]

;qtt_2_bin = vexb_ave_y
;qtt_2_title = 'vexb_y_ave [mV/m]'
;qtt_2_range = [-20,20]

;;;;;;; selective of events
;;; ny sepe
;;i_bin = where(abs(ny) lt 0.3)
;;; x separate
;i_bin = where((x gt -10) and (x lt -9))
;
;qtt_bin = qtt_bin(*, i_bin)
;qtt_2_bin = qtt_2_bin(*, i_bin)

;;;;;;;;;;;;;;; plot the bin medians and means of quantities ;;;;;;;;;;;;;
qtt_bin = transpose(qtt_bin)
qtt_2_bin = transpose(qtt_2_bin)
stat_plot, qtt_bin, qtt_2_bin, k_c = k_c, bin_range = bin_range, binsize = binsize, qtt_2_range = qtt_2_range, qtt_range = qtt_range, qtt_2_title = qtt_2_title, qtt_title = qtt_title, kinbin = kinbin, bincntrs_out = bincenters, pm = pm, ratio_pm = ratio_pm, vertical = vertical, bin_boundaries = bin_boundaries, title = 'DFB '+title, avrg = avrg, std = std, med = med
print, med

;;;;;;;;;;;;;;;;;;;; make circle plot ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; change here for quantities to plot ;;;;;
;;; axis
;dim1 = x
;dim1_title = 'X [RE]'
;dim2 = y
;dim2_title = 'Y [RE]'
;;; quantity as circle
;circ_q = v_peak
;circ_title = 'V_peak [km/s]'
;circ_legend = [100, 200]
;;; quantity as color
;color_q = df_i
;color_title = 'i [nA/m]'
;;;; make the plot
;quantity_loc_plot, dim1, dim2, circ_q = circ_q, color_q = color_q, qtt_2_title = dim1_title, qtt_title = dim2_title, circ_title = circ_title, color_title = color_title, circ_legend = circ_legend, /isotropic

stop
end
