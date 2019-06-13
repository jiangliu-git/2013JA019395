pro main_relations_bbf
;;; examine all kinds of relations (flow, field strength, etc)

thm_init

;;; the event
computer = 'I:'
;computer = '/home/jliu'

@declare

listname = 'bbf_ps_list_original2.txt'
;listname = 'dfb_bbf_ps_list_original2.txt' ;; the embeded dfbs (for future)
;listname = 'bbf_ps_list_lead.txt'

events_all = load_list(listname, folder = listfolder)
events = events_all
;events = events_all(*,30:40)

bbf_ranges = time_double(events(0:1, *))

;;;;;;;;;; run
;;;;;;;; get the event location and seperate
;pos_list = event_status_instant(events, 'pos', vtrange = bbf_ranges, datafolder=pos_folder)
;x = pos_list(2,*)/6371.
;y = pos_list(3,*)/6371.
;z = pos_list(4,*)/6371.

;;;;;;;; get event scale height and Blobeq
h = scale_height(events, b_folder = fgs_folder, p_folder = Pall_folder, time_length = 1., pre_time = 2.5, Blobe = Blobe)

;dataout_simple, 'bbf_h', h
;;;; load
;h = datain_simple(save_folder+'/bbf_h'+'.dat', dim=1, type='double')

;;;;;;;;;;;;;;;  ground signatures
;;;;;;;; get ae index
;kyoto_ae_peak_list = event_status_instant(events, 'kyoto_ae', pre_time = 0, time_length=3., datafolder=kyoto_ae_folder, /peak)
;kyoto_ae_peak = kyoto_ae_peak_list(2,*)

;;;;;;;; get ae index
;kyoto_ae_ave_list = event_status_instant(events, 'kyoto_ae', pre_time = 0, time_length=3., datafolder=kyoto_ae_folder)
;kyoto_ae_ave = kyoto_ae_ave_list(2,*)

;;;;;;;; get ae index
;pseudo_ae_peak_list = event_status_instant(events, 'pseudo_ae', pre_time = -1.5, time_length=3., datafolder=pseudo_ae_folder, /peak)
;pseudo_ae_peak = pseudo_ae_peak_list(2,*)

;;;;;;;; get ae index
;pseudo_ae_ave_list = event_status_instant(events, 'pseudo_ae', pre_time = -1.5, time_length=3., datafolder=pseudo_ae_folder)
;pseudo_ae_ave = pseudo_ae_ave_list(2,*)

;;;;;;;; get ae index
;pseudo_ae_ave_list = event_status_instant(events, 'pseudo_ae', pre_time = -1.5, time_length=3., datafolder=pseudo_ae_folder)
;pseudo_ae_ave = pseudo_ae_ave_list(2,*)

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
;;; peak velocity
;v_peak_list = event_status_instant(events, 'vi', pre_time = -0.25, time_length=0.5, datafolder=vi_folder, /peak)
;v_peak_x = v_peak_list(2,*)
;v_peak_y = v_peak_list(3,*)
;v_peak_z = v_peak_list(4,*)
;v_peak = v_peak_list(5,*)

;;; average velocity
;v_ave_list = event_status_instant(events, 'vi', time_length=4, datafolder=vi_folder)
;v_ave_x = v_ave_list(2,*)
;v_ave_y = v_ave_list(3,*)
;v_ave_z = v_ave_list(4,*)
;v_ave = v_ave_list(5,*)

;;; quiet time average velocity
v_ave_list = event_status_instant(events, 'vi', time_length=1, pre_time = 2.5, datafolder=vi_folder)
v_ave_x = v_ave_list(2,*)
v_ave_y = v_ave_list(3,*)
v_ave_z = v_ave_list(4,*)
v_ave = v_ave_list(5,*)

v_std_list = event_status_instant(events, 'vi', time_length=1, pre_time = 2.5, datafolder=vi_folder, /std)
v_std_x = v_std_list(2,*)
v_std_y = v_std_list(3,*)
v_std_z = v_std_list(4,*)
v_std = v_std_list(5,*)

stop

;;;; perp velocity
;;;; peak velocity
;vperp_peak_list = event_status_instant(events, 'viperp', pre_time = -0.25, time_length=0.5, datafolder=viperp_folder, /peak)
;vperp_peak_x = vperp_peak_list(2,*)
;vperp_peak_y = vperp_peak_list(3,*)
;vperp_peak_z = vperp_peak_list(4,*)
;vperp_peak = vperp_peak_list(5,*)

;;; average perp velocity
;vperp_ave_list = event_status_instant(events, 'viperp', time_length=4, datafolder=viperp_folder)
;vperp_ave_x = vperp_ave_list(2,*)
;vperp_ave_y = vperp_ave_list(3,*)
;vperp_ave_z = vperp_ave_list(4,*)
;vperp_ave = vperp_ave_list(5,*)

;;; quiet time average perp velocity
;vperp_ave_list = event_status_instant(events, 'viperp', pre_time = 2.5, time_length=1, datafolder=viperp_folder)
;vperp_ave_x = vperp_ave_list(2,*)
;vperp_ave_y = vperp_ave_list(3,*)
;vperp_ave_z = vperp_ave_list(4,*)
;vperp_ave = vperp_ave_list(5,*)
;
;;;; quiet time perp velocity std
;vperp_std_list = event_status_instant(events, 'viperp', pre_time = 2.5, time_length=1, datafolder=viperp_folder, /std)
;vperp_std_x = vperp_std_list(2,*)
;vperp_std_y = vperp_std_list(3,*)
;vperp_std_z = vperp_std_list(4,*)
;vperp_std = vperp_std_list(5,*)
;
;stop

;;;; flux transport
;intEy_dsl_list = event_status_instant_efi(events, 'intEy_dsl', pre_time = -0.5, time_length=4., rpre_time = 2.5, rtime_length=1., /inc, e_folder=efs_folder, /reverse_bc)
;intEy_dsl = intEy_dsl_list(2,*)

;;;; peak Ey DSL
;Epeak_dsl_list = event_status_instant_efi(events, 'efs_dsl', vtrange = bbf_ranges, rpre_time = 2.5, rtime_length=1., /peak, e_folder=efs_folder, /reverse_bc)
;;;; use vexb for offset removal
;;Epeak_dsl_list = event_status_instant_efi(events, 'efs_dsl', vtrange = bbf_ranges, rpre_time = 2.5, rtime_length=1., /peak, e_folder=efs_folder, /reverse_bc, vexb_dsl_folder = vexb_dsl_folder)
;Ey_peak_dsl = Epeak_dsl_list(3,*)

;;;; all types
fluxtypes = ['vxbz', 'vzbx']
flux = bflux_cal(bbf_ranges, probe_list = events(3, *), fluxtypes = fluxtypes, v_folder = vi_folder, fgs_folder = fgs_folder) 
flux_vxbz = flux(0,*)
flux_vzbx = flux(1,*)
vxbz_o_vzbx = abs(flux_vxbz/flux_vzbx)

;;;; average Ey DSL
;;;;; 30s
;Eave_dsl_list = event_status_instant_efi(events, 'efs_dsl', vtrange = bbf_ranges, rpre_time = 2.5, rtime_length=1., e_folder=efs_folder, /reverse_bc)
;;Eave_dsl_list = event_status_instant_efi(events, 'efs_dsl', vtrange = bbf_ranges, rpre_time = 2.5, rtime_length=1., e_folder=efs_folder, /reverse_bc, vexb_dsl_folder = vexb_dsl_folder)
;Ey_ave_dsl = Eave_dsl_list(3,*)

;; peak vexb
;vexb_peak_list = event_status_instant(events, 'vexb', pre_time=-0.25, time_length=0.5, datafolder=vexb_folder, /peak)
;vexb_peak_y = vexb_peak_list(3,*)

;;; average vexb
;vexb_ave_list = event_status_instant(events, 'vexb', pre_time=-0.25, time_length=0.5, datafolder=vexb_folder)
;vexb_ave_y = vexb_ave_list(3,*)
;;;;;;;;;;;;;;;;;;;; make traditional statistical plot

;;;;;;;;;; qtt bin

;qtt_bin = x
;qtt_title = 'X [RE]'
;;qtt_range = [-12, -6]
;qtt_range = [-25, -6]
;binsize = 2.
;;bin_boundaries = [-25,-22,-20,-18,-16,-14,-12,-11,-10,-9,-8,-7,-6]

;qtt_bin = ny
;qtt_title = 'ny'
;vertical = 1
;binsize = 0.25

qtt_bin = h
qtt_title = 'h'
vertical = 1
binsize = 0.2

;qtt_bin = kyoto_ae_peak
;qtt_title = 'kyoto AE peak [nT]'
;binsize = 100

;qtt_bin = kyoto_ae_ave
;qtt_title = 'kyoto AE average [nT]'
;binsize = 100

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

;qtt_bin = pseudo_ae_peak
;qtt_title = 'pseudo AE peak [nT]'
;binsize = 100

;qtt_bin = pseudo_ae_ave
;qtt_title = 'pseudo AE average [nT]'
;binsize = 100

;qtt_bin = v_peak 
;qtt_title = 'Vpeak [km/s]'

;qtt_bin = vperp_peak_x
;qtt_title = 'Vperpx_peak [km/s]'
;binsize = 80

;qtt_bin = df_i
;qtt_title = 'i [nA/m2]'

;;;;;;;;;;; qtt 2 bin

;qtt_2_bin = df_i
;qtt_2_title = 'i [nA/m2]'
;qtt_2_range = [0,3e7]
;k_c = 5

;qtt_2_bin = intEy_dsl
;qtt_2_title = 'Fluxtran [nWb/m]'

;qtt_2_bin = v_ave
;qtt_2_title = 'Vave [km/s]'
;
;qtt_2_bin = v_peak
;qtt_2_title = 'Vpeak [km/s]'

;qtt_2_bin = v_peak_x
;qtt_2_title = 'Vx_peak [km/s]'
;pm = 1
;k_c = 5

;qtt_2_bin = vperp_peak_x
;qtt_2_title = 'Vperpx_peak [km/s]'
;pm = 1
;k_c = 5

qtt_2_bin = vxbz_o_vzbx
qtt_2_title = '|VxBz|/|VzBx|'
qtt_2_range = [0.,100.]

;qtt_2_bin = Ey_peak_dsl
;qtt_2_title = 'Ey_peak_dsl [mV/m]'
;qtt_2_range = [-20,30]
;pm = 1
;k_c = 5
                
;qtt_2_bin = Ey_ave_dsl
;qtt_2_title = 'Ey_ave_dsl [mV/m]'
;qtt_2_range = [-5,5]
;pm = 1
;k_c = 5

;qtt_2_bin = vexb_peak_y
;qtt_2_title = 'vexb_y_peak [mV/m]'
;qtt_2_range = [-20,20]

;qtt_2_bin = vexb_ave_y
;qtt_2_title = 'vexb_y_ave [mV/m]'
;qtt_2_range = [-20,20]

;;;;;;; selective of events
;i_bin = where(abs(ny) lt 0.3)
;qtt_bin = qtt_bin(*, i_bin)
;qtt_2_bin = qtt_2_bin(*, i_bin)

;;;;;;;;;;;;;; begin binning
qtt_bin = transpose(qtt_bin)
qtt_2_bin = transpose(qtt_2_bin)
stat_plot, qtt_bin, qtt_2_bin, k_c = k_c, bin_range = bin_range, binsize = binsize, qtt_2_range = qtt_2_range, qtt_range = qtt_range, qtt_2_title = qtt_2_title, qtt_title = qtt_title, kinbin = kinbin, bincntrs_out = bincenters, pm = pm, vertical = vertical, bin_boundaries = bin_boundaries, title = 'BBF'

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
