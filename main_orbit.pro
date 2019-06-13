pro main_orbit
;;; get the orbit of DFB-observing in different conditions.
thm_init
;; generate bin plots, like occurence rates at different locations, etc.
computer='I:'
;computer='/home/jliu'

@declare

listname = 'dfb_list_lead_tail_fgs.txt'
events = load_list(listname, folder = listfolder)
;events = events(*, 10:20) ;; for diagnose

xyz_ae = [0D,0,0,0,0.]

probes_less = ['a', 'd', 'e']
probes_more = ['a', 'b', 'c', 'd', 'e']

event_days = list_days(events(0,*), /compile, /double) ; get the days to begin the examination, (no repeat)
for i = 0, n_elements(event_days)-1 do begin
  day = event_days(i)
  if time_double(day) gt time_double('2010 1 1') then probes = probes_less else probes = probes_more
  del_data, 'th?_state_pos*'
  load_bin_data, trange = [day, day+24.*3600.], probes = probes, datatype = 'pos', datafolder = pos_folder, /tclip
  del_data, 'thg_idx_ae*'
  load_bin_data, trange = [day, day+24.*3600.], datatype = 'pseudo_ae', datafolder = pseudo_ae_folder, /tclip
  get_data, 'thg_idx_ae_tclip', data = data_ae
  ae = transpose(data_ae.y)
  ;;; for future
  ;del_data, 'th?_beta*'
  ;load_bin_data, trange = [day, day+24.*3600.], probes = probes, datatype = 'beta', datafolder = beta_folder, /tclip
  for j = 0, n_elements(probes)-1 do begin
    sc = probes(j)
    tinterpol_mxn, 'th'+sc+'_state_pos_tclip', 'thg_idx_ae_tclip', suff='_int'
    get_data, 'th'+sc+'_state_pos_tclip_int', data = data_pos
	t = transpose(data_pos.x)
	xyz = transpose(data_pos.y/RE)
	xyz_ae_this = [t, xyz, ae]
	;;; refine to range
	i_good = where((xyz(0,*) gt -30) and (xyz(0,*) lt -6) and (sqrt(xyz(1,*)^2+xyz(2,*)^2) lt 12), n_good)
	if n_good gt 0 then begin
		xyz_ae_this = xyz_ae_this(*, i_good)
		xyz_ae = [[xyz_ae], [xyz_ae_this]]
	endif
  endfor ; for of probes
endfor ; for of days

xyz_ae = xyz_ae(*,1:*)

dataout_simple, save_folder+'/orbit_xyz_ae_fgs', xyz_ae
stop

end
