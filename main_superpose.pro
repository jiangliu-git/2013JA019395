pro main_superpose
;; do superposed epoch analysis
thm_init
;; generate bin plots, like occurence rates at different locations, etc.
;computer='I:'
computer='/home/jliu'

@declare

;split = 'x'
;split = 'y'
;split = 'h'
;split = 'ny'
split = 'ae'

;quantity = 'bz'
;quantity = 'dbz'
;quantity = 'vi'
;quantity = 'vx'
;quantity = 'vz'
;quantity = 'ey_dsl'
quantity = 'intEy_dsl'
;quantity = 'int_vxbz'
;quantity = 'int_vzbx'
;quantity = 'int_vxbz_dtd'
;quantity = 'int_vzbx_dtd'
;quantity = 'int_plasmaEy'
;quantity = 'int_plasmaEy_dtd'

;;; whether use only positive/negative transport events
sign_use_suf = ''
;sign_use_suf = '_positive'
;sign_use_suf = '_negative'

;;; choose seasons to avoid solar minimum
season_suf = ''
;season_suf = '_0809'

;;;; choose whether to use only 9-10 RE events
x_limit_suf = '_9to10'
;x_limit_suf = ''

;;;; choose whether to use 0-14UT only 
if strcmp(split, 'ae') then ut_suf = '_014' else ut_suf = '' ;; 0-14 UT only

;;;; choose whether to follow full plasma Ey for positive value 
if strmatch(quantity, 'int_v*') then vxb_suf = '_followfull' else vxb_suf = ''
vxb_suf = ''

if strcmp_or(quantity, ['vi']) then sign_use_suf = ''

;;; whether to use vexb for removeing electric field offset or 0 to remove
dtrd_suf = ''
;dtrd_suf = '_vexb'

if ~strcmp_or(quantity, ['ey_dsl', 'int_plasmaEy']) then dtrd_suf= ''
if strcmp(dtrd_suf, 'vexb') then edtrd_folder = vexb_dsl_folder else edtrd_folder = ''

;;; choose the plot type
plot_type = 'together'
;plot_type = 'split'

t_type = 't0'
;t_type = 'tout'
;t_type = 'tin'

;; for after value
aft_range = [2., 3.]

;; for average
time_length = 1.
pre_time = 2.5

;; for superpose range
pretime = 3.
afttime = 3.

;; x and h locations
case split of
'x': locs = [-20., -16., -12., -9., -6.] ;; works, use this for paper
;'x': locs = [-60., -1.]  ;; for all distances
;'y': locs = [-12., 0., 12.]
'y': locs = [-12., -4., 0., 4., 12.]
'h': locs = [0., 0.2, 0.4, 0.6, 1.]
;'h': locs = [0., 1.] ;; for all
'ny': locs = [0., 0.25, 0.5, 0.75, 1.]
'ae': locs = [0., 100, 200, 400, 1000000000.]
;'ae': locs = [0., 100, 200, 350, 1000000000.]
endcase

;; for smoothing
if n_elements(locs) le 2 then smooth_width = 1 else smooth_width = 3

;btype = 'fgl'
btype = 'fgs'

list_suf = ''
if strcmp_or(split, ['x', 'ae', 'y']) and strcmp(btype, 'fgs') and (strcmp(quantity, 'b', 1) or strcmp(quantity, 'db', 2) or strmatch(quantity, '*dsl*')) then list_suf = '_fgs'
if strcmp(split, 'ny') then list_suf = '_earthward_df_binxbout'
listname = 'dfb_list_lead_tail'+list_suf+'.txt'

case btype of
'fgl': b_folder = fgl_folder
'fgs': b_folder = fgs_folder
endcase

;;;; load events
events = load_list(listname, folder = listfolder)
;;; for test
;events = events(*,100:110)

;;; get 0809 events only
if strcmp(season_suf, '_0809') then begin
	i_use = where(time_double(events(0,*)) lt time_double('2010 1 1'))
	events = events(*, i_use)
endif

;;; get 0-14UT events only
time = time_double(events(0,*))
date_str = strmid(events(0,*), 0, 10)
date = time_double(date_str)
ut=(time-date)/3600.
if strcmp(ut_suf, '_014') then begin
	i_good = where(ut le 14.)
	events = events(*, i_good)
endif

;;; try events at same X
if strcmp(x_limit_suf, '_9to10') then begin
	pos_list = event_status_instant(events, 'pos', time_length=3, datafolder=pos_folder)
	x = pos_list(2,*)/6371.
	i_good = where((x gt -10) and (x lt -9))
	events = events(*, i_good)
endif

;;;; split positive and negative events
if ~strcmp(sign_use_suf, '') then begin
	if strmatch(quantity, 'v[^i]') or strmatch(quantity, '*bz') then begin
		vi_peak_list = event_status_instant(events, 'vi', pre_time = -0.25, time_length=0.5, /peak, datafolder=vi_folder)
		vx_peak = vi_peak_list(2,*)
		vy_peak = vi_peak_list(3,*)
		vz_peak = vi_peak_list(4,*)
		case quantity of
			'bz': quantity_check = vx_peak
			'dbz': quantity_check = vx_peak
			'vx': quantity_check = vx_peak
			'vy': quantity_check = vy_peak
			'vz': quantity_check = vz_peak
		endcase
	endif
	if strcmp(quantity, 'ey_dsl') then begin
		Epeak_dsl_list = event_status_instant_efi(events, 'efs_dsl', pre_time = -0.25, time_length=0.5, rpre_time = 2.5, rtime_length=1., /peak, e_folder=efs_folder, /reverse_bc, vexb_dsl_folder = edtrd_folder)
		Ey_peak_dsl = Epeak_dsl_list(3,*)
		quantity_check = Ey_peak_dsl
	endif
	;;; for flux: the integration of t0 to t0+2min
	if strcmp(quantity, 'intEy_dsl') then begin
		intEy_dsl_list = event_status_instant_efi(events, 'intEy_dsl', pre_time = -1., time_length=2., rpre_time = 2.5, rtime_length=1., /inc, e_folder=efs_folder, /reverse_bc, vexb_dsl_folder = edtrd_folder)
		intEy_dsl = intEy_dsl_list(2,*)
		quantity_check = intEy_dsl
	endif
	if strmatch(quantity, 'int_v*') or strmatch(quantity, 'int_plasmaEy*') then begin
		if strmatch(quantity, 'int_v*') then prefix = strmid(quantity, 4, 4)
		if strmatch(quantity, 'int_plasma*') then prefix = 'vxbz-vzbx'
		if strmatch(quantity, '*dtd') then begin
			vxb_quiet_all = event_status_instant_vxb(events, prefix, time_length = time_length, pre_time = pre_time, fgs_folder = fgs_folder, v_folder=vi_folder)
			vxb_q_all = vxb_quiet_all(2,*)
		endif else begin
			vxb_q_all = 0
		endelse
		if strmatch(quantity, 'int_v*') and strcmp(vxb_suf, '_followfull') then prefix = 'vxbz-vzbx'
		int_vxb_list = event_status_instant_vxb(events, prefix, pre_time = -1., time_length=2., /int, fgs_folder = fgs_folder, v_folder=vi_folder, detrend_v = vxb_q_all)
		int_vxb = int_vxb_list(2,*)
		quantity_check = int_vxb
	endif
	;;;; split events
	i_pos = where(quantity_check gt 0)
	i_neg = where(quantity_check lt 0)
	case sign_use_suf of
	'_positive': i_use = i_pos
	'_negative': i_use = i_neg
	endcase
	if i_use(0) ne -1 then events = events(*, i_use) else begin
		print, 'No pos/neg events'
		stop
	endelse
endif

;;;;;;;;;;;;;;;;;;;;;;;; end of editabel part ;;;;;;;;;;;;;;;;;;
if strcmp(t_type, 't0') then begin
	t_use = events(0,*)
endif
if strcmp(t_type, 'tin') then begin
	b_in_list = event_status_instant(events, 'fgl', time_length = 15./60., pre_time = -0.5*15./60., datafolder=fgl_folder, /max, comp_v = 2, time_interest = t_max_arr)
	t_in = t_max_arr(1,*)
	t_use = time_string(t_in)
endif
if strcmp(t_type, 'tout') then begin
	b_in_list = event_status_instant(events, 'fgl', time_length = 15./60., pre_time = -0.5*15./60., datafolder=fgl_folder, /max, comp_v = 2, time_interest = t_max_arr)
	b_out_list = event_status_instant(events, 'fgl', vtrange = [time_double(events(0,*))-15.,t_max_arr(1,*)], datafolder=fgl_folder, /min, comp_v = 2, time_interest = t_min_arr)
	t_out = t_min_arr(1,*)
	t_use = time_string(t_out)
endif
t0p_list = [t_use, events(3, *)]

del_data, '*'

if strcmp(split, 'x') then begin
	pos_list = event_status_instant(events, 'pos', time_length=3, datafolder=pos_folder)
	x = pos_list(2,*)/6371.
	loc_qtt = x
endif
if strcmp(split, 'y') then begin
	pos_list = event_status_instant(events, 'pos', time_length=3, datafolder=pos_folder)
	y = pos_list(3,*)/6371.
	loc_qtt = y
endif
if strcmp(split, 'h') then begin
	h = scale_height(events, b_folder = fgs_folder, p_folder = Pall_folder)
	loc_qtt = abs(h)
endif
if strcmp(split, 'ny') then begin
	n = df_normal(events, method=3, datatype='fgl', bfolder=fgl_folder, c_angle = 30)
	loc_qtt = abs(n(1,*))
endif
if strcmp(split, 'ae') then begin
	pseudo_ae_ave_list = event_status_instant(events, 'pseudo_ae', pre_time = -3, time_length=3., datafolder=pseudo_ae_folder)
	pseudo_ae_ave = pseudo_ae_ave_list(2,*)
	loc_qtt = pseudo_ae_ave
endif

;;; get quiet time values for detrending
if strcmp(quantity, 'dbz') then begin
	b_quiet = event_status_instant(events, 'fgs', time_length = time_length, pre_time = pre_time, datafolder=fgs_folder)
	bz_q = b_quiet(4,*)
	bt_q = b_quiet(5,*)
endif

if strmatch(quantity, 'int_v*dtd') or strcmp(quantity, 'int_plasmaEy_dtd') then begin
	if strmatch(quantity, 'int_v*dtd') then prefix = strmid(quantity, 4, 4)
	if strmatch(quantity, 'int_plasmaEy_dtd') then prefix = 'vxbz-vzbx'
	vxb_quiet = event_status_instant_vxb(events, prefix, time_length = time_length, pre_time = pre_time, fgs_folder = fgs_folder, v_folder=vi_folder)
	vxb_q = vxb_quiet(2,*)
endif

;; begin superpose
after_v = dblarr(n_elements(locs)-1)

loc_meds = dblarr(n_elements(locs))
for i = 0, n_elements(locs)-2 do begin
	j_this = where((loc_qtt gt locs(i)) and (loc_qtt lt locs(i+1)))
	loc_meds(i) = median(loc_qtt(j_this), /even)
	t0p_list_this = t0p_list(*, j_this)
	t0_this = time_double(t0p_list_this(0,*))
	if strcmp(quantity, 'dbz') then begin
		prefix = btype+'_gsm_z'
		bz_q_this = bz_q(*, j_this)
		superpose_data_fast, btype, components='z', t0p_list=t0p_list_this, pre=3., after=3., datafolder=b_folder, detrend_v = bz_q_this, smooth_width = smooth_width
	endif
	if strcmp(quantity, 'vi') then begin
		prefix = 'ptix_velocity_gsm_ttl'
		superpose_data_fast, quantity, /total, t0p_list=t0p_list_this, pre=3., after=3., datafolder=vi_folder, smooth_width = smooth_width
	endif
	if strcmp(quantity, 'vx') then begin
		prefix = 'ptix_velocity_gsm_x'
		superpose_data_fast, 'vi', components='x', t0p_list=t0p_list_this, pre=3., after=3., datafolder=vi_folder, smooth_width = smooth_width
	endif
	if strcmp(quantity, 'ey_dsl') then begin
		prefix = 'efs_dsl_y'
		superpose_data_fast_efi, 'efs_dsl', components='y', t0p_list=t0p_list_this, pre=3., after=3., e_folder=efs_folder, rtrange = [t0_this-180., t0_this-120.], smooth_width = smooth_width, /reverse_bc
	endif
	if strcmp(quantity, 'intEy_dsl') then begin
		prefix = quantity
		superpose_data_fast_efi, quantity, t0p_list=t0p_list_this, pre=3., after=3., e_folder=efs_folder, rtrange = [t0_this-180., t0_this-120.], smooth_width = smooth_width, /reverse_bc
	endif
	if strmatch(quantity, 'int_v*') or strmatch(quantity, 'int_plasmaEy*') then begin
		if strmatch(quantity, 'int_v*') then prefix = strmid(quantity, 4, 4)
		if strmatch(quantity, 'int_plasmaEy*') then prefix = 'vxbz-vzbx'
		if ~strmatch(quantity, '*dtd') then vxb_q = 0 else begin
			vxb_q_this = vxb_q(*, j_this)
		endelse
		superpose_data_fast_vxb, prefix, t0p_list=t0p_list_this, pre=3., after=3., v_folder=vi_folder, fgs_folder=fgs_folder, smooth_width = smooth_width, /int, detrend_v = vxb_q_this
	endif
	split_vec, prefix+'_superpose_quartiles'
	copy_data, prefix+'_superpose_quartiles_y', prefix+'_superpose_quartiles_med_'+strcompress(string(i), /remove)
	time_clip, prefix+'_superpose_quartiles_y', time_double('2000 1 1')+aft_range(0)*60., time_double('2000 1 1')+aft_range(1)*60.
	get_data,  prefix+'_superpose_quartiles_y_tclip', data = after
	after_v(i) = mean(after.y, /nan)
endfor

if strcmp(plot_type, 'split') then begin
	tplot, prefix+'_superpose_quartiles_med_*'
	for i = 0, n_elements(after_v)-1 do begin
		;timebar, after_v(i), varname = prefix+'_superpose_quartiles_med_'+strcompress(string(i), /remove), /databar
	endfor
endif

if strcmp(plot_type, 'together') then begin
	;;; make all medians together
	get_data, prefix+'_superpose_quartiles_med_0', data=data
	t = data.x
	n_pts = n_elements(data.y)
	n_dists = n_elements(after_v)
	median_arr = dblarr(n_pts, n_dists)
	;; make the color array
	color_dark = 20
	color_bright = 240
	color_arr_auto = floor(color_bright+dindgen(n_dists)/(n_dists-1)*(color_dark-color_bright))
	;; specific settings
	case n_dists of
		1: color_arr = 0
		4: color_arr = [40, 59, 78, 97]
		else: color_arr = color_arr_auto
	endcase
	for i = 0, n_dists-1 do begin
		get_data, prefix+'_superpose_quartiles_med_'+strcompress(string(i), /remove), data = data_this
		median_arr(*,i) = data_this.y
	endfor
	store_data, prefix+'_superpose_quartiles_med_all', data = {x:t, y:median_arr}
	options, prefix+'_superpose_quartiles_med_all', colors = color_arr
	tplot, prefix+'_superpose_quartiles_med_all', trange = ['1999 12 31 23 57', '2000 1 1 0 3']
	for i = 0, n_dists-1 do begin
		;timebar, after_v(i), /databar, varname = prefix+'_superpose_quartiles_med_all'
	endfor
	if n_dists lt 2 then split = ''
	dataout, prefix+'_superpose_quartiles_med_all', filename = save_folder+'/'+quantity+'_superpose_'+dtrd_suf+split+sign_use_suf+vxb_suf+season_suf+ut_suf+x_limit_suf
endif


timebar, '2000 1 1'
;makepng, quantity+'_'+split+'_'+t_type+sign_use_suf+'_'+plot_type
print, 'List used: '+'dfb_list_lead_tail'+list_suf+'.txt'

print, loc_meds

stop
end
