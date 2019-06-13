pro main_superpose
;; do superposed epoch analysis
thm_init
;; generate bin plots, like occurence rates at different locations, etc.
computer='I:'
;computer='/home/jliu'

@declare

split = 'x'
;split = 'h'

quantity = 'bz'
;quantity = 'ey_dsl'
;quantity = 'intEy_dsl'
;quantity = 'int_vxbz'
;quantity = 'int_vzbx'

;btype = 'fgl'
btype = 'fgs'

t_type = 't0'
;t_type = 'tout'
;t_type = 'tin'

case split of
'h': listname = 'dfb_list_lead_tail.txt'
'x': listname = 'dfb_list_lead_tail_fgs.txt'
;'x': listname = 'dfb_list_lead_tail.txt'
endcase

case btype of
'fgl': b_folder = fgl_folder
'fgs': b_folder = fgs_folder
endcase
events = load_list(listname, folder = listfolder)
;;; for test
;events = events(*,0:50)

;; for smoothing
;smooth_width = 1
smooth_width = 3

;; for average
time_length = 1.
pre_time = 2.5

;; for after value
aft_range = [2., 3.]

;; for superpose
pretime = 3.
afttime = 3.

case split of
;'x': locs = [-20., -15., -12., -8., -6.] ;; last one not work
'x': locs = [-20., -16., -12., -9., -6.] ;; works, use this for paper
;'x': locs = [-20., -15., -10., -6.] ;; works
;'x': locs = [-60., -1.] ;; for all distances
'h': locs = [0., 0.2, 0.4, 0.6, 1.]
endcase

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
if strcmp(split, 'h') then begin
	h = scale_height(events, b_folder = fgs_folder, p_folder = Pall_folder)
	loc_qtt = abs(h)
endif

;;; get quiet time values for detrending
if strcmp(quantity, 'bz') then begin
	b_quiet = event_status_instant(events, 'fgs', time_length = time_length, pre_time = pre_time, datafolder=fgs_folder)
	bz_q = b_quiet(4,*)
	bt_q = b_quiet(5,*)
endif

;; begin superpose
after_v = dblarr(n_elements(locs)-1)

for i = 0, n_elements(locs)-2 do begin
	j_this = where((loc_qtt gt locs(i)) and (loc_qtt lt locs(i+1)))
	t0p_list_this = t0p_list(*, j_this)
	t0_this = time_double(t0p_list_this(0,*))
	if strcmp(quantity, 'bz') then begin
		prefix = btype+'_gsm_z'
		bz_q_this = bz_q(*, j_this)
		superpose_data_fast, btype, components='z', t0p_list=t0p_list_this, pre=3., after=3., datafolder=b_folder, detrend_v = bz_q_this, smooth_width = smooth_width
	endif
	if strcmp(quantity, 'ey_dsl') then begin
		prefix = 'efs_dsl_y'
		superpose_data_fast_efi, 'efs_dsl', components='y', t0p_list=t0p_list_this, pre=3., after=3., e_folder=efs_folder, rtrange = [t0_this-180., t0_this-120.], smooth_width = smooth_width, /reverse_bc
	endif
	if strcmp(quantity, 'intEy_dsl') then begin
		prefix = quantity
		superpose_data_fast_efi, quantity, t0p_list=t0p_list_this, pre=3., after=3., e_folder=efs_folder, rtrange = [t0_this-180., t0_this-120.], smooth_width = smooth_width, /reverse_bc
	endif
	if strmatch(quantity, 'int_v*') then begin
		prefix = strmid(quantity, 4)
		superpose_data_fast_vxb, prefix, t0p_list=t0p_list_this, pre=3., after=3., v_folder=vi_folder, fgs_folder=fgs_folder, smooth_width = smooth_width, /int
	endif
	split_vec, prefix+'_superpose_quartiles'
	copy_data, prefix+'_superpose_quartiles_y', prefix+'_superpose_quartiles_dtd_med_'+strcompress(string(i), /remove)
	time_clip, prefix+'_superpose_quartiles_y', time_double('2000 1 1')+aft_range(0)*60., time_double('2000 1 1')+aft_range(1)*60.
	get_data,  prefix+'_superpose_quartiles_y_tclip', data = after
	after_v(i) = mean(after.y, /nan)
endfor

tplot, prefix+'_superpose_quartiles_dtd_med_*'
for i = 0, n_elements(after_v)-1 do begin
	timebar, after_v(i), varname = prefix+'_superpose_quartiles_dtd_med_'+strcompress(string(i), /remove), /databar
endfor
timebar, '2000 1 1'
makepng, quantity+'_'+split+'_'+t_type

stop
end
