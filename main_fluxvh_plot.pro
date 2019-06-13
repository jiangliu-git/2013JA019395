pro main_fluxvh_plot
;;;; plot dfb flux transport versus h, the paper's Figure 4.
;;;; savefiles come from main_relations_dfb.pro and main_superpose.pro
thm_init
@declare

;;;;; choose which detrend method to use
dtrd_suf = ''
;dtrd_suf = '_vexb'

if strcmp(dtrd_suf, '_vexb') then begin
	pref_e = ''
endif else begin 
	pref_e = delta_letter
endelse

suffix = ''

;;;;;;;;;;;;;;;;;;;;;;; the median/bar stat plots for ey peaks ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;dur_suf = '_40s'
;binsize_h = 0.23
;k_c = 5
;symsize = 0.8
;;;; load data
;h = datain_simple(save_folder+'/dfb_h'+suffix+'.dat', dim=1, type='double')
;vxbz_peak = datain_simple(save_folder+'/dfb_vxbz_peak'+dur_suf+suffix+'.dat', dim=1, type='double')
;vzbx_peak = datain_simple(save_folder+'/dfb_vzbx_peak'+dur_suf+suffix+'.dat', dim=1, type='double')
;Ey_peak = datain_simple(save_folder+'/dfb_ey_peak'+dur_suf+suffix+dtrd_suf+'.dat', dim=1, type='double')
;
;qtts_2_bin = [ey_peak, vxbz_peak, vzbx_peak]
;ey_types = ['efi', 'vxbz', 'vzbx']
;ey_ranges = [[-15., 20.], [0,0], [0,0]]
;ytitle = 'B!dqx!n/B!dlobe,q!n'
;xtitles = [pref_e+'E!dy,DSL'+"'"+'!n [mV/m]', '','']
;titles = ['Peak Values', '', '']
;abc = ['(a)','(b)','(c)']
;i_abc = 0
;
;limit = (ceil(1./binsize_h-0.5)+0.5)*binsize_h
;
;;;; more precise: bin boundaries
;bins_pos = [-1.,-0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 1]
;bins_neg = [-0.85, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.85, 1.]
;bin_boundaries = combine_arrays(bins_pos, bins_neg)
;
;;for i = 0, n_elements(ey_types)-1 do begin
;for i = 0, 0 do begin
;  if strcmp(i, 2) then begin
;	  xticknames = ''
;  endif else begin
;	  xticknames = [' ', ' ', ' ']
;  endelse
;
;  popen, pic_folder+'/ey_v_h'+dtrd_suf
;  print_options,xsize=4.5, ysize=3. ;; use this for single plot
;  thin_margin = 0.2 ;; character size
;  ;;;; write the letter
;;;; for old mean-median plot
;;	stat_plot, h, qtts_2_bin(i,*), k_c = k_c, bin_range = [-limit, limit], binsize = binsize_h, qtt_2_range = ey_ranges(*,i), qtt_range = [-1., 1.], qtt_2_title = xtitles(i), qtt_title = ytitle, kinbin = kinbin, bincntrs_out = bincenters, pm = 1, title = titles(i), bar_thick = l_thick, symsize = symsize, n_pst = n_pst, n_neg = n_neg, /no_write_pm, /vertical, /no_mean, color_med = 6, color_quar = 2, type_med = 'square'
;;;; quartiles
;	stat_plot, h, qtts_2_bin(i,*), k_c = k_c, bin_boundaries = bin_boundaries, qtt_2_range = ey_ranges(*,i), qtt_range = [-1., 1.], qtt_2_title = xtitles(i), qtt_title = ytitle, kinbin = kinbin, bincntrs_out = bincenters, pm = 1, title = titles(i), bar_thick = l_thick, symsize = symsize, n_pst = n_pst, n_neg = n_neg, /no_write_pm, /vertical, /no_mean, color_med = 6, color_quar = 2, type_med = 'square'
;	xyouts, !x.crange(0)+0.98*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.5*(!y.crange(1)-!y.crange(0)), 'Positive '+strcompress(string(n_pst), /remove), /data, alignment = 1, charsize = 0.6
;	xyouts, !x.crange(0)+0.02*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.5*(!y.crange(1)-!y.crange(0)), 'Negative '+strcompress(string(n_neg), /remove), /data, alignment = 0, charsize = 0.6
;	;xyouts, !x.crange(0)+0.32*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.9*(!y.crange(1)-!y.crange(0)), abc(i), /data, alignment = 0, charsize = 1.2 ;; for old median-bar plot
;	xyouts, !x.crange(0)+0.03*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.85*(!y.crange(1)-!y.crange(0)), abc(i), /data, alignment = 0, charsize = 1.2 ;; for new quartiles plot
;  ;;; write the right label
;  ;xyouts, 0.953, 0.63, method_s, orientation = -90., /normal, charsize = 1.4
;  pclose
;endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; the superpose of flux transport ;;;;;;;;;;;;;;;;;;;;;;;;;;
del_data, '*'

;;;; choose whether to follow full for vxbz and vzbx
vxb_suf = ''
;vxb_suf = '_followfull'

case dtrd_suf of
'': vbsuf = '_dtd'
'_vexb': vbsuf = ''
endcase

names = [$
'dbz_superpose_h', $
'intEy_dsl_superpose_h_positive'+dtrd_suf, $
'int_plasmaEy'+vbsuf+'_superpose_h_positive', $
'int_vxbz'+vbsuf+'_superpose_h_positive', $
'int_vzbx'+vbsuf+'_superpose_h_positive'$
]

;color_arr = [40, 62, 82, 97]
color_arr = [60, 132, 197, 250]

for i = 0, n_elements(names)-1 do begin
  name = names(i)
  if strcmp_or(strmid(name, 4, 4), ['vxbz', 'vzbx']) then name_file = name+vxb_suf else name_file = name
  datain, name, filename = save_folder+'/'+name_file+'.dat', dim = 4
	options, name, colors = color_arr
endfor

;;; create Bz vertical and horizontal lines
;;; the range for average
aft_range = [2., 3.]
;;; begin
get_data, 'dbz_superpose_h', data = dbz
time_clip, 'dbz_superpose_h', time_double('2000 1 1')+aft_range(0)*60., time_double('2000 1 1')+aft_range(1)*60.
get_data, 'dbz_superpose_h_tclip', data = after
horiz_lines = strarr(n_elements(dbz.y(0,*)))
vert_lines = strarr(n_elements(dbz.y(0,*)))
t_end = dbz.x(n_elements(dbz.x)-1)
for i = 0, n_elements(dbz.y(0,*))-1 do begin
	horiz_this = 'horiz_line_'+strcompress(string(i),/remove)
	vert_this = 'vert_line_'+strcompress(string(i),/remove)
	no_use = max(dbz.y(*,i), i_max, /nan)
	dbz_behind = dbz.y(i_max:*, i)
	t_behind = dbz.x(i_max:*)
	after_v = mean(after.y(*,i), /nan)
	i_neg = where(dbz_behind-after_v lt 0)
	i_first_neg = i_neg(0)
	i_last_pos = i_neg(0)-1
	t_cut = interpol(t_behind([i_last_pos, i_first_neg]), dbz_behind([i_last_pos, i_first_neg]), after_v)
	print, t_cut-time_double('2000 1 1')
	store_data, horiz_this, data = {x:[t_cut-0.06*(t_end-t_cut), t_end], y:[after_v, after_v]}
	store_data, vert_this, data = {x:[t_cut, t_cut], y:[-2, after_v+(after_v)*0.08]}
	options, horiz_this, thick = 0.7, color = color_arr(i)
	options, vert_this, thick = 0.7, color = color_arr(i), line = 1
	horiz_lines(i) = horiz_this
	vert_lines(i) = vert_this
endfor

;;; mark the durations to the flux plot
get_data, 'intEy_dsl_superpose_h_positive', data = flux 
store_data, 'intEy_dsl_superpose_h_positive', data = {x:flux.x, y:flux.y*6371000./1e9} ;; convert to MWb/RE
get_data, 'intEy_dsl_superpose_h_positive', data = flux 
horiz_lines_flux = strarr(n_elements(flux.y(0,*)))
vert_lines_flux = strarr(n_elements(flux.y(0,*)))
t_flux = flux.x
t_end = t_flux(n_elements(t_flux)-1)
for i = 0, n_elements(flux.y(0,*))-1 do begin
	horiz_fx_this = 'horiz_line_fx_'+strcompress(string(i),/remove)
	vert_fx_this = 'vert_line_fx_'+strcompress(string(i),/remove)
	flux_this = flux.y(*,i)
	get_data, vert_lines(i), data = vert
	t_dur = vert.x(0)
	flux_value = interpol(flux_this, t_flux, t_dur)
	print, 'flux of this DFB:'
	print, flux_value
	store_data, horiz_fx_this, data = {x:[t_dur-0.06*(t_end-t_cut), t_end], y: [flux_value, flux_value]}
	store_data, vert_fx_this, data = {x:[t_dur, t_dur], y: [-5, flux_value+flux_value*0.08]}
	options, horiz_fx_this, thick = 0.9, color = color_arr(i)
	options, vert_fx_this, thick = 0.9, color = color_arr(i)
	horiz_lines_flux(i) = horiz_fx_this
	vert_lines_flux(i) = vert_fx_this
endfor

;;; change units
get_data, 'int_plasmaEy'+vbsuf+'_superpose_h_positive',data = flux
store_data, 'int_plasmaEy'+vbsuf+'_superpose_h_positive',data = {x:flux.x, y:flux.y*6371000./1e9}
get_data, 'int_vxbz'+vbsuf+'_superpose_h_positive',data = flux
store_data, 'int_vxbz'+vbsuf+'_superpose_h_positive',data = {x:flux.x, y:flux.y*6371000./1e9}
get_data, 'int_vzbx'+vbsuf+'_superpose_h_positive',data = flux
store_data, 'int_vzbx'+vbsuf+'_superpose_h_positive',data = {x:flux.x, y:flux.y*6371000./1e9}

ylim, 'intEy_dsl_superpose_h_positive'+dtrd_suf, -20.*6371000./1e9, 400.*6371000./1e9
ylim, 'int_plasmaEy'+vbsuf+'_superpose_h_positive', -20.*6371000./1e9, 200.*6371000./1e9
ylim, 'int_vxbz'+vbsuf+'_superpose_h_positive', -8.*6371000./1e9, 180.*6371000./1e9
ylim, 'int_vzbx'+vbsuf+'_superpose_h_positive', -5.*6371000./1e9, 90.*6371000./1e9

;; generate the plots for MAIN_SUPERPOSE, run it first
options, 'int_plasmaEy'+vbsuf+'_superpose_h_positive'+dtrd_suf, ytitle=intg_sign+pref_e+'(V!dx!nB!dz!n-V!dz!nB!dx!n)!ndt';, ysubtitle = '[mWb/m]'
options, 'intEy_dsl_superpose_h_positive'+dtrd_suf, ytitle=intg_sign+pref_e+'E!dy,DSL'+"'"+'!ndt';, ysubtitle = '[mWb/m]'
options, 'int_vxbz'+vbsuf+'_superpose_h_positive', ytitle=intg_sign+pref_e+'(V!dx!nB!dz!n)dt';, ysubtitle = '[mWb/m]'
options, 'int_vzbx'+vbsuf+'_superpose_h_positive', ytitle=intg_sign+pref_e+'(-V!dz!nB!dx!n)dt';, ysubtitle = '[mWb/m]'

options, '*superpose_h*', xtickname=['-2','0','2'], thick=l_thick

;;; create the Bz horizontal lines plot
;store_data, 'dbz_superpose_h_lines', data=['dbz_superpose_h', horiz_lines, vert_lines]
;options, vert_lines, line = 1
;names(0) = 'dbz_superpose_h_lines'
;
;;; create the flux horizontal line plot
;store_data, 'intEy_dsl_superpose_h_positive_lines', data=['intEy_dsl_superpose_h_positive', horiz_lines_flux, vert_lines_flux]
;names(1) = 'intEy_dsl_superpose_h_positive_lines'
;ylim, 'intEy_dsl_superpose_h_positive_lines', -5, 269

;;; do not plot Bz and bars
names = names(1:*)

popen, pic_folder+'/superpose_h'+dtrd_suf+vxb_suf
print_options,xsize=5.,ysize=6
tplot_options,'vtitle',''
tplot, names, trange = ['1999 12 31 23 57', '2000 1 1 0 3'], title = 'Medians of Superposed Fluxes!d[MWb/R!IE!d]'
timebar, '2000 1 1 0 0', line = 1
timebar_mass, 0, varname=names, /databar, line=3
xyouts, 0.42,0.05, 'Minutes to t!d0!n', /normal
;;; label panels
xs=0.23*(1+dblarr(4,1))
ys=[0.28, 0.49, 0.7, 0.91]
ys_inc = -0.03
ss=['(e)', '(d)', '(c)','(b)']
xyouts, xs, ys, ss, charsize=1, /normal
xyouts, xs, ys(0:2)+ys_inc, 'Colors represent!cthe same regions!cas in panel b.', charsize = 0.5, /normal

;;; make the legend
charsize = 0.7
x_leg = 0.23
y_leg = 0.88
y_leg_inc = -0.025

xyouts, x_leg, y_leg,  '|B!dqx!n/B!dlobe,q!n|<0.2', charsize = charsize, color = color_arr(0), /normal
xyouts, x_leg, y_leg+y_leg_inc,  '0.2<|B!dqx!n/B!dlobe,q!n|<0.4', charsize = charsize, color = color_arr(1), /normal
xyouts, x_leg, y_leg+2*y_leg_inc,  '0.4<|B!dqx!n/B!dlobe,q!n|<0.6', charsize = charsize, color = color_arr(2), /normal
xyouts, x_leg, y_leg+3*y_leg_inc, '|B!dqx!n/B!dlobe,q!n|>0.6', charsize = charsize, color = color_arr(3), /normal

;;; make the label
ver=0.78
hor=0.02
pos_right=[0.82, 0.54]
make_px, pos_right, ver=ver, hor=hor, charname='Positive Only', lchar=0.21, orientation=-90, charsize = 1.05

pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end
