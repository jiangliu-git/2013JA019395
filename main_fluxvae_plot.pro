pro main_fluxvae_plot
;;; plot flux distributions versus AE
;;;; savefiles come from main_relations_dfb.pro and main_superpose.pro
thm_init
@declare

;;;;; choose which detrend method to use
dtrd_suf = ''
;dtrd_suf = '_vexb'

;;;; choose whether to use 0-14UT only 
;ut_suf = ''
ut_suf = '_014' ;; 0-14 UT only

if strcmp(dtrd_suf, '_vexb') then begin
	pref_e = ''
endif else begin 
	pref_e = delta_letter
endelse

suffix = '_fgs'
;;;;;;;;;;;;;;;;;;;;;;;; the median/bar stat plots for ey peaks ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;dur_suf = '_40s'
;binsize_h = 0.23
;k_c = 5
;symsize = 0.8
;;;; load data
;ae_ave = datain_simple(save_folder+'/pseudo_ae_ave'+suffix+'.dat', dim=1, type='double')
;ey_peak = datain_simple(save_folder+'/dfb_ey_peak'+dur_suf+suffix+dtrd_suf+'.dat', dim=1, type='double')
;dbz_peak = datain_simple(save_folder+'/dfb_dbz_peak'+dur_suf+'_fgs'+'.dat', dim=1, type='double')
;;;; load list
;listname = 'dfb_list_lead_tail'+suffix+'.txt'
;events = load_list(listname, folder = listfolder)
;time = time_double(events(0,*))
;date_str = strmid(events(0,*), 0, 10)
;date = time_double(date_str)
;ut=(time-date)/3600.
;;;; get 0-14 UT events ;;;
;if strcmp(ut_suf, '_014') then begin
;	i_good = where(ut le 14.)
;	events = events(*, i_good)
;	ae_ave = ae_ave(*, i_good)
;	ey_peak = ey_peak(*, i_good)
;	dbz_peak = dbz_peak(*, i_good)
;endif
;
;store_data, 'ey', data=ey_peak
;store_data, 'dbz', data=dbz_peak
;
;vars_2_bin = ['dbz', 'ey']
;;yranges = [[7,14], [-17., 18.]] ;; for median-mean-std plot
;yranges = [[6.1, 15], [-13., 18.]] ;; for quartiles plot
;xrange = [0., 1199.] ;; for quartiles plot
;pm = [0, 1]
;ytitles = [delta_letter+'B!dz!n [nT]', pref_e+'E!dy,DSL'+"'"+'!n [mV/m]']
;titles = ['Peak Values', '']
;abc = ['(a)','(b)']
;;title = 'All events'
;title = ''
;
;;;; bin boudnarys of AE
;;;;;; for median-mean-std plot
;;store_data, 'bin_ey', data = [[0., 70., 140, 210, 280, 360, 470, 630, 800, 1200], [0, 70, 140, 220, 360, 500, 650, 1200, !values.f_nan, !values.f_nan]] ;; positive and negative
;;store_data, 'bin_dbz', data = [0., 70., 140, 210, 290, 380, 480, 630, 780, 980, 1200]
;
;;;;;; for quartile plot, all events
;;bin_ey_pos = [0., 70., 140, 260, 400, 800, 1200]
;;bin_ey_neg = [0., 70., 140, 260, 500, 1200]
;;store_data, 'bin_ey', data = combine_arrays(bin_ey_pos, bin_ey_neg) ;; positive and negative
;;store_data, 'bin_dbz', data = [0., 70., 140, 250, 370, 540, 720, 940, 1200]
;
;;;;;; for quartile plot, 0-14UT only events
;bin_ey_pos = [0., 70., 140, 260, 400, 800, 1200]
;bin_ey_neg = [0., 70., 140, 250, 500, 1200]
;store_data, 'bin_ey', data = combine_arrays(bin_ey_pos, bin_ey_neg) ;; positive and negative
;store_data, 'bin_dbz', data = [0., 70., 140, 250, 390, 540, 720, 940, 1200]
;
;positions = panel_positions([1, 2], space = [0., 0.01])
;
;popen, pic_folder+'/bin_vs_ae'+dtrd_suf+ut_suf
;print_options,xsize=3.6, ysize=5. ;; use this for single plot
;for i = 0, n_elements(vars_2_bin)-1 do begin
;	if i eq n_elements(vars_2_bin)-1 then begin
;	    xticknames = ''
;		xtitle = 'THEMIS Pseudo AE [nT]'
;	endif else begin
;	    xticknames = replicate(' ', 59)
;		xtitle = ''
;	endelse
;	get_data, vars_2_bin(i), data = qtt_this
;	get_data, 'bin_'+vars_2_bin(i), data = ae_bin_boundaries
;
;	stat_plot, transpose(ae_ave), transpose(qtt_this), k_c = k_c, bin_boundaries = ae_bin_boundaries, qtt_2_range = yranges(*,i), qtt_range = xrange, qtt_2_title = ytitles(i), qtt_title = xtitle, qtt_tickname = xticknames, kinbin = kinbin, bincntrs_out = bincenters, pm = pm(i), title = titles(i), bar_thick = l_thick, symsize = symsize, n_pst = n_pst, n_neg = n_neg, /no_write_pm, position = positions(*,n_elements(vars_2_bin)-1-i), n_ttl = n_ttl, /noerase, /no_mean, color_med = 6, color_quar = 2, type_med = 'square'
;	;;;;; write words
;	;;; old
;	;if ~(strcmp(vars_2_bin(i),'dbz')) then begin
;	;	xyouts, !x.crange(0)+0.98*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.7*(!y.crange(1)-!y.crange(0)), 'Positive '+strcompress(string(n_pst), /remove), /data, alignment = 1, charsize = 0.6
;	;	xyouts, !x.crange(0)+0.98*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.1*(!y.crange(1)-!y.crange(0)), 'Negative '+strcompress(string(n_neg), /remove), /data, alignment = 1, charsize = 0.6
;	;	xyouts, !x.crange(0)+0.86*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.86*(!y.crange(1)-!y.crange(0)), abc(i), /data, alignment = 0, charsize = 1.2
;	;endif else begin
;	;	xyouts, !x.crange(0)+0.9*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.22*(!y.crange(1)-!y.crange(0)), strcompress(string(n_ttl), /remove), /data, alignment = 0.5, charsize = 0.6
;	;	xyouts, !x.crange(0)+0.86*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.12*(!y.crange(1)-!y.crange(0)), abc(i), /data, alignment = 0, charsize = 1.2
;	;endelse
;	;;; new
;	if ~(strcmp(vars_2_bin(i),'dbz')) then begin
;		xyouts, !x.crange(0)+0.98*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.63*(!y.crange(1)-!y.crange(0)), 'Positive '+strcompress(string(n_pst), /remove), /data, alignment = 1, charsize = 0.6
;		xyouts, !x.crange(0)+0.98*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.1*(!y.crange(1)-!y.crange(0)), 'Negative '+strcompress(string(n_neg), /remove), /data, alignment = 1, charsize = 0.6
;		xyouts, !x.crange(0)+0.86*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.2*(!y.crange(1)-!y.crange(0)), abc(i), /data, alignment = 0, charsize = 1.2
;	endif else begin
;		xyouts, !x.crange(0)+0.9*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.22*(!y.crange(1)-!y.crange(0)), strcompress(string(n_ttl), /remove), /data, alignment = 0.5, charsize = 0.6
;		xyouts, !x.crange(0)+0.86*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.12*(!y.crange(1)-!y.crange(0)), abc(i), /data, alignment = 0, charsize = 1.2
;	endelse
;endfor
;pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; the superpose of flux transport ;;;;;;;;;;;;;;;;;;;;;;;;;;
del_data, '*'

names = [$
'dbz_superpose_ae',$
'intEy_dsl_superpose'+dtrd_suf+'_ae_positive'$
]

;color_arr = [40, 62, 82, 97]
;color_arr = [100, 82, 62, 40]
color_arr = [60, 132, 197, 250]

for i = 0, n_elements(names)-1 do begin
  name = names(i)
  datain, name, filename = save_folder+'/'+name+ut_suf+'.dat', dim = 4
	options, name, colors = color_arr
endfor

;;; create Bz vertical and horizontal lines
;;; the range for average
aft_range = [2., 3.]
;;; begin
get_data, 'dbz_superpose_ae', data = dbz
time_clip, 'dbz_superpose_ae', time_double('2000 1 1')+aft_range(0)*60., time_double('2000 1 1')+aft_range(1)*60.
get_data, 'dbz_superpose_ae_tclip', data = after
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
	options, horiz_this, thick = 0.9, color = color_arr(i)
	options, vert_this, thick = 0.9, color = color_arr(i), line = 1
	horiz_lines(i) = horiz_this
	vert_lines(i) = vert_this
endfor

;;; mark the durations to the flux plot
get_data, 'intEy_dsl_superpose_ae_positive', data = flux 
store_data, 'intEy_dsl_superpose_ae_positive', data = {x:flux.x, y:flux.y*6371000./1e9} ;; convert to MWb/RE
get_data, 'intEy_dsl_superpose_ae_positive', data = flux 
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

;; generate the plots for MAIN_SUPERPOSE, run it first
options, 'intEy_dsl_superpose_ae_positive', ytitle=intg_sign+pref_e+'E!dy,DSL'+"'"+'!ndt [MWb/R!dE!n]'
options, 'dbz_superpose_ae', ytitle=delta_letter+'B!dz!n [nT]'

options, '*superpose_ae*', xtickname=['-2','0','2'], thick=l_thick

;; create the Bz horizontal lines plot
store_data, 'dbz_superpose_ae_lines', data=['dbz_superpose_ae', horiz_lines, vert_lines]
names(0) = 'dbz_superpose_ae_lines'

;; create the flux horizontal line plot
store_data, 'intEy_dsl_superpose_ae_positive_lines', data=['intEy_dsl_superpose_ae_positive', horiz_lines_flux, vert_lines_flux]
names(n_elements(names)-1) = 'intEy_dsl_superpose_ae_positive_lines'
ylim, 'intEy_dsl_superpose_ae_positive_lines', -5*6371000./1e9, 310*6371000./1e9

popen, pic_folder+'/superpose_ae'+dtrd_suf+ut_suf
print_options,xsize=5,ysize=4
tplot_options,'vtitle',''
tplot, names, trange = ['1999 12 31 23 57', '2000 1 1 0 3'], title = 'Medians of Superposed Profiles'
timebar, '2000 1 1 0 0', line = 1
timebar_mass, 0, varname=names, /databar, line=3
xyouts, 0.42,0.08, 'Minutes to t!d0!n', /normal
;;; label panels
xs=0.23*(1+dblarr(2))
ys=[0.87, 0.48]
ys_inc = -0.04
ss=['(c)', '(d)']
xyouts, xs, ys, ss, charsize=1, /normal
xyouts, xs, ys(1)+ys_inc, 'Colors represent!cthe same regions!cas in panel c.', charsize = 0.5, /normal

;;; manually add the title
;xyouts, 0.8, 0.73, 'All Events', /normal, align = 0.5, charsize = 1.0, orientation = -90
xyouts, 0.8, 0.35, 'Positive Only', /normal, align = 0.5, charsize = 1.0, orientation = -90

;;; make the legend
charsize = 0.5
x_leg = 0.225
y_leg = 0.83
y_leg_inc = -0.03

if strcmp(ut_suf, '_014') then begin
	;;;; 0-14UT events ;;;
	xyouts, x_leg, y_leg,  'AE>400nT (X!dmed!n=-9.41R!dE!n)', charsize = charsize, color = color_arr(3), /normal
	xyouts, x_leg, y_leg+y_leg_inc,  '200nT<AE<400nT (X!dmed!n=-9.42R!dE!n)', charsize = charsize, color = color_arr(2), /normal
	xyouts, x_leg, y_leg+2*y_leg_inc,  '100nT<AE<200nT (X!dmed!n=-9.59R!dE!n)', charsize = charsize, color = color_arr(1), /normal
	xyouts, x_leg, y_leg+3*y_leg_inc,  'AE<100nT (X!dmed!n=-10.06R!dE!n)', charsize = charsize, color = color_arr(0), /normal
endif else begin
	;;;; all events ;;;
	xyouts, x_leg, y_leg,  'AE>400nT (X!dmed!n=-9.30R!dE!n)', charsize = charsize, color = color_arr(3), /normal
	xyouts, x_leg, y_leg+y_leg_inc,  '200nT<AE<400nT (X!dmed!n=-9.24R!dE!n)', charsize = charsize, color = color_arr(2), /normal
	xyouts, x_leg, y_leg+2*y_leg_inc,  '100nT<AE<200nT (X!dmed!n=-9.51R!dE!n)', charsize = charsize, color = color_arr(1), /normal
	xyouts, x_leg, y_leg+3*y_leg_inc,  'AE<100nT (X!dmed!n=-9.97R!dE!n)', charsize = charsize, color = color_arr(0), /normal
endelse

;;; 200-350 X_med: 9.22
;;; 350+ X_med: 9.36

pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
stop
end
