pro main_fluxvy_plot
;;;; plot dfb flux transport versus x, the paper's Figure 3.
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

;;;;;;;;;;;;;;;;;;;;;; the median/bar stat plots for peaks of things ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; panel specifics
xsize = 3
ysize = 5.3
left_margin = 0.15
right_margin = 0.01
top_margin = 0.05
bot_margin = 0.05
vspace = 0.007
n_panels = 3

;;; statistical specifics
dur_suf = '_40s'
binsize_x = 1.
k_c = 5
symsize = 0.8
xrange = [-12, 12]

;;; load data
events_fgs = load_list('dfb_list_lead_tail_fgs.txt', folder = listfolder)
events = load_list('dfb_list_lead_tail.txt', folder = listfolder)
t_fgs = time_double(events_fgs(0,*))
t = time_double(events(0,*))
y = datain_simple(save_folder+'/dfb_y.dat', dim=1, type='double')
y_fgs = datain_simple(save_folder+'/dfb_y_fgs.dat', dim=1, type='double')
Ey_peak = datain_simple(save_folder+'/dfb_ey_peak'+dtrd_suf+dur_suf+'_fgs.dat', dim=1, type='double')
dbz_peak = datain_simple(save_folder+'/dfb_dbz_peak'+dur_suf+'_fgs'+'.dat', dim=1, type='double')
df_i = datain_simple(save_folder+'/df_i'+'.dat', dim=1, type='double')
Vx_peak = datain_simple(save_folder+'/dfb_vx_peak'+dur_suf+'.dat', dim=1, type='double')
stop

store_data, 'ey', data=ey_peak
store_data, 'dbz', data=dbz_peak
store_data, 'i', data=df_i
store_data, 'vx', data=vx_peak

vars_2_bin = ['dbz', 'vx', 'ey']
;yranges = [[6,13], [-160, 400], [-15., 20.]] ;;; old: for median-mean plot of all events
yranges = [[6,13], [-190, 370], [-11.9, 14.]] ;;; new: for quartile plot of 0809 events
ytitles = [delta_letter+'B!dz!n [nT]', 'V!dx!n [km/s]', pref_e+'E!dy,DSL'+"'"+'!n [mV/m]']
titles = ['Peak Values', '','']
abc = ['(a)','(b)','(c)']

;;; bin boundaries
store_data, 'bin_dbz', data = [-12.,-6, -2.5, 0, 2.5, 5, 7.8, 12]
vx_bins_pos = [-12.,-6.5, -2.5, 0, 2.5, 4.6, 7.3, 12]
vx_bins_neg = [-12., -2.8, 0, 2.5, 5, 9, 12]
store_data, 'bin_vx', data = combine_arrays(vx_bins_pos, vx_bins_neg)
ey_bins_pos = [-12., -6, -2.7, 0, 2.5, 5, 7.8, 12]
ey_bins_neg = [-12., -3, 0, 2.5, 5, 7.8, 12]
store_data, 'bin_ey', data = combine_arrays(ey_bins_pos, ey_bins_neg)

;;; get the positions of the panels
positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [top_margin, bot_margin], space = [0., vspace], height = height)

popen, pic_folder+'/bin_vs_y'+dtrd_suf
print_options,xsize=xsize, ysize=ysize ;; use this for single plot
for i = 0, n_elements(vars_2_bin)-1 do begin
;for i = 0, 0 do begin
	if i eq n_elements(vars_2_bin)-1 then begin
	    xticknames = ''
		xtitle = 'Y!dGSM!n [R!dE!n]'
	endif else begin
	    xticknames = replicate(' ', 59)
		xtitle = ''
	endelse
	if strcmp(vars_2_bin(i), 'i') or strcmp(vars_2_bin(i), 'vx') then y_this = y else y_this = y_fgs 
	get_data, vars_2_bin(i), data = qtt_this
	get_data, 'bin_'+vars_2_bin(i), data = x_bin_boundaries

	stat_plot, transpose(y_this), transpose(qtt_this), k_c = k_c, bin_boundaries = x_bin_boundaries, qtt_2_range = yranges(*,i), qtt_range = xrange, qtt_2_title = ytitles(i), qtt_title = xtitle, qtt_tickname = xticknames, kinbin = kinbin, bincntrs_out = bincenters, pm = 1, title = titles(i), bar_thick = l_thick, symsize = symsize, n_pst = n_pst, n_neg = n_neg, /no_write_pm, position = positions(*,n_elements(vars_2_bin)-1-i), n_ttl = n_ttl, /noerase, /no_mean, color_med = 6, color_quar = 2, type_med = 'square'
	oplot, [0, 0], !y.crange, thick = 2
	if ~(strcmp(vars_2_bin(i),'dbz') or strcmp(vars_2_bin(i),'i')) then begin
		if strcmp(vars_2_bin(i), 'vx') then begin
			x_positive = 0.98
			y_positive = 0.4
			y_negative = 0.03
		endif
		if strcmp(vars_2_bin(i), 'ey') then begin
			x_positive = 0.26
			y_positive = 0.76
			y_negative = 0.1
		endif
		xyouts, !x.crange(0)+x_positive*(!x.crange(1)-!x.crange(0)), !y.crange(0)+y_positive*(!y.crange(1)-!y.crange(0)), 'Positive '+strcompress(string(n_pst), /remove), /data, alignment = 1, charsize = 0.6
		xyouts, !x.crange(0)+0.04*(!x.crange(1)-!x.crange(0)), !y.crange(0)+y_negative*(!y.crange(1)-!y.crange(0)), 'Negative '+strcompress(string(n_neg), /remove), /data, alignment = 0, charsize = 0.6
	endif else begin
		xyouts, !x.crange(0)+0.05*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.6*(!y.crange(1)-!y.crange(0)), strcompress(string(n_ttl), /remove), /data, alignment = 0., charsize = 0.6
	endelse
	xyouts, !x.crange(0)+0.04*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.86*(!y.crange(1)-!y.crange(0)), abc(i), /data, alignment = 0, charsize = 1.
endfor
pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; the superpose of flux transport, run main_superpose first ;;;;;;;;;;;;;;;;;;;;;;;;;;
;del_data, '*'
;
;;names = [$
;;'dbz_superpose_y', $
;;'vx_superpose_y_positive', $
;;'ey_dsl_superpose_'+dtrd_suf+'y_positive', $
;;'intEy_dsl_superpose'+dtrd_suf+'_y_positive'$
;;]
;
;names = [$
;'dbz_superpose_y', $
;'intEy_dsl_superpose'+dtrd_suf+'_y_positive'$
;]
;
;dim = 2
;;dim = 4
;
;if dim eq 2 then begin
;	color_arr = [2, 6]
;endif else begin
;	color_arr = [97, 82, 62, 40]
;endelse
;
;for i = 0, n_elements(names)-1 do begin
;  name = names(i)
;  datain, name, filename = save_folder+'/'+name+'.dat', dim = dim
;	options, name, colors = color_arr
;endfor
;
;;;; create Bz vertical and horizontal lines
;;;; the range for average
;aft_range = [2., 3.]
;;;; begin
;get_data, 'dbz_superpose_y', data = dbz
;time_clip, 'dbz_superpose_y', time_double('2000 1 1')+aft_range(0)*60., time_double('2000 1 1')+aft_range(1)*60.
;get_data, 'dbz_superpose_y_tclip', data = after
;horiz_lines = strarr(n_elements(dbz.y(0,*)))
;vert_lines = strarr(n_elements(dbz.y(0,*)))
;t_end = dbz.x(n_elements(dbz.x)-1)
;for i = 0, n_elements(dbz.y(0,*))-1 do begin
;	horiz_this = 'horiz_line_'+strcompress(string(i),/remove)
;	vert_this = 'vert_line_'+strcompress(string(i),/remove)
;	no_use = max(dbz.y(*,i), i_max, /nan)
;	dbz_behind = dbz.y(i_max:*, i)
;	t_behind = dbz.x(i_max:*)
;	after_v = mean(after.y(*,i), /nan)
;	i_neg = where(dbz_behind-after_v lt 0)
;	i_first_neg = i_neg(0)
;	i_last_pos = i_neg(0)-1
;	t_cut = interpol(t_behind([i_last_pos, i_first_neg]), dbz_behind([i_last_pos, i_first_neg]), after_v)
;	print, t_cut-time_double('2000 1 1')
;	store_data, horiz_this, data = {x:[t_cut-0.06*(t_end-t_cut), t_end], y:[after_v, after_v]}
;	store_data, vert_this, data = {x:[t_cut, t_cut], y:[-2, after_v+(after_v)*0.08]}
;	options, horiz_this, thick = 0.9, color = color_arr(i)
;	options, vert_this, thick = 0.9, color = color_arr(i), line = 1
;	horiz_lines(i) = horiz_this
;	vert_lines(i) = vert_this
;endfor
;
;;;; mark the durations to the flux plot
;get_data, 'intEy_dsl_superpose_y_positive', data = flux 
;store_data, 'intEy_dsl_superpose_y_positive', data = {x:flux.x, y:flux.y*6371000./1e9} ;; convert to MWb/RE
;get_data, 'intEy_dsl_superpose_y_positive', data = flux 
;horiz_lines_flux = strarr(n_elements(flux.y(0,*)))
;vert_lines_flux = strarr(n_elements(flux.y(0,*)))
;t_flux = flux.x
;t_end = t_flux(n_elements(t_flux)-1)
;for i = 0, n_elements(flux.y(0,*))-1 do begin
;	horiz_fx_this = 'horiz_line_fx_'+strcompress(string(i),/remove)
;	vert_fx_this = 'vert_line_fx_'+strcompress(string(i),/remove)
;	flux_this = flux.y(*,i)
;	get_data, vert_lines(i), data = vert
;	t_dur = vert.x(0)
;	flux_value = interpol(flux_this, t_flux, t_dur)
;	print, 'flux of this DFB:'
;	print, flux_value
;	store_data, horiz_fx_this, data = {x:[t_dur-0.06*(t_end-t_cut), t_end], y: [flux_value, flux_value]}
;	store_data, vert_fx_this, data = {x:[t_dur, t_dur], y: [-5, flux_value+flux_value*0.08]}
;	options, horiz_fx_this, thick = 0.9, color = color_arr(i)
;	options, vert_fx_this, thick = 0.9, color = color_arr(i)
;	horiz_lines_flux(i) = horiz_fx_this
;	vert_lines_flux(i) = vert_fx_this
;endfor
;
;;; generate the plots for MAIN_SUPERPOSE, run it first
;options, 'intEy_dsl_superpose_y_positive', ytitle=intg_sign+pref_e+'E!dy,DSL'+"'"+'!ndt [MWb/R!dE!n]'
;options, 'vx_superpose_y_positive', ytitle='V!dx!n [km/s]'
;options, 'dbz_superpose_y', ytitle=delta_letter+'B!dz!n [nT]'
;
;options, '*superpose_y*', xtickname=['-2','0','2'], thick=l_thick
;
;;; create the Bz horizontal lines plot
;store_data, 'dbz_superpose_y_lines', data=['dbz_superpose_y', horiz_lines, vert_lines]
;names(0) = 'dbz_superpose_y_lines'
;
;;; create the flux horizontal line plot
;store_data, 'intEy_dsl_superpose_y_positive_lines', data=['intEy_dsl_superpose_y_positive', horiz_lines_flux, vert_lines_flux]
;names(n_elements(names)-1) = 'intEy_dsl_superpose_y_positive_lines'
;ylim, 'intEy_dsl_superpose_y_positive_lines', -5*6371000./1e9, 209*6371000./1e9
;
;popen, pic_folder+'/superpose_y'+dtrd_suf
;;print_options,xsize=5,ysize=6 ;; for 3-panel plot
;print_options,xsize=6.5,ysize=3.7 ;; for 2-panel plot
;tplot, names, trange = ['1999 12 31 23 57', '2000 1 1 0 3'], title = 'Medians of Superposed Profiles'
;timebar, '2000 1 1 0 0', line = 1
;timebar_mass, 0, varname=names, /databar, line=3
;xyouts, 0.43,0.08, 'Minutes to t!d0!n', /normal
;
;;;; label panels
;;;; three-panel plot
;;xs=0.23*(1+dblarr(4,1))
;;ys=[0.35, 0.62,0.92]
;;ys_inc = -0.03
;;ss=['(f)', '(e)','(d)']
;;size_abc = 1.
;;;; two-panel plot
;xs=0.18*(1+dblarr(2,1))
;ys=[0.47, 0.84]
;ys_inc = -0.05
;ss=['(e)','(d)']
;size_abc = 1.4
;xyouts, xs, ys, ss, charsize=size_abc, /normal
;xyouts, xs, ys(0:n_elements(ys)-2)+ys_inc, 'Colors represent!cthe same regions!cas in panel d.', charsize = 0.5*size_abc, /normal
;
;;;; make the legend
;;;; three-panel plot
;;charsize = 0.7
;;x_leg = 0.23
;;y_leg = 0.88
;;y_leg_inc = -0.03
;;; two-panel plot
;charsize = 0.7*size_abc
;x_leg = 0.23
;y_leg = 0.82
;y_leg_inc = -0.05
;
;xyouts, x_leg, y_leg,  '-12R!dE!n<Y<0!n', charsize = charsize, color = color_arr(0), /normal
;xyouts, x_leg, y_leg+y_leg_inc, '0<Y<12R!dE!n', charsize = charsize, color = color_arr(1), /normal
;
;;ver=0.5
;;hor=0.02
;pos_right=[0.85, 0.39]
;;make_px, pos_right, ver=ver, hor=hor, charname='Positive Only', lchar=0.18, orientation=-90, charsize = 1
;;xyouts, pos_right(0), 0.85, 'All Events', charsize = 1.1, /normal, orientation=-90
;xyouts, pos_right(0), 0.51, 'Positive Only', charsize = 1.1, /normal, orientation=-90
;pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
stop
end
