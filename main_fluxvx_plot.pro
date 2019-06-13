pro main_fluxvx_plot
;;;; plot dfb flux transport versus x, the paper's Figure 3.
;;;; savefiles come from main_relations_dfb.pro and main_superpose.pro
thm_init
@declare

;;;;; choose which detrend method to use
dtrd_suf = ''
;dtrd_suf = '_vexb'

;season_suf = ''
season_suf = '_0809'

if strcmp(dtrd_suf, '_vexb') then begin
	pref_e = ''
endif else begin 
	pref_e = delta_letter
endelse

;;;;;;;;;;;;;;;;;;;;;; the median/bar stat plots for peaks of things ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; panel specifics
;xsize = 3
;ysize = 6.1
;left_margin = 0.15
;right_margin = 0.01
;top_margin = 0.05
;bot_margin = 0.05
;vspace = 0.007
;n_panels = 3
;
;;;; statistical specifics
;dur_suf = '_40s'
;binsize_x = 1.
;k_c = 5
;symsize = 0.8
;xrange = [-5, -24]
;
;;;; load data
;events_fgs = load_list('dfb_list_lead_tail_fgs.txt', folder = listfolder)
;events = load_list('dfb_list_lead_tail.txt', folder = listfolder)
;t_fgs = time_double(events_fgs(0,*))
;t = time_double(events(0,*))
;x = datain_simple(save_folder+'/dfb_x.dat', dim=1, type='double')
;x_fgs = datain_simple(save_folder+'/dfb_x_fgs.dat', dim=1, type='double')
;Ey_peak = datain_simple(save_folder+'/dfb_ey_peak'+dtrd_suf+dur_suf+'_fgs.dat', dim=1, type='double')
;dbz_peak = datain_simple(save_folder+'/dfb_dbz_peak'+dur_suf+'_fgs'+'.dat', dim=1, type='double')
;df_i = datain_simple(save_folder+'/df_i'+'.dat', dim=1, type='double')
;Vx_peak = datain_simple(save_folder+'/dfb_vx_peak'+dur_suf+'.dat', dim=1, type='double')
;
;;;;; for use of seting selection criterion
;print, median(dbz_peak, /even)
;
;if strcmp(season_suf, '_0809') then begin
;	i_0809 = where(t lt time_double('2010 1 1'))
;	x = x(*, i_0809)
;	Vx_peak = Vx_peak(*, i_0809)
;	df_i = df_i(*, i_0809)
;	i_0809_fgs = where(t_fgs lt time_double('2010 1 1'))
;	x_fgs = x_fgs(*, i_0809_fgs)
;	Ey_peak = Ey_peak(*, i_0809_fgs)
;	dbz_peak = dbz_peak(*, i_0809_fgs)
;endif
;
;store_data, 'ey', data=ey_peak
;store_data, 'dbz', data=dbz_peak
;store_data, 'i', data=df_i
;store_data, 'vx', data=vx_peak
;
;vars_2_bin = ['dbz', 'vx', 'ey']
;;yranges = [[6,13], [-160, 400], [-15., 20.]] ;;; old: for median-mean plot of all events
;yranges = [[6,13], [-190, 490], [-15., 20.]] ;;; new: for quartile plot of 0809 events
;ytitles = [delta_letter+'B!dz!n [nT]', 'V!dx!n [km/s]', pref_e+'E!dy,DSL'+"'"+'!n [mV/m]']
;titles = ['Peak Values', '','']
;abc = ['(a)','(b)','(c)']
;
;;;; new: for quartile plot of 0809 events
;store_data, 'bin_dbz', data = [-24, -17.8, -15.2,-13,-11.1,-9.5, -8,-6]
;vx_bins_pos = [-22.2, -19., -16.25,-14,-12,-10.3,-9.2, -7.5,-6]
;vx_bins_neg = [-14.5,-11.2,-9.9,-8.6, -7.3,-6.]
;store_data, 'bin_vx', data = combine_arrays(vx_bins_pos, vx_bins_neg)
;ey_bins_pos = [-22.5, -20,-16.5,-13.8,-11.3,-9.2, -7.5,-6]
;ey_bins_neg = [-20, -15, -11.5,-10.1,-9,-8,-7,-6]
;store_data, 'bin_ey', data = combine_arrays(ey_bins_pos, ey_bins_neg)
;
;;;; get the positions of the panels
;positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [top_margin, bot_margin], space = [0., vspace], height = height)
;
;;popen, pic_folder+'/bin_vs_x'+dtrd_suf+season_suf
;popen, 'bin_vs_x'+dtrd_suf+season_suf
;print_options,xsize=xsize, ysize=ysize ;; use this for single plot
;for i = 0, n_elements(vars_2_bin)-1 do begin
;;for i = 0, 0 do begin
;	if i eq n_elements(vars_2_bin)-1 then begin
;	    xticknames = ''
;		xtitle = 'X!dGSM!n [R!dE!n]'
;	endif else begin
;	    xticknames = replicate(' ', 59)
;		xtitle = ''
;	endelse
;	if strcmp(vars_2_bin(i), 'i') or strcmp(vars_2_bin(i), 'vx') then x_this = x else x_this = x_fgs 
;	get_data, vars_2_bin(i), data = qtt_this
;	get_data, 'bin_'+vars_2_bin(i), data = x_bin_boundaries
;
;	print, 'This quantity: '+vars_2_bin(i)
;	stat_plot, transpose(x_this), transpose(qtt_this), k_c = k_c, bin_boundaries = x_bin_boundaries, qtt_2_range = yranges(*,i), qtt_range = xrange, qtt_2_title = ytitles(i), qtt_title = xtitle, qtt_tickname = xticknames, kinbin = kinbin, bincntrs_out = bincenters, med = med, pm = 1, title = titles(i), bar_thick = l_thick, symsize = symsize, n_pst = n_pst, n_neg = n_neg, /no_write_pm, position = positions(*,n_elements(vars_2_bin)-1-i), n_ttl = n_ttl, /noerase, /no_mean, color_med = 6, color_quar = 2, type_med = 'square'
;	;oplot, [-20, -20], !y.crange, thick = 0.5
;	;print, bincenters
;	;print, med
;	if ~(strcmp(vars_2_bin(i),'dbz') or strcmp(vars_2_bin(i),'i')) then begin
;		case vars_2_bin(i) of
;		'ey': y_positive = 0.7
;		'vx': y_positive = 0.4
;		endcase
;		xyouts, !x.crange(0)+0.98*(!x.crange(1)-!x.crange(0)), !y.crange(0)+y_positive*(!y.crange(1)-!y.crange(0)), 'Positive '+strcompress(string(n_pst), /remove), /data, alignment = 1, charsize = 0.6
;		xyouts, !x.crange(0)+0.98*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.1*(!y.crange(1)-!y.crange(0)), 'Negative '+strcompress(string(n_neg), /remove), /data, alignment = 1, charsize = 0.6
;	endif else begin
;		xyouts, !x.crange(0)+0.9*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.76*(!y.crange(1)-!y.crange(0)), strcompress(string(n_ttl), /remove), /data, alignment = 0.5, charsize = 0.6
;	endelse
;	xyouts, !x.crange(0)+0.86*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.86*(!y.crange(1)-!y.crange(0)), abc(i), /data, alignment = 0, charsize = 1.
;endfor
;pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; the superpose of flux transport, run main_superpose first ;;;;;;;;;;;;;;;;;;;;;;;;;;
del_data, '*'

names = [$
'dbz_superpose_x', $
'vx_superpose_x_positive', $
;'ey_dsl_superpose_'+dtrd_suf+'x_positive', $
'intEy_dsl_superpose'+dtrd_suf+'_x_positive'$
]

;color_arr = [97, 82, 62, 40]
color_arr = [60, 132, 197, 250]

for i = 0, n_elements(names)-1 do begin
  name = names(i)
  datain, name, filename = save_folder+'/'+name+season_suf+'.dat', dim = 4
	options, name, colors = color_arr
endfor

;;; create Bz vertical and horizontal lines
;;; the range for average
aft_range = [2., 3.]
;;; begin
get_data, 'dbz_superpose_x', data = dbz
time_clip, 'dbz_superpose_x', time_double('2000 1 1')+aft_range(0)*60., time_double('2000 1 1')+aft_range(1)*60.
get_data, 'dbz_superpose_x_tclip', data = after
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
get_data, 'intEy_dsl_superpose_x_positive', data = flux 
store_data, 'intEy_dsl_superpose_x_positive', data = {x:flux.x, y:flux.y*6371000./1e9} ;; convert to MWb/RE
get_data, 'intEy_dsl_superpose_x_positive', data = flux 
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
options, 'intEy_dsl_superpose_x_positive', ytitle=intg_sign+pref_e+'E!dy,DSL'+"'"+'!ndt [MWb/R!dE!n]'
options, 'vx_superpose_x_positive', ytitle='V!dx!n [km/s]'
options, 'dbz_superpose_x', ytitle=delta_letter+'B!dz!n [nT]'

options, '*superpose_x*', xtickname=['-2','0','2'], thick=l_thick

;; create the Bz horizontal lines plot
store_data, 'dbz_superpose_x_lines', data=['dbz_superpose_x', horiz_lines, vert_lines]
options, vert_lines, line = 1
names(0) = 'dbz_superpose_x_lines'

;; create the flux horizontal line plot
store_data, 'intEy_dsl_superpose_x_positive_lines', data=['intEy_dsl_superpose_x_positive', horiz_lines_flux, vert_lines_flux]
names(n_elements(names)-1) = 'intEy_dsl_superpose_x_positive_lines'
ylim, 'intEy_dsl_superpose_x_positive_lines', -5*6371000./1e9, 269*6371000./1e9

;popen, pic_folder+'/superpose_x'+dtrd_suf+season_suf
popen, 'superpose_x'+dtrd_suf+season_suf
print_options,xsize=5,ysize=6
tplot_options,'vtitle',''
tplot, names, trange = ['1999 12 31 23 57', '2000 1 1 0 3'], title = 'Medians of Superposed Profiles'
timebar, '2000 1 1 0 0', line = 1
timebar_mass, 0, varname=names, /databar, line=3
xyouts, 0.42,0.05, 'Minutes to t!d0!n', /normal
;;; label panels
xs=0.23*(1+dblarr(4,1))
ys=[0.35, 0.62,0.92]
ys_inc = -0.03
ss=['(f)', '(e)','(d)']
xyouts, xs, ys, ss, charsize=1, /normal
xyouts, xs, ys(0:1)+ys_inc, 'Colors represent!cthe same regions!cas in panel d.', charsize = 0.5, /normal

;;; make the legend
charsize = 0.7
x_leg = 0.23
y_leg = 0.88
y_leg_inc = -0.03

xyouts, x_leg, y_leg,  '-20R!dE!n<X<-16R!dE!n', charsize = charsize, color = color_arr(0), /normal
xyouts, x_leg, y_leg+y_leg_inc, '-16R!dE!n<X<-12R!dE!n', charsize = charsize, color = color_arr(1), /normal
xyouts, x_leg, y_leg+2*y_leg_inc, '-12R!dE!n<X<-9R!dE!n', charsize = charsize, color = color_arr(2), /normal
xyouts, x_leg, y_leg+3*y_leg_inc, '-9R!dE!n<X<-6R!dE!n', charsize = charsize, color = color_arr(3), /normal

ver=0.5
hor=0.02
pos_right=[0.82, 0.39]
make_px, pos_right, ver=ver, hor=hor, charname='Positive Only', lchar=0.18, orientation=-90, charsize = 1
xyouts, pos_right(0), 0.88, 'All Events', charsize = 1, /normal, orientation=-90
pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;; the median/bar stat plots for average and total transport ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; panel specifics
xsize = 3.2
ysize = 6.9
left_margin = 0.15
right_margin = 0.01
top_margin = 0.05
top_panel = 0.16
bot_margin = 0.05
vspace = 0.007
n_panels = 3

;;; statistical specifics
dur_suf = '_40s'
binsize_x = 1.
k_c = 5
symsize = 0.8
xrange = [-5, -24]

;;; load data
events_fgs = load_list('dfb_list_lead_tail_fgs.txt', folder = listfolder)
t_fgs = time_double(events_fgs(0,*))
x_fgs = datain_simple(save_folder+'/dfb_x_fgs.dat', dim=1, type='double')
ey_ave_fgs = datain_simple(save_folder+'/dfb_ey_ave'+dtrd_suf+'_xdur_fgs.dat', dim=1, type='double')

events = load_list('dfb_list_lead_tail.txt', folder = listfolder)
t = time_double(events(0,*))
x = datain_simple(save_folder+'/dfb_x.dat', dim=1, type='double')
length = datain_simple(save_folder+'/dfb_length_xdur.dat', dim=1, type='double')

;;; load DFB duration
sec_dfb = datain_simple(save_folder+'/dfb_duration_fgs.dat', dim=1, type='double')

if strcmp(season_suf, '_0809') then begin
	i_0809_fgs = where(t_fgs lt time_double('2010 1 1'))
	x_fgs = x_fgs(*, i_0809_fgs)
	ey_ave_fgs = ey_ave_fgs(*, i_0809_fgs)
	sec_dfb = sec_dfb(*, i_0809_fgs)
	i_0809 = where(t lt time_double('2010 1 1'))
	x = x(*, i_0809)
	length = length(*, i_0809)
endif

flux = ey_ave_fgs*sec_dfb*6371000./1e9 ;; in MWb/RE
length = length/6371.

store_data, 'length', data=length
store_data, 'ey', data=ey_ave_fgs
store_data, 'flux', data=flux

vars_2_bin = ['length', 'ey', 'flux']
yranges = [[-0.8, 4.5], [-3.5, 5.99], [-99*6371000./1e9, 220*6371000./1e9]]
ytitles = ['DFB X-length [R!dE!n]', jiao_l+pref_e+'E!dy,DSL'+"'"+'!n'+jiao_r+' [mV/m]', phi_letter+' [MWb/R!dE]']
titles = ['', '', '']
abc = ['(b)','(c)','(d)']

;;;; for quartiles plot
length_bins_pos = [-22.2, -18.4, -16.3,-14.1,-12,-10.3,-9.2, -7.5,-6]
length_bins_neg = [-14.5,-11.2,-9.9,-8.9, -7.8,-6.]
store_data, 'bin_length', data = combine_arrays(length_bins_pos, length_bins_neg)
ey_bins_pos = [-22.5, -20.02,-16.8,-14.25,-11.3,-9.6, -7.8,-6]
ey_bins_neg = [-20.1, -12, -10.5, -9.2, -8., -6]
store_data, 'bin_ey', data = combine_arrays(ey_bins_pos, ey_bins_neg)
flux_bins_pos = [-23,-15,-11.2, -9.7, -8.,-6]
flux_bins_neg = [-20, -12,-10,-8.3,-6]
store_data, 'bin_flux', data = combine_arrays(flux_bins_pos, flux_bins_neg)

;;; get the positions of the panels
positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [top_margin, bot_margin+top_panel+vspace], space = [0., vspace], height = height)

;popen, pic_folder+'/bin_ave_vs_x'+dtrd_suf+season_suf
popen, 'bin_ave_vs_x'+dtrd_suf+season_suf
print_options,xsize=xsize, ysize=ysize ;; use this for single plot

;; first plot the durations
plot, x_fgs, sec_dfb, xtitle = '', ytitle = 'DFB duration [s]', xtickname = replicate(' ', 59), xrange = xrange, xstyle = 1, position = [left_margin, 1.-(top_panel+top_margin), 1.-right_margin, 1.-top_margin], yrange = [0, 199], ystyle = 1, psym = 3;, symsize = 0.2
x_marks = [-17.722568, -14.013473, -10.338759, -7.9557382] ;;; in RE for all events
dur_marks = [97.32, 47.81, 34.555398, 27.695121] ;;; in seconds
oplot, x_marks, dur_marks, psym = 6, symsize = 0.7
y_positive = 0.65
xyouts, !x.crange(0)+0.15*(!x.crange(1)-!x.crange(0)), !y.crange(0)+y_positive*(!y.crange(1)-!y.crange(0)), strcompress(string(n_elements(where(finite(x_fgs) and finite(sec_dfb)))), /remove), /data, alignment = 0, charsize = 0.6
xyouts, !x.crange(0)+0.05*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.8*(!y.crange(1)-!y.crange(0)), '(a)', /data, alignment = 0, charsize = 1.

for i = 0, n_elements(vars_2_bin)-1 do begin
	if i eq n_elements(vars_2_bin)-1 then begin
	    xticknames = ''
		xtitle = 'X!dGSM!n [R!dE!n]'
	endif else begin
	    xticknames = replicate(' ', 59)
		xtitle = ''
	endelse
	if i eq 0 then x_this = x else x_this = x_fgs 
	get_data, vars_2_bin(i), data = qtt_this
	get_data, 'bin_'+vars_2_bin(i), data = x_bin_boundaries

	print, 'This quantity: '+vars_2_bin(i)
	stat_plot, transpose(x_this), transpose(qtt_this), k_c = k_c, bin_boundaries = x_bin_boundaries, qtt_2_range = yranges(*,i), qtt_range = xrange, qtt_2_title = ytitles(i), qtt_title = xtitle, qtt_tickname = xticknames, kinbin = kinbin, bincntrs_out = bincenters, pm = 1, title = titles(i), bar_thick = l_thick, symsize = symsize, n_pst = n_pst, n_neg = n_neg, /no_write_pm, position = positions(*,n_elements(vars_2_bin)-1-i), n_ttl = n_ttl, /noerase, /no_mean, color_med = 6, color_quar = 2, type_med = 'square'
	;oplot, [-20, -20], !y.crange, thick = 0.5
	;;; write event numbers
	case vars_2_bin(i) of
	'length': y_positive = 0.7
	'ey': y_positive = 0.7
	'flux': y_positive = 0.4
	endcase
	case vars_2_bin(i) of
	'length': y_negative = 0.05
	'ey': y_negative= 0.1
	'flux': y_negative = 0.1
	endcase
	case vars_2_bin(i) of
	'length': x_positive = 0.25
	'ey': x_positive= 0.98
	'flux': x_positive = 0.98
	endcase
	case vars_2_bin(i) of
	'length': x_abc = 0.05
	'ey': x_abc= 0.86
	'flux': x_abc = 0.86
	endcase
	xyouts, !x.crange(0)+x_positive*(!x.crange(1)-!x.crange(0)), !y.crange(0)+y_positive*(!y.crange(1)-!y.crange(0)), 'Positive '+strcompress(string(n_pst), /remove), /data, alignment = 1, charsize = 0.6
	xyouts, !x.crange(0)+0.98*(!x.crange(1)-!x.crange(0)), !y.crange(0)+y_negative*(!y.crange(1)-!y.crange(0)), 'Negative '+strcompress(string(n_neg), /remove), /data, alignment = 1, charsize = 0.6
	xyouts, !x.crange(0)+x_abc*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.86*(!y.crange(1)-!y.crange(0)), abc(i), /data, alignment = 0, charsize = 1.
endfor
pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
stop
end
