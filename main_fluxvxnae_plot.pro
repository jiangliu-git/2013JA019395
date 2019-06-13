pro main_fluxvxnae_plot
;;;; plot dfb flux transport versus x in separated AE conditions, the paper's Figure 3.
;;;; savefiles come from main_relations_dfb.pro and main_superpose.pro
thm_init
@declare

;;;;; choose which detrend method to use
dtrd_suf = ''
;dtrd_suf = '_vexb'

;season_suf = ''
season_suf = '_0809'
;season_suf = '_both' ;; plot 0809 as solid and 1011 as dash

;;;; choose whether to use 0-14UT only 
;ut_suf = ''
ut_suf = '_014' ;; 0-14 UT only

if strcmp(dtrd_suf, '_vexb') then begin
	pref_e = ''
endif else begin 
	pref_e = delta_letter
endelse

ae_ranges = 'fix' ;; fixed ae values
;ae_ranges = 'divide' ;; determine from number of events to be even

;;; ranges of AE indices (fixed value)
ae_4 = [0., 100, 200, 400, 1000000000.] ;; four levels
color_ae_4 = [100, 82, 62, 40]
ae_3 = [0., 150, 300, 1000000000.] ;; two levels
;color_ae_3 = [92, 73, 42]
color_ae_3 = [78, 148, 230]
;color_ae_3 = [6, 4, 2]
ae_2 = [0., 200, 1000000000.] ;; two levels
color_ae_2 = [82, 40]
ae_1 = [0., 1000000000.] ;; one level only
color_ae_1 = 0

;;; whether use only positive/negative transport events, told by ey_dsl
;sign_use_suf = ''
sign_use_suf = '_positive'
;sign_use_suf = '_negative'

;;;;;;;;;;;;;;;;;;;;;; the median/bar stat plots for peaks of things ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; panel specifics
;;;;; four panel plot
;xsize = 3
;ysize = 7.6
;n_panels = 4
;;;;; two panel plot
xsize = 3
ysize = 4.2
n_panels = 2
;;;;;; one panel plot
;xsize = 3
;ysize = 3
;n_panels = 1

left_margin = 0.15
right_margin = 0.01
top_margin = 0.05
bot_margin = 0.05
vspace = 0.007

;;; statistical specifics
binsize_x = 1.
k_c = 5
symsize = 0.8
xrange = [-5, -29.99]

;;; load value data
events_fgs = load_list('dfb_list_lead_tail_fgs.txt', folder = listfolder)
t_fgs = time_double(events_fgs(0,*))
x_fgs = datain_simple(save_folder+'/dfb_x_fgs.dat', dim=1, type='double')
ey_peak = datain_simple(save_folder+'/dfb_ey_peak'+dtrd_suf+'_40s_fgs.dat', dim=1, type='double')
ey_ave = datain_simple(save_folder+'/dfb_ey_ave'+dtrd_suf+'_xdur_fgs.dat', dim=1, type='double')
ae_ave = datain_simple(save_folder+'/pseudo_ae_ave_fgs.dat', dim=1, type='double')
;;; load DFB duration
sec_dfb = datain_simple(save_folder+'/dfb_duration_fgs.dat', dim=1, type='double')
;;; load occurance data
orbit_pos_ae = datain_simple(save_folder+'/orbit_xyz_ae_fgs.dat', dim=5, type='double')
t_orbit = orbit_pos_ae(0,*)

;;; get 0-14 UT DFB events ;;;
date_str = strmid(events_fgs(0,*), 0, 10)
date = time_double(date_str)
ut=(t_fgs-date)/3600.
if strcmp(ut_suf, '_014') then begin
	i_good = where(ut le 14.)
	events_fgs = events_fgs(*, i_good)
	t_fgs = t_fgs(*, i_good)
	x_fgs = x_fgs(*, i_good)
	ey_peak =  ey_peak(*, i_good)
	ey_ave = ey_ave(*, i_good)
	ae_ave = ae_ave(*, i_good)
	sec_dfb = sec_dfb(*, i_good)
endif

;;; get 0-14 UT orbits
t_orbit_str = time_string(t_orbit)
date_orbit_str = strmid(t_orbit_str, 0, 10)
date_orbit = time_double(date_orbit_str)
ut_orbit=(t_orbit-date_orbit)/3600.
if strcmp(ut_suf, '_014') then begin
	i_good = where(ut_orbit le 14.)
	orbit_pos_ae = orbit_pos_ae(*, i_good)
	t_orbit = t_orbit(*, i_good)
endif

if strcmp_or(season_suf, ['_0809', '_both']) then begin
	i_0809_fgs = where(t_fgs lt time_double('2010 1 1'))
	i_1011_fgs = where(t_fgs gt time_double('2010 1 1'))
	x_fgs_1011 = x_fgs(*, i_1011_fgs)
	x_fgs = x_fgs(*, i_0809_fgs)
	ey_peak_1011 = ey_peak(*, i_1011_fgs)
	ey_peak = ey_peak(*, i_0809_fgs)
	ey_ave_1011 = ey_ave(*, i_1011_fgs)
	ey_ave = ey_ave(*, i_0809_fgs)
	sec_dfb_1011 = sec_dfb(*, i_1011_fgs)
	sec_dfb = sec_dfb(*, i_0809_fgs)
	ae_ave_1011 = ae_ave(*, i_1011_fgs)
	ae_ave = ae_ave(*, i_0809_fgs)
	i_0809_orbit = where(t_orbit lt time_double('2010 1 1'))
	i_1011_orbit = where(t_orbit gt time_double('2010 1 1'))
	orbit_pos_ae_1011 = orbit_pos_ae(*, i_1011_orbit)
	orbit_pos_ae = orbit_pos_ae(*, i_0809_orbit)
endif

if strcmp(season_suf, '_both') then begin
	pos_orbit_1011 = orbit_pos_ae_1011(1:3,*)
	ae_orbit_1011 = orbit_pos_ae_1011(4,*)
	;;; for transport
	ey_ave_1011 = ey_ave_1011*sec_dfb_1011
endif
pos_orbit = orbit_pos_ae(1:3,*)
ae_orbit = orbit_pos_ae(4,*)
;;; for transport
ey_ave = ey_ave*sec_dfb


;;; split events
case sign_use_suf of
'': i_use_peak = intarr(n_elements(ey_peak))
'_positive': i_use_peak = where(ey_peak gt 0)
'_negative': i_use_peak = where(ey_peak lt 0)
endcase
case sign_use_suf of
'': i_use_ave = intarr(n_elements(ey_ave))
'_positive': i_use_ave = where(ey_ave gt 0)
'_negative': i_use_ave = where(ey_ave lt 0)
endcase
x_peak = x_fgs(i_use_peak)
x_ave = x_fgs(i_use_ave)
ey_peak = ey_peak(i_use_peak)
ey_ave = ey_ave(i_use_ave)
ae_ave_peak = ae_ave(i_use_peak)
ae_ave_ave = ae_ave(i_use_ave)
sec_dfb = sec_dfb(i_use_peak)

store_data, 'ey_peak', data={x:x_peak, qtt:ey_peak, ae:ae_ave_peak}
store_data, 'ey_ave', data={x:x_ave, qtt:ey_ave, ae:ae_ave_ave}

if strcmp(season_suf, '_both') then begin
	case sign_use_suf of
	'': i_use_peak = intarr(n_elements(ey_peak_1011))
	'_positive': i_use_peak = where(ey_peak_1011 gt 0)
	'_negative': i_use_peak = where(ey_peak_1011 lt 0)
	endcase
	case sign_use_suf of
	'': i_use_ave = intarr(n_elements(ey_ave_1011))
	'_positive': i_use_ave = where(ey_ave_1011 gt 0)
	'_negative': i_use_ave = where(ey_ave_1011 lt 0)
	endcase
	x_peak_1011 = x_fgs_1011(i_use_peak)
	x_ave_1011 = x_fgs_1011(i_use_ave)
	ey_peak_1011 = ey_peak_1011(i_use_peak)
	ey_ave_1011 = ey_ave_1011(i_use_ave)
	ae_ave_peak_1011 = ae_ave_1011(i_use_peak)
	ae_ave_ave_1011 = ae_ave_1011(i_use_ave)
	sec_dfb_1011 = sec_dfb_1011(i_use_peak)

	store_data, 'ey_peak', data={x:x_peak_1011, qtt:ey_peak_1011, ae:ae_ave_peak_1011}
	store_data, 'ey_ave', data={x:x_ave_1011, qtt:ey_ave_1011, ae:ae_ave_ave_1011}
endif


;;; get the positions of the panels
positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [top_margin, bot_margin], space = [0., vspace], height = height)

case sign_use_suf of
'': sign_title = 'All Events'
'_positive': sign_title = 'Positive '+delta_letter+'E!S!dy,DSL'+"'"+'!R!Upeak!n DFBs'
'_negative': sign_title = 'Negative Events'
endcase

popen, pic_folder+'/occurrence'+dtrd_suf+season_suf+ut_suf
print_options,xsize=xsize, ysize=ysize ;; use this for single plot

;;;;;; plot the two bin plots
;;; settings
;;;; choose number of AE levels
;ae_levels = ae_1
;color_ae = color_ae_1
;vars_2_bin = ['ey_peak', 'ey_ave']
;;yranges = [[0,25], [0, 10]]
;yranges = [[0,25], [0, 400]] ;; for flux transport
;;ytitles = [pref_e+'E!S!dy,DSL'+"'"+'!R!Upeak!n [mV/m]', pref_e+'E!S!dy,DSL'+"'"+'!R!Uaverage!n [mV/m]']
;ytitles = [pref_e+'E!S!dy,DSL'+"'"+'!R!Upeak!n [mV/m]', 'Transport [mWb/m]']
;titles = [sign_title, '']
;abc = ['(a)','(b)']
;xticknames = replicate(' ', 59)
;xtitle = ''
;;;; set x bin boundaries
;store_data, 'bin_x_ey_peak', data = [$
;[-20,-16.5,-13.8,-11.3,-9.2, -7.5,-6], $
;[-20,-16.5,-13.8,-11.3,-9.2, -7.5,-6], $
;[-20,-16.5,-13.8,-11.3,-9.2, -7.5,-6], $
;[-20,-16.5,-13.8,-11.3,-9.2, -7.5,-6]]
;store_data, 'bin_x_ey_ave', data = [$
;[-20,-16.5,-13.8,-11.3,-9.2, -7.5,-6], $
;[-20,-16.5,-13.8,-11.3,-9.2, -7.5,-6], $
;[-20,-16.5,-13.8,-11.3,-9.2, -7.5,-6], $
;[-20,-16.5,-13.8,-11.3,-9.2, -7.5,-6]]
;;; plot the 2 panels
;for i = 0, n_elements(vars_2_bin)-1 do begin
;	get_data, vars_2_bin(i), data = data
;	x_use = data.x
;	qtt_use = data.qtt
;	ae_use = data.ae
;	get_data, 'bin_x_'+vars_2_bin(i), data = x_bin_boundaries
;	if strcmp(ae_ranges, 'divide') then begin
;		;;; get ae_levels with divided values
;	endif
;	;; plot lines for different AEs
;	for j = 0, n_elements(ae_levels)-2 do begin
;		i_this = where((ae_use gt ae_levels(j)) and (ae_use lt ae_levels(j+1)))
;		x_this = x_use(i_this)
;		qtt_this = qtt_use(i_this)
;		if j eq 0 then add = 0 else add = 1
;		stat_plot, x_this, qtt_this, k_c = k_c, bin_boundaries = x_bin_boundaries(*,j), qtt_2_range = yranges(*,i), qtt_range = xrange, qtt_2_title = ytitles(i), qtt_title = xtitle, qtt_tickname = xticknames, kinbin = kinbin, bincntrs_out = bincenters, pm = 1, title = titles(i), bar_thick = l_thick, symsize = symsize, n_pst = n_pst, n_neg = n_neg, /no_write_pm, position = positions(*,n_elements(positions(0,*))-i-1), n_ttl = n_ttl, /noerase, add = add, color_all = color_ae(j)
;	endfor
;	xyouts, !x.crange(0)+0.86*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.86*(!y.crange(1)-!y.crange(0)), abc(i), /data, alignment = 0, charsize = 1.
;endfor

;;;;;; plot the two occurence rate plots
;; settings 
;;; choose number of AE levels
ae_levels = ae_3
color_ae = color_ae_3
vars_2_bin = ['num', 'dur']
yranges = [[0,5.9], [0, 249]] ;; for 0809 only
ytitles = ['# of DFBs per!c1000min orbit time', '']
;abc = ['(c)','(d)']
abc = ['(a)','(b)']
bin_range = [-30., 0.]
titles = ['DFB occurrence rate', '']
usersym, 1.*[-1,0,1,0], 1.*[0,1,0,-1], /fill, thick = l_thick ;; diamond
;;;;;;; set x bin boundaries
;;;; for all seasons
;store_data, 'bin_x_num', data = [$
;[-33, -30, -26, -21, -17.8, -14.5, -11.5,-10.5,-9.5, -8.5,-7.5,-6], $
;[-35, -33, -30, -22, -17.8, -14.5, -11.8,-10.5,-9.5, -8.5,-7.5,-6], $
;[-38, -35, -33, -30, -17, -14.5, -11.8,-10.5,-9.5, -8.5,-7.5,-6], $
;[-26, -23, -21, -17.8, -15,-13,-11.5,-10.5,-9.5, -8.5,-7.5,-6]]
;store_data, 'bin_x_dur', data = [$
;[-33, -30, -26, -21, -17.8, -14.5, -11.5,-10.5,-9.5, -8.5,-7.5,-6], $
;[-35, -33, -30, -22, -17.8, -14.5, -11.8,-10.5,-9.5, -8.5,-7.5,-6], $
;[-38, -35, -33, -30, -17, -14.5, -11.8,-10.5,-9.5, -8.5,-7.5,-6], $
;[-26, -23, -21, -17.8, -15,-13,-11.5,-10.5,-9.5, -8.5,-7.5,-6]]
;;;; for 0809 only
store_data, 'bin_x_num', data = [$
[-33, -30, -24, -21, -17.8, -14.5, -12.5,-10.5,-9.5, -8.5,-7.5,-6], $
[-40, -38, -35, -33, -30, -22, -17.8, -13.5, -11.2, -9.5, -8., -6], $
[-40,-38, -35, -33, -30, -17, -15., -13.8,-11.,-9.5, -8,-6], $
[-26, -23, -21, -17.8, -15,-13,-12.5,-11.,-9.5, -8.5,-7.5,-6]]
store_data, 'bin_x_dur', data = [$
[-30, -28, -24, -21, -17.8, -14.5, -12.5,-10.5,-9.5, -8.5,-7.5,-6], $
[-40, -38, -35, -33, -30, -22, -17.8, -13.5, -11.2, -9.5, -8., -6], $
[-40,-38, -35, -33, -30, -17, -15., -13.8,-11.,-9.5, -8,-6], $
[-26, -23, -21, -17.8, -15,-13,-11.5,-10.5,-9.5, -8.5,-7.5,-6]]
;;;; for 1011 only
store_data, 'bin_x_num_1011', data = [$
[-12,-10.5,-9.5, -8.5,-7.5,-6], $
[-12,-10.5,-9.5, -8.5,-7.5,-6], $
[-12,-10.5,-9.5, -8.5,-7.5,-6], $
[-11.5,-10.5,-9.5, -8.5,-7.5,-6]]
store_data, 'bin_x_dur_1011', data = [$
[-12,-10.5,-9.5, -8.5,-7.5,-6], $
[-12,-10.5,-9.5, -8.5,-7.5,-6], $
[-12,-10.5,-9.5, -8.5,-7.5,-6], $
[-11.5,-10.5,-9.5, -8.5,-7.5,-6]]
;; plot the 2 panels
for i = 0, n_elements(vars_2_bin)-1 do begin
	if i eq n_elements(vars_2_bin)-1 then begin
	    xticknames = ''
		xtitle = 'X!dGSM!n [R!dE!n]'
	endif else begin
	    xticknames = replicate(' ', 59)
		xtitle = ''
	endelse
	get_data, vars_2_bin(i), data = qtt_this
	get_data, 'bin_x_'+vars_2_bin(i), data = x_bin_boundaries
	if strcmp(season_suf, '_both') then get_data, 'bin_x_'+vars_2_bin(i)+'_1011', data = x_bin_boundaries_1011

	;plot, [0,0], [0,0], /nodata, xtitle = xtitle, ytitle = ytitles(i), title = titles(i), xstyle = 1, ystyle = 1, xrange = xrange, yrange = yranges(*,i), position = positions(*,n_elements(positions(0,*))-i-3), /noerase, xtickname = xticknames
	plot, [0,0], [0,0], /nodata, xtitle = xtitle, ytitle = ytitles(i), title = titles(i), xstyle = 1, ystyle = 1, xrange = xrange, yrange = yranges(*,i), position = positions(*,n_elements(positions(0,*))-i-1), /noerase, xtickname = xticknames
	for j = 0, n_elements(ae_levels)-2 do begin
		i_this_dfb = where((ae_ave_peak gt ae_levels(j)) and (ae_ave_peak lt ae_levels(j+1)))
		i_this_orbit = where((ae_orbit gt ae_levels(j)) and (ae_orbit lt ae_levels(j+1)))
		x_dfb_this = x_peak(i_this_dfb)
		sec_dfb_this = sec_dfb(i_this_dfb)
		x_orbit_this = pos_orbit(0, i_this_orbit)
		bin1d, x_dfb_this, sec_dfb_this, bin_range(0), bin_range(1), binsize, kinbin_dfb, bincntrs, mean_sec, std_th, med_th, flag4nodata = !values.f_nan, bin_boundaries = x_bin_boundaries(*,j)
		bin1d, transpose(x_orbit_this), transpose(x_orbit_this), bin_range(0), bin_range(1), binsize, kinbin_orbit, bincntrs, avrg_th, std_th, med_th, flag4nodata = !values.f_nan, bin_boundaries = x_bin_boundaries(*,j)
		print, kinbin_dfb
		print, mean_sec
		print, kinbin_orbit
		rate = double(kinbin_dfb)/kinbin_orbit*1000
		i_few = where(kinbin_dfb lt k_c, j_few)
		if j_few gt 0 then rate(i_few) = !values.f_nan
		print, rate
		ttl_sec = rate*mean_sec
		if strcmp(vars_2_bin(i), 'num') then plt_qtt = rate else plt_qtt = ttl_sec
		oplot, bincntrs, plt_qtt, color = color_ae(j)
		oplot, bincntrs, plt_qtt, psym = 8, color = color_ae(j)
		;;; plot other two season events
		if strcmp(season_suf, '_both') and strcmp(vars_2_bin(i), 'num') then begin
			i_this_dfb = where((ae_ave_peak_1011 gt ae_levels(j)) and (ae_ave_peak_1011 lt ae_levels(j+1)))
			i_this_orbit = where((ae_orbit_1011 gt ae_levels(j)) and (ae_orbit_1011 lt ae_levels(j+1)))
			x_dfb_this = x_peak_1011(i_this_dfb)
			sec_dfb_this = sec_dfb_1011(i_this_dfb)
			x_orbit_this = pos_orbit_1011(0, i_this_orbit)
			bin1d, x_dfb_this, sec_dfb_this, bin_range(0), bin_range(1), binsize, kinbin_dfb, bincntrs, mean_sec, std_th, med_th, flag4nodata = !values.f_nan, bin_boundaries = x_bin_boundaries_1011(*,j)
			bin1d, transpose(x_orbit_this), transpose(x_orbit_this), bin_range(0), bin_range(1), binsize, kinbin_orbit, bincntrs, avrg_th, std_th, med_th, flag4nodata = !values.f_nan, bin_boundaries = x_bin_boundaries_1011(*,j)
			print, kinbin_dfb
			print, mean_sec
			print, kinbin_orbit
			rate = double(kinbin_dfb)/kinbin_orbit*1000
			i_few = where(kinbin_dfb lt k_c, j_few)
			if j_few gt 0 then rate(i_few) = !values.f_nan
			print, rate
			ttl_sec = rate*mean_sec
			if strcmp(vars_2_bin(i), 'num') then plt_qtt = rate else plt_qtt = ttl_sec
			oplot, bincntrs, plt_qtt, color = color_ae(j), line = 0
			oplot, bincntrs, plt_qtt, psym = 4, color = color_ae(j), symsize = 0.6
			y_0809 = 0.45
			y_1011 = 0.38
			x_leg = 0.7
			x_leg_inc = 0.025
			y_leg_inc = -0.02
			oplot, [!x.crange(0)+x_leg*(!x.crange(1)-!x.crange(0)), !values.f_nan], [!y.crange(0)+y_0809*(!y.crange(1)-!y.crange(0)), !values.f_nan], psym = 8
			oplot, [!x.crange(0)+x_leg*(!x.crange(1)-!x.crange(0)), !values.f_nan], [!y.crange(0)+y_1011*(!y.crange(1)-!y.crange(0)), !values.f_nan], psym = 4, symsize = 0.6
			xyouts, !x.crange(0)+(x_leg+x_leg_inc)*(!x.crange(1)-!x.crange(0)), !y.crange(0)+(y_0809+y_leg_inc)*(!y.crange(1)-!y.crange(0)), '08-09', charsize = 0.7
			xyouts, !x.crange(0)+(x_leg+x_leg_inc)*(!x.crange(1)-!x.crange(0)), !y.crange(0)+(y_1011+y_leg_inc)*(!y.crange(1)-!y.crange(0)), '10-11', charsize = 0.7
		endif
	endfor
	xyouts, !x.crange(0)+0.86*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.86*(!y.crange(1)-!y.crange(0)), abc(i), /data, alignment = 0, charsize = 1.
	if strcmp(vars_2_bin(i), 'dur') then begin
		;; after dur because using normal
		x_leg = 0.96
		y_leg = 0.83
		y_leg_inc = -0.04
		charsize = 0.8
		xyouts, x_leg, y_leg,  'AE>300nT', charsize = charsize, color = color_ae(2), /normal, align = 1
		xyouts, x_leg, y_leg+y_leg_inc,  '150nT<AE<300nT', charsize = charsize, color = color_ae(1), /normal, align = 1
		xyouts, x_leg, y_leg+2*y_leg_inc,  'AE<150nT', charsize = charsize, color = color_ae(0), /normal, align = 1
	endif
	if strcmp(vars_2_bin(i), 'dur') then begin
		;oplot, [-17.7, -17.7], !y.crange, line = 1
		;xyouts, -18.5, 0.02*(!y.crange(1)-!y.crange(0)), 'DFB seconds!care unrealistic!cto this side of!cthe dotted line', orientation = 90., charsize = 0.5, alignment = 0, /data
		if strcmp(season_suf, '0809') then x_color = 0.03 else x_color = 0.55
		xyouts, !x.crange(0)+x_color*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.15*(!y.crange(1)-!y.crange(0)), 'Colors & symbols!chave the same!cmeanings!cas in panel a.', charsize = 0.5, /data
	endif
endfor
;;; write the ytitle of panel b
xyouts, -0.015, 0.27, 'Seconds of DFBs per !c1000min orbit time', orientation = 90, alignment = 0.5, /normal
pclose

stop
end
