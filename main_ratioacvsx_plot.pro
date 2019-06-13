pro main_ratioacvsx_plot
;;; plot for dfb/bbf ratio vs x. Run main_ratioacvsx.pro first.

thm_init

@declare

listname = 'bbf_ps_list_original2'

;;; for x dependence one should only use 0809 events.
;season_suf = ''
season_suf = '_0809'

;;; set the minimum bin number
k_c = 5

;;;; panel specifics
l_thick = 2.
xsize = 3.5
ysize = 4
left_margin = 0.15
right_margin = 0.01
top_margin = 0.05
bot_margin = 0.05
vspace = 0.017
n_panels = 2
rate_bar = 0.1
abc = ['(a)', '(b)']
charsize_name = 0.8
xdur = [0.2,0.8]
ydur = [0.3,0.3]
charsize_met = 0.7
xf = [0.5,0.8]
yf = [0.64,0.62]
method_name = ['fixed', 'X-dependent']
xmet = [0.85,0.5]
ymet = [0.73,0.1]
xrange = [-5, -29.99]
usersym, 1.1*[-1,0,1,0], 1.1*[0,1,0,-1], /fill, thick = l_thick ;; diamond

data_original = datain_simple(save_folder+'/ratio_t_flux_vsx_bz'+season_suf+'.dat', dim=4, type='double')
data_assume = datain_simple(save_folder+'/ratio_t_flux_vsx_bz_trend_xdur'+season_suf+'.dat', dim=4, type='double')

bin_bound_original = [-30., -23, -18.2, -16.,-13.4,-11,-9,-6]
bin_bound_assume = [-30., -24, -20., -16,-13,-11,-9,-6]

store_data, 'original', data = {data:data_original, bin:bin_bound_original}
store_data, 'assume', data = {data:data_assume, bin:bin_bound_assume}

data_plot = ['original', 'assume']

positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [top_margin, bot_margin], space = [0., vspace], height = height)

;;;;;;; plot the dfb-bbf ratio dependence on cri_ddt

popen, pic_folder+'/ratio_vs_x'+season_suf
print_options,xsize=xsize, ysize=ysize ;; use this for single plot
for i = 0, n_elements(data_plot)-1 do begin
	if i eq 0 then title = 'DFB flux transport importance' else title = ''
	if i eq n_elements(data_plot)-1 then begin
	    xticknames = ''
		xtitle = 'X!dGSM!n [R!dE!n]'
		method_str = 'From '+method_name(i)+' DFB criteria & durations'
	endif else begin
	    xticknames = replicate(' ', 59)
		xtitle = ''
		method_str = 'From '+method_name(i)+'!cDFB criterion!c& duration'
	endelse

	get_data, data_plot(i), data = this
	data = this.data
	binb = this.bin

	;;; remove the 30 RE bin for the bz_trend plot
	if i eq n_elements(data_plot)-1 then begin
		data(*, 0) = !values.f_nan
		binb(0) = !values.f_nan
	endif
	;;; remove unpopulated bins
	kinbin = data(3,*)
	print, kinbin
	i_bad = where(kinbin le k_c, j_bad)
	if j_bad gt 0 then begin
		data(*,i_bad) = !values.f_nan
	endif
	ratio_secs = data(1,*)
	ratio_flux = data(2,*)
	loc_values = data(0,*)

	plot, loc_values, ratio_secs, xrange = xrange, yrange = [0., 1.], xtitle = xtitle, xstyle = 1, ystyle = 1, title = title, xtickname = xticknames, ytitle = 'DFB/BBF Ratio', yticklen = 1, ygridstyle = 1, thick = l_thick, position = positions(*,n_elements(data_plot)-1-i), /noerase
	oplot, loc_values, ratio_secs, psym = 8
	oplot, loc_values, ratio_flux, color = 6, thick = l_thick
	oplot, loc_values, ratio_flux, psym = 8, color = 6
	;; plot the bin boundaries
	for j = 0, n_elements(binb)-1 do begin
		oplot, [binb(j), binb(j)], [!y.crange(1), !y.crange(1)-rate_bar*(!y.crange(1)-!y.crange(0))], thick = 0.1
	endfor

	xyouts, !x.crange(0)+0.86*(!x.crange(1)-!x.crange(0)), !y.crange(0)+0.86*(!y.crange(1)-!y.crange(0)), abc(i), /data, alignment = 0, charsize = 1.
	xyouts, !x.crange(0)+xdur(i)*(!x.crange(1)-!x.crange(0)), !y.crange(0)+ydur(i)*(!y.crange(1)-!y.crange(0)), charsize = charsize_name, 'duration', /data, align = 0.5
	xyouts, !x.crange(0)+xf(i)*(!x.crange(1)-!x.crange(0)), !y.crange(0)+yf(i)*(!y.crange(1)-!y.crange(0)), charsize = charsize_name, 'flux transport', /data, align = 0.5, color = 6
	xyouts, !x.crange(0)+xmet(i)*(!x.crange(1)-!x.crange(0)), !y.crange(0)+ymet(i)*(!y.crange(1)-!y.crange(0)), charsize = charsize_met, method_str, /data, align = 0.5
endfor
pclose

stop
end
