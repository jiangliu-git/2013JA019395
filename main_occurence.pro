pro main_occurence
;;; plot the occurence rate of DFB events.

thm_init
;; generate bin plots, like occurence rates at different locations, etc.
computer='I:'
;computer='/home/jliu'

@declare

suffix = '_fgs'
dur_suf = '_40s'
dtrd_suf = ''

;;; whether use only positive/negative transport events, told by Ey_dsl
;sign_use_suf = ''
sign_use_suf = '_positive'
;sign_use_suf = '_negative'

;;; x binsizes of the ae regions
store_data, 'bin0', data = [-21, -17.8, -15,-13,-11.5,-10.5,-9.5, -8.5,-7.5,-6]
store_data, 'bin1', data = [-21, -18, -15, -13, -11.5,-10.5,-9.5, -8.5,-7.5,-6]
store_data, 'bin2', data = [-21, -18, -15, -13, -11.5,-10.5,-9.5, -8.5,-7.5,-6]
store_data, 'bin3', data = [-21, -17.8, -15,-13,-11.5,-10.5,-9.5, -8.5,-7.5,-6]

ae_sep = [0., 100, 200, 400, 1000000000.]

;; plot specifics
xsize = 4
ysize = 3
xplotrange = [-6, -22]
yplotrange = [0, 6]
xtitle = 'X [R!dE!n]'
ytitle = '# of DFBs/1000min orbit time'
colors = [97, 82, 62, 40]
usersym, 1.*[-1,0,1,0], 1.*[0,1,0,-1], /fill, thick = l_thick ;; diamond

case sign_use_suf of
'': title = 'All Events'
'_positive': title = 'Positive '+delta_letter+'E!S!dy,DSL'+"'"+'!R!Upeak!n Events'
'_negative': title = 'Negative '+delta_letter+'E!S!dy,DSL'+"'"+'!R!Upeak!n Events'
endcase

listname = 'dfb_list_lead_tail'+suffix+'.txt'
events = load_list(listname, folder = listfolder)

ae_ave = datain_simple(save_folder+'/pseudo_ae_ave'+suffix+'.dat', dim=1, type='double')
pos = datain_simple(save_folder+'/pos'+suffix+'.dat', dim=3, type='double')
ey_peak = datain_simple(save_folder+'/dfb_ey_peak'+dur_suf+suffix+dtrd_suf+'.dat', dim=1, type='double')
orbit_pos_ae = datain_simple(save_folder+'/orbit_xyz_ae'+suffix+'.dat', dim=4, type='double')
pos_orbit = orbit_pos_ae(0:2,*)
ae_orbit = orbit_pos_ae(3,*)

;;; begin working
case sign_use_suf of
'': i_use = intarr(n_elements(ae_ave))
'_positive': i_use = where(ey_peak gt 0)
'_negative': i_use = where(ey_peak lt 0)
endcase

ae_ave = ae_ave(*,i_use)
pos = pos(*,i_use)

bin_range = [0, -30]

popen, pic_folder+'/dfb_occurence'+dtrd_suf
print_options,xsize=xsize, ysize=ysize ;; use this for single plot
plot, [0,0], [0,0], /nodata, xtitle = xtitle, ytitle = ytitle, xstyle = 1, title = title, xrange = xplotrange, yrange = yplotrange
for i = 0, n_elements(ae_sep)-2 do begin
	j_this_dfb = where((ae_ave gt ae_sep(i)) and (ae_ave lt ae_sep(i+1)))
	j_this_orbit = where((ae_orbit gt ae_sep(i)) and (ae_orbit lt ae_sep(i+1)))
	x_dfb_this = pos(0, j_this_dfb)
	x_orbit_this = pos_orbit(0, j_this_orbit)
	get_data, 'bin'+strcompress(string(i), /remove), data = x_bins
	bin1d, transpose(x_dfb_this), transpose(x_dfb_this), bin_range(0), bin_range(1), binsize, kinbin_dfb, bincntrs, avrg_th, std_th, med_th, flag4nodata = !values.f_nan, bin_boundaries = x_bins
	bin1d, transpose(x_orbit_this), transpose(x_orbit_this), bin_range(0), bin_range(1), binsize, kinbin_orbit, bincntrs, avrg_th, std_th, med_th, flag4nodata = !values.f_nan, bin_boundaries = x_bins
	print, kinbin_dfb
	print, kinbin_orbit
	rate = double(kinbin_dfb)/kinbin_orbit*1000
	print, rate
	oplot, bincntrs, rate, color = colors(i)
	oplot, bincntrs, rate, psym = 8, color = colors(i)
endfor
pclose

end
