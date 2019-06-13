pro main_ratio_plot
;;; plot dfb/bbf ratios.
;;; run main_ratio_indiv and main_ratio_acum first

thm_init
@declare

;;; plot specifics
xsize = 4
ysize = 10
left_margin = 0.15
right_margin = 0.01
top_margin = 0.05
bot_margin = 0.05
vspace = 0.007
n_panels = 5
xrange = [0., 1.]
;;; the increment of last panel
bottom_inc = 0.05

loc_suf = ''

lim_suf = ''

dtrd_suf = ''
;dtrd_suf = '_vexb'
if strcmp(dtrd_suf, '_vexb') then pref_e = '' else pref_e = delta_letter

ratio_t_flux_all = datain_simple(save_folder+'/ratio_t_flux'+dtrd_suf+loc_suf+lim_suf+'.dat', dim=9, type='double')
number_all  = datain_simple(save_folder+'/number'+loc_suf+lim_suf+'.dat', dim=9, type='double')
ratio_bz_all  = datain_simple(save_folder+'/ratio_bz'+loc_suf+lim_suf+'.dat', dim=9, type='double')
ratio_vx_all  = datain_simple(save_folder+'/ratio_vx'+loc_suf+lim_suf+'.dat', dim=9, type='double')
ratio_ey_all  = datain_simple(save_folder+'/ratio_ey'+dtrd_suf+loc_suf+lim_suf+'.dat', dim=9, type='double')

cri_ddts = ratio_t_flux_all(*,0)
ratio_secs = ratio_t_flux_all(*,1)
ratio_flux = ratio_t_flux_all(*,2)
n_dfbs_in_bbf = number_all(*, 1:*)
ratio_bz = ratio_bz_all(*, 1:*)
ratio_vx = ratio_vx_all(*, 1:*)
ratio_ey = ratio_ey_all(*, 1:*)

store_data, 'ratio_n', data = n_dfbs_in_bbf 
store_data, 'ratio_bz', data = ratio_bz
store_data, 'ratio_vx', data = ratio_vx
store_data, 'ratio_ey', data = ratio_ey

;;;; the attributes of indiv plots
indiv_ratios = ['ratio_n', 'ratio_bz', 'ratio_vx', 'ratio_ey']
ytitles = ['# of DFBs in each BBF', jiao_l+'B!S!dz!R!UDFB!n'+jiao_r+'/'+jiao_l+'B!S!dz!R!UBBF!n'+jiao_r, jiao_l+'V!S!dx!R!UDFB!n'+jiao_r+'/'+jiao_l+'V!S!dx!R!UBBF!n'+jiao_r, jiao_l+pref_e+'E!S!dy,DSL'+"'"+'!R!UDFB   !n'+jiao_r+'/'+jiao_l+pref_e+'E!S!dy,DSL'+"'"+'!R!UBBF   !n'+jiao_r]
abc = ['(a)', '(b)', '(c)', '(d)'] 
titles = ['Attributes of DFBs Embedded in BBFs', '', '', '']
;yranges = [[0., 17.], [0., 3.], [0., 5.], [0., 15.]] ;; for mean-median plot
yranges = [[0., 14.99], [0., 2.4], [0., 3.3], [0., 5.99]] ;; for quartiles plot
y_abc = [0.83, 0.15, 0.15, 0.15]
error_scales = [1., 10., 1., 10.]
std_colors = [2,6,2,6]

;;; get the positions of the panels
positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [top_margin, bot_margin], space = [0., vspace], height = height)
inc_each = bottom_inc/n_elements(positions(0,*))
height_new = height-inc_each
positions(3,0) = positions(3,0)+bottom_inc
positions(1,1) = positions(1,0)+height+bottom_inc+vspace 
positions(3,1) = positions(1,1)+height_new
for i = 2, 4 do begin
	positions([1,3],i) = positions([1,3],i-1)+height_new+vspace
endfor

popen, pic_folder+'/ratios'+dtrd_suf+loc_suf+lim_suf
print_options,xsize=xsize, ysize=ysize ;; use this for single plot
;;;;;; time and flux ratio
usersym, 1.3*[-1,0,1,0], 1.3*[0,1,0,-1], /fill, thick = l_thick ;; diamond
plot, cri_ddts, ratio_secs, position = positions(*,0), xrange = xrange, yrange = [0., 1.], xtitle = 'DFB threshold of dB!S!dz!R!Usm!n/dt [nT/s]', ytitle = 'DFB/BBF Ratio', yticklen = 1, ygridstyle = 1, thick = l_thick
oplot, cri_ddts, ratio_secs, psym = 8;, color = 2
oplot, cri_ddts, ratio_flux, color = 6, thick = l_thick
oplot, cri_ddts, ratio_flux, psym = 8, color = 6
;;; add a 0.5 line for future use
oplot, [0.5, 0.5], !y.crange, line = 1
xyouts, 0.38, 0.3, 'duration', /data, align = 1
xyouts, 0.57, 0.64, 'flux transport', /data, color = 6
xyouts, 0.87, 0.87, '(e)', charsize = 1.2, /data

;;;; plot "individial ratios"
for i = 0, n_elements(indiv_ratios)-1 do begin
	get_data, indiv_ratios(i), data = qtt_this
	position_this = positions(*, 4-i)
	stat_plot_disc, cri_ddts, qtt_this, qtt_range = xrange, qtt_tickname = replicate(' ',6), position = position_this, qtt_2_range = yranges(*,i), qtt_2_title = ytitles(i), /noerase, error_scale = error_scales(i), color_std = std_colors(i), title = titles(i), /no_dots, color_quar = 2, color_med = 6, /no_mean, sym_thick = 2
	xyouts, 0.9, y_abc(i)*!y.crange(1), abc(i), align = 0.5, charsize = 1.2
endfor
pclose

stop
end
