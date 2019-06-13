pro bbf_flux_cal, bbf_trange_input, dfb_tranges_input, probe = probe, b_type = b_type, b_folder = b_folder, fgs_folder = fgs_folder, e_folder = e_folder, gsm = gsm, v_type = v_type, v_folder = v_folder, flux_bbf = flux_bbf, flux_dfb = flux_dfb, x_only = x_only, eflux_bbf = eflux_bbf, eflux_dfb= eflux_dfb, secs_dfb = secs_dfb, percent_t_dfb = percent_t_dfb, dfb_tranges_actual = dfb_tranges_actual, rtrange = rtrange, rt_pre = rt_pre, rt_length = rt_length, vexb_dsl_folder = vexb_dsl_folder
;; input:
;; bbf_trange_input: the trange of the bbf, for a single event
;; dfb_trange_input: the tranges of the DFBs inside the BBF, |t_start|t_end|, can be all, the routine will pick those inside BBFs, can be overlapping, the routine will solve overlapping.
;; gsm: only useful when calculating from EFI data, which has DSL as default; if use other data gsm is default and no need to set.
;; e_folder: the folder of EFS_DSL, if set will be calculated using efi data
;; output:
;; flux_bbf/dfb: the magnetic flux transport by the bbf and the DFBs inside, in Wb/RE, is integral of VxBz-VzBx, if x_only is set, then VxBz only
;; eflux: energy flux.
;; dfb_tranges_actual: the actual DFB tranges from input that are within the BBF
;; secs_dfb: the total duration of the DFB tranges
;; rt_pre, rt_length: in minutes, if rtrange set, no need to set
;; How to calculate flux is on Angelopoulos et al. [1994] Pg 16
;; To be precise, if one point of the Data is NaN then the entire flux transport is NaN

if keyword_set(fgs_folder) then begin
	b_type = 'fgs'
	b_folder = fgs_folder
endif

if ~keyword_set(b_type) then b_type = 'fgs'
if ~keyword_set(v_type) then v_type = 'vi'

b_str = b_type
if strcmp(v_type, 'vi') then v_str = 'ptix_velocity'
if strcmp(v_type, 'vperp') or strcmp(v_type, 'viperp') then v_str = 'ptix_vperp'

bbf_trange = time_double(bbf_trange_input)
dfb_tranges = time_double(dfb_tranges_input)
secs_bbf = bbf_trange(1)-bbf_trange(0)

if (~keyword_set(rtrange)) and keyword_set(rt_length) then rtrange = [bbf_trange(0)-60.*(rt_pre+0.5*rt_length), bbf_trange(0)-60.*(rt_pre-0.5*rt_length)]

;;;;;;;;;;;;;;; calculate the transport by the BBF ;;;;;;;;;;;;;;;
if keyword_set(e_folder) then begin
	;;;;;; use EFI
	del_data, 'th'+probe+'_fgs*'
	del_data, 'th'+probe+'_efs*'
	del_data, 'th'+probe+'_intEy*'
	if keyword_set(gsm) then begin
		e_name = 'th'+probe+'_efs_gsm'
		load_efi_data, trange = bbf_trange, probe = probe, rtrange = rtrange, e_folder = e_folder, b_folder = fgs_folder, /tclip, vexb_dsl_folder = vexb_dsl_folder
	endif else begin
		e_name = 'th'+probe+'_efs_dsl'
		load_efi_data, trange = bbf_trange, probe = probe, rtrange = rtrange, e_folder = e_folder, /dsl, /tclip, vexb_dsl_folder = vexb_dsl_folder
	endelse
	if tv_exist(e_name+'_tclip') then begin
		get_data, e_name+'_tclip', data=e_data
		t = e_data.x
		Ey = e_data.y(*,1)
		;; correct for up-side-down spacecrafts thb and thc for DSL
		if ~keyword_set(gsm) then begin
			if strcmp(probe, 'b') or strcmp(probe, 'c') then Ey=-Ey
		endif
	endif else begin
		print, 'BBF_FLUX_CAL: No data!'
		flux_bbf = !values.f_nan
		flux_dfb = !values.f_nan
		eflux_bbf = !values.f_nan
		eflux_dfb = !values.f_nan
		secs_dfb = !values.f_nan
		percent_t_dfb = !values.f_nan
		dfb_tranges_actual = [!values.f_nan, !values.f_nan]
		return
	endelse
endif else begin
	;;;;;; use VxB
	del_data, 'th'+probe+'_'+b_str+'_gsm*'
	del_data, 'th'+probe+'_'+v_str+'_gsm*'
	del_data, 'th'+probe+'_intEy*'
	load_bin_data, trange = bbf_trange, probe = probe, datatype = b_type, datafolder = b_folder, /tclip
	load_bin_data, trange = bbf_trange, probe = probe, datatype = v_type, datafolder = v_folder, /tclip
	if tv_exist('th'+probe+'_'+b_str+'_gsm_tclip') and tv_exist('th'+probe+'_'+v_str+'_gsm_tclip') then begin
		tinterpol_mxn,'th'+probe+'_'+b_str+'_gsm_tclip','th'+probe+'_'+v_str+'_gsm_tclip',newname='th'+probe+'_'+b_str+'_gsm_int'
		get_data, 'th'+probe+'_'+b_str+'_gsm_int', data = b_data
		get_data, 'th'+probe+'_'+v_str+'_gsm_tclip', data = v_data
		t = v_data.x
		bx = b_data.y(*,0)
		by = b_data.y(*,1)
		bz = b_data.y(*,2)
		vx = v_data.y(*,0)
		vy = v_data.y(*,1)
		vz = v_data.y(*,2)
		if keyword_set(x_only) then Ey = 0.001*vx*Bz else Ey = 0.001*(vx*Bz-vz*Bx) ;; in mV/m
	endif else begin
		print, 'BBF_FLUX_CAL: No data!'
		flux_bbf = !values.f_nan
		flux_dfb = !values.f_nan
		eflux_bbf = !values.f_nan
		eflux_dfb = !values.f_nan
		secs_dfb = !values.f_nan
		percent_t_dfb = !values.f_nan
		dfb_tranges_actual = [!values.f_nan, !values.f_nan]
		return
	endelse
endelse


;; calculate the flux of the BBF
dt = t(1:*)-t(0:n_elements(t)-2)
flux_tran = [0,total(Ey(1:*)*dt, /cum)]
store_data, 'th'+probe+'_intEy', data = {x:t, y:flux_tran}
flux_bbf = flux_tran(n_elements(flux_tran)-1)
;eflux_tran = 
;eflux_bbf = 

;;;;;;;;;;;;;;; look at DFBs ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; manipulate tranges
dfb_tranges_actual = ranges_compress(dfb_tranges, limit = bbf_trange)
if ~finite(dfb_tranges_actual(0)) then begin
	flux_dfb = 0.
	eflux_dfb = 0.
	secs_dfb = 0.
	percent_t_dfb = 0.
	return
end

secs_dfb = total(dfb_tranges_actual(1,*)-dfb_tranges_actual(0,*))
percent_t_dfb = secs_dfb/secs_bbf

;;;;;; calculate the flux transport within the tranges
flux_dfb = 0.
eflux_dfb = 0.
for i = 0, n_elements(dfb_tranges_actual(0,*))-1 do begin
	;; time clip the already loaded data when calculating transport of BBF
	del_data, 'th'+probe+'_intEy_tclip'
	time_clip, 'th'+probe+'_intEy', dfb_tranges_actual(0,i), dfb_tranges_actual(1,i)
	if tv_exist('th'+probe+'_intEy_tclip') then begin
		get_data, 'th'+probe+'_intEy_tclip', data = data
		intEy = data.y
		flux_dfb = flux_dfb+intEy(n_elements(intEy)-1)-intEy(0)
	endif
	;;;; eflux: same round, to be built
endfor
end
