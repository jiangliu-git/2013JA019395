function tr_flow_burst, bbf_trange_input, probe = probe, cri_fb = cri_fb, v_type = v_type, v_folder = v_folder, b_folder = b_folder, e_folder = e_folder, secs_fb = secs_fb, percent_t_fb = percent_t_fb, rtrange = rtrange, rt_pre = rt_pre, rt_length = rt_length, gsm = gsm
;; get the timeranges of the flow burst
;; bbf_trange_input: the trange of the bbf, for a single event
;; cri_fb: criterion for flow burst, in km/h, default is 400
;; b_folder: the folder of fgs, if set and e_folder not set, use vxB for E.
;; e_folder: the folder of EFS_DSL, if set use EFI measurement for E
;; gsm: only effective when both efield and e_folder are set, indicates transforming the EFI data into gsm.
;; rtrange: the time range to remove the efi offset, only useful when e_folder is set.
;; output:
;; secs_fb: the total duration of the FB tranges
;;; note: When using efield and not setting GSM, only EyDSL is used to select

if ~keyword_set(v_type) then v_type = 'vi'
if ~keyword_set(cri_fb) then begin
	if keyword_set(e_folder) or keyword_set(b_folder) then cri_fb = 4. else cri_fb = 400. ;; mV/m and km/s
endif
if strcmp(v_type, 'vi') then v_str = 'ptix_velocity'
if strcmp(v_type, 'vperp') or strcmp(v_type, 'viperp') then v_str = 'ptix_vperp'
b_str = 'fgs'

bbf_trange = time_double(bbf_trange_input)
if (~keyword_set(rtrange)) and keyword_set(rt_length) then rtrange = [bbf_trange(0)-60.*(rt_pre+0.5*rt_length), bbf_trange(0)-60.*(rt_pre-0.5*rt_length)]

fb_tranges = [!values.f_nan, !values.f_nan]

;;;; load data
if keyword_set(e_folder) then begin
	del_data, 'th'+probe+'_efs*'
	if keyword_set(gsm) then begin
		e_name = 'th'+probe+'_efs_gsm'
		load_efi_data, trange = bbf_trange, probe = probe, rtrange = rtrange, e_folder = e_folder, b_folder = b_folder, /tclip
	endif else begin
		e_name = 'th'+probe+'_efs_dsl'
		load_efi_data, trange = bbf_trange, probe = probe, rtrange = rtrange, e_folder = e_folder, /dsl, /tclip
	endelse
	if tv_exist(e_name+'_tclip') then begin
		get_data, e_name+'_tclip', data=e_data
		t = e_data.x
		Ex = e_data.y(*,0)
		Ey = e_data.y(*,1)
		Ez = e_data.y(*,2)
		Ettl = sqrt(Ex^2+Ey^2+Ez^2)
		;; correct for up-side-down spacecrafts thb and thc for DSL
		if ~keyword_set(gsm) then begin
			if strcmp(probe, 'b') or strcmp(probe, 'c') then Ey=-Ey
		endif
		;; the data to select from
		if keyword_set(gsm) then data_select = Ettl else data_select = abs(Ey)
	endif else begin
		print, 'TR_FLOW_BURST: No data!'
		secs_fb = !values.f_nan
		percent_t_fb = !values.f_nan
		return, !values.f_nan
	endelse
endif else begin
	del_data, 'th'+probe+'_'+v_str+'_gsm*'
	load_bin_data, trange = bbf_trange, probe = probe, datatype = v_type, datafolder = v_folder, /tclip
	if tv_exist('th'+probe+'_'+v_str+'_gsm_tclip') then begin
		get_data, 'th'+probe+'_'+v_str+'_gsm_tclip', data = v_data
		t = v_data.x
		vx = v_data.y(*,0)
		vy = v_data.y(*,1)
		vz = v_data.y(*,2)
		vttl = sqrt(vx^2+vy^2+vz^2)
		if keyword_set(b_folder) then begin
			del_data, 'th'+probe+'_'+b_str+'_gsm*'
			load_bin_data, trange = bbf_trange, probe = probe, datatype = b_str, datafolder = b_folder, /tclip
			if tv_exist('th'+probe+'_'+b_str+'_gsm_tclip') then begin
				tinterpol_mxn,'th'+probe+'_'+b_str+'_gsm_tclip','th'+probe+'_'+v_str+'_gsm_tclip',newname='th'+probe+'_'+b_str+'_gsm_int'
				get_data, 'th'+probe+'_'+b_str+'_gsm_int', data = b_data
				bx = b_data.y(*,0)
				by = b_data.y(*,1)
				bz = b_data.y(*,2)
				Ex = 0.001*(vz*By-vy*Bz)
				Ey = 0.001*(vx*Bz-vz*Bx)
				Ez = 0.001*(vy*Bx-vx*By)
				Ettl = sqrt(Ex^2+Ey^2+Ez^2)
				data_select = Ettl
			endif else begin
				print, 'TR_FLOW_BURST: No data!'
				secs_fb = !values.f_nan
				percent_t_fb = !values.f_nan
				return, !values.f_nan
			endelse
		endif else data_select = vttl
	endif else begin
		print, 'TR_FLOW_BURST: No data!'
		secs_fb = !values.f_nan
		percent_t_fb = !values.f_nan
		return, !values.f_nan
	endelse
endelse

;; find the flow bursts
i_over = where(data_select gt cri_fb, ni_over)
if ni_over gt 1 then begin
	j_duan = where(i_over(1:ni_over-1)-i_over(0:ni_over-2) gt 1, nj_duan)
	if nj_duan gt 0 then begin
		i_fb_start = [i_over(0), i_over(j_duan+1)]
		i_fb_end = [i_over(j_duan), i_over(ni_over-1)]
		fb_tranges = t([transpose(i_fb_start), transpose(i_fb_end)])
	endif else begin
		i_fb_start = i_over(0)
		i_fb_end = i_over(ni_over-1)
		fb_tranges = t([i_fb_start, i_fb_end])
	endelse
endif else begin
	if ni_over eq 1 then begin
		fb_tranges = t([i_over, i_over])
	endif else begin
		print, 'No flow burst in the given time range!'
	endelse
endelse
secs_fb = total(fb_tranges(1,*)-fb_tranges(0,*))
percent_t_fb = secs_fb/(bbf_trange(1)-bbf_trange(0))
return, fb_tranges
end
