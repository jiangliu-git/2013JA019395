function bflux_cal, tranges_input, probe_list = probe_list, fluxtypes = fluxtypes, b_type = b_type, b_folder = b_folder, fgs_folder = fgs_folder, e_folder = e_folder, gsm = gsm, v_type = v_type, v_folder = v_folder, flux_dfb = flux_dfb, x_only = x_only, rtranges = rtranges, rt_pre = rt_pre, rt_length = rt_length, peak_ey = peak_ey
;;; calculate the bflux of specific types of the input tranges 
;;; fluxtypes: an arry specifying the types of flux you want, including: ['vxbz', 'vzbx', 'plasma', 'ey']. plasma = vxbz-vzbx
;;; efolder, gsm, rtrangs only useful for ey
;;; peak_ey: return the peak value of the "ey" used for flux transport calculation; same dimension as the return value
;;; return value is the fluxes, the same length of tranges_input
if keyword_set(fgs_folder) then begin
	b_type = 'fgs'
	b_folder = fgs_folder
endif

if ~keyword_set(b_type) then b_type = 'fgs'
if ~keyword_set(v_type) then v_type = 'vi'

b_str = b_type
if strcmp(v_type, 'vi') then v_str = 'ptix_velocity'
if strcmp(v_type, 'vperp') or strcmp(v_type, 'viperp') then v_str = 'ptix_vperp'

tranges = time_double(tranges_input)

if (~keyword_set(rtranges)) and keyword_set(rt_length) then rtranges = [tranges(0)-60.*(rt_pre+0.5*rt_length), tranges(0)-60.*(rt_pre-0.5*rt_length)]

bflux_out = dblarr(n_elements(fluxtypes), n_elements(tranges_input(0,*)))
peak_ey = dblarr(n_elements(fluxtypes), n_elements(tranges_input(0,*)))

for i = 0, n_elements(fluxtypes)-1 do begin
	for j = 0, n_elements(tranges_input(0,*))-1 do begin
		probe = probe_list(j)
		trange = tranges(*,j)
		if strcmp(fluxtypes(i), 'ey') then begin
			;;;;;; use EFI
			if keyword_set(e_folder) then begin
				del_data, 'th'+probe+'_fgs*'
				del_data, 'th'+probe+'_efs*'
				del_data, 'th'+probe+'_intEy*'
				if keyword_set(gsm) then begin
					e_name = 'th'+probe+'_efs_gsm'
					load_efi_data, trange = trange, probe = probe, rtrange = rtranges(*,j), e_folder = e_folder, b_folder = fgs_folder, /tclip
				endif else begin
					e_name = 'th'+probe+'_efs_dsl'
					load_efi_data, trange = trange, probe = probe, rtrange = rtrange(*,j), e_folder = e_folder, /dsl, /tclip
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
					print, 'FLUX_CAL: No data!'
					t = !values.f_nan
				endelse
			endif else begin
				print, 'FLUX_CAL: Efolder not set!'
				t = !values.f_nan
			endelse
		endif else begin
			;;;;;; use VxB
			del_data, 'th'+probe+'_'+b_str+'_gsm*'
			del_data, 'th'+probe+'_'+v_str+'_gsm*'
			del_data, 'th'+probe+'_intEy*'
			load_bin_data, trange = trange, probe = probe, datatype = b_type, datafolder = b_folder, /tclip
			load_bin_data, trange = trange, probe = probe, datatype = v_type, datafolder = v_folder, /tclip
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
				case fluxtypes(i) of
				'vxbz': Ey = vx*Bz
				'vzbx': Ey = -vz*Bx
				'plasma': Ey = vx*Bz-vz*Bx
				endcase
				Ey = 0.001*Ey ;; in mV/m
			endif else begin
				print, 'FLUX_CAL: No data!'
				t = !values.f_nan
			endelse
		endelse ;; else of flux type
		;; calculate the flux transport
		if finite(t(0)) then begin
			dt = t(1:*)-t(0:n_elements(t)-2)
			bflux_out(i,j) = total(Ey(1:*)*dt)
			min_this = min(Ey, /nan)
			max_this = max(Ey, /nan)
			if abs(max_this) ge abs(min_this) then peak_this = max_this else peak_this = min_this
			peak_ey(i,j) = peak_this
		endif else bflux_out(i,j) = !values.f_nan
	endfor ;; for of j, events
endfor ;; for of i, flux type
return, bflux_out
end
