function dfb_criteria, pre_b_data, window_b_data, btotal = btotal, ignore = ignore, bz = bz, dbz = dbz, m_bz = m_bz, m_dbz = m_dbz, mm_dbz = mm_dbz
;; the criteria are specified in this procedure
;; criteria about dBz/dt is in dfb_select_bottom2.pro
;; pre_data: the data previous to the window
;; window_data: the data inside the window
;; btotal: use btotal instead of Bz
;; ignore: if set, disregard all criteria
;; format of data: [[Bx],[By],[Bz],[time_double]]
;; all behind dbz are thresholds
;; Jiang Liu, 3/23/2011

if keyword_set(ignore) then return, 1

;;; thresholds, defaults are for DFCS paper
; average values
if ~ keyword_set(bz) then bz = 0. ; in nT
if ~ keyword_set(thetab) then thetab = -1000. ; in degrees
if ~ keyword_set(dbz) then dbz = 0. ; in nT
if ~ keyword_set(dthetab) then dthetab = -1000. ; in degrees
; peak values / peak to pre average values
if ~ keyword_set(m_bz) then m_bz = 5. ; in nT
if ~ keyword_set(m_thetab) then m_thetab = -1000. ; in degrees
if ~ keyword_set(m_dbz) then m_dbz = 0. ; in nT
if ~ keyword_set(m_dthetab) then m_dthetab = -1000. ; in degrees
; peak to peak values
;mm_dbz = 3. ;; originally used
if ~ keyword_set(mm_dbz) then mm_dbz = 5. ;; used for DFCS statistic paper
;mm_dbz = -100. ;; disable this criteria, for ignoring criterion 2

if keyword_set(btotal) then begin
	comp_use = sqrt(total(window_b_data^2, 2))
	pre_comp_use = sqrt(total(pre_b_data(*,2)^2, 2))
endif else begin
	comp_use = window_b_data(*,2)
	pre_comp_use = pre_b_data(*,2)
endelse

; begin to work
pre_thetab_data = elevation_angle(pre_b_data)
window_thetab_data = elevation_angle(window_b_data)
; averages
pre_bz = mean(pre_comp_use)
pre_thetab = mean(pre_thetab_data)

window_bz = mean(comp_use)
window_thetab = mean(window_thetab_data)

; peaks
m_window_bz = max(comp_use)
m_pre_bz = max(pre_comp_use)
m_window_thetab = max(window_thetab_data)

; criteria
c_bz = (window_bz gt bz)
c_thetab = (window_thetab gt thetab)
c_dbz = (window_bz-pre_bz gt dbz)
c_dthetab = (window_thetab-pre_thetab gt dthetab)
c_maxb = (m_window_bz gt m_bz)
c_maxthetab = (m_window_thetab gt m_thetab)
c_maxb2pre = (m_window_bz-pre_bz gt m_dbz)
c_maxthetab2pre = (m_window_thetab-pre_thetab gt m_dthetab)
c_maxb2maxb = (m_window_bz-m_pre_bz gt mm_dbz)

if c_bz and c_thetab and c_dbz and c_dthetab and c_maxb and c_maxthetab and c_maxb2pre and c_maxthetab2pre and c_maxb2maxb then yes = 1 else yes = 0
return, yes
end
