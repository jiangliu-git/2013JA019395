pro main_sublists_bbf
;; create sublists for use
computer = 'I:'
;computer = '/home/jliu'

@declare

;;;;;;;;; 1. make another original list with more strict sc location requirements, not in the shadow, correct for sun pulse, and in fast survey ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;; make sure the loaded list is 'bbf_list_original.txt' ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; get the sc positions of the events
;;listname = 'bbf_list_original.txt'
;listname = 'bbf_ps_list_original.txt'
;events = load_list(listname, folder = listfolder)
;pos_list = event_status_instant(events, 'pos', vtrange = time_string(events(0:1,*)), datafolder=pos_folder)
;pos_list(2:*,*) = pos_list(2:*,*)/6371. ; ...|x|y|z|
;roi_list =event_status_instant(events, 'roi', vtrange = time_string(events(0:1,*)), datafolder=roi_folder)
;;; must be fast survey
;mode_list =event_status_instant(events, 'mode', vtrange = time_string(events(0:1,*)), datafolder=mode_folder)
;
;;; create a new poslist with ...|x|rho|sun shadow|moon shadow|
;pos_shadow_new = [pos_list(0:2,*), sqrt(pos_list(3,*)^2+pos_list(4,*)^2), fix(roi_list(2,*), type = 3) and 1, fix(roi_list(2,*), type = 3) and 2, mode_list(2,*)]
;;; x range: -6 to -30, rho<12
;good_list = list_split(pos_shadow_new, criteria=[[-30., 12., 1, 1, 0.5], [-6., 0, 0 ,0, 0]], method = ['between','lt','lt','lt','gt'])
;
;events_newpos = events(*, good_list(0,*))
;events_good = events_newpos
;;output_txt, events_good, filename = listfolder+'/bbf_list_original2.txt'
;output_txt, events_good, filename = listfolder+'/bbf_ps_list_original2.txt'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;; 2. make a list with quiet pre-condition ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
listname = 'bbf_ps_list_original2.txt'
events = load_list(listname, folder = listfolder)
vam_peak_list = event_status_instant(events, 'vi', pre_time = 3.5, time_length=3, datafolder=vi_folder, /peak)
vam_peak = vam_peak_list(5,*)
i_good = where(vam_peak lt 50)
events_good = events(*, i_good)
output_txt, events_good, filename = listfolder+'/bbf_ps_list_lead.txt'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
