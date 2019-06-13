pro main_test
;;; test things
;thm_init
;
;;; the event
;computer = 'I:'
;;computer = '/home/jliu'
;
;;; folders
;@declare
;
;listfolder='../../lists'
;picfolder = '../../pics_test'
;picfolder_events = '../../pics_test/events'
;
;listname = 'bbf_ps_list_original2.txt'
;
;events = load_list(listname, folder = listfolder)
;
;i_event = 1
;
;bbf_trange = time_double(events(0:1, i_event))
;sc = events(3, i_event)
;
;fb_tranges = tr_flow_burst(bbf_trange, probe = sc, v_folder = vi_folder, secs_fb = secs_fb, percent_t_fb = percent_t_fb)
;
;load_bin_data, trange = bbf_trange+[-120., 120.], probe = sc, /tclip, datatype = 'vi', datafolder = vi_folder
;get_data, 'th'+sc+'_ptix_velocity_gsm_tclip', data = data
;vx = data.y(*,0)
;vy = data.y(*,1)
;vz = data.y(*,2)
;vttl = sqrt(vx^2+vy^2+vz^2)
;store_data, 'th'+sc+'_vttl', data = {x:data.x, y:vttl}
;tplot, 'th'+sc+'_vttl', trange = bbf_trange+[-120., 120.]
;timebar, 400., varname = 'th'+sc+'_vttl', /databar
;timebar_mass, fb_tranges(0,*), varname = 'th'+sc+'_vttl'
;timebar_mass, fb_tranges(1,*), varname = 'th'+sc+'_vttl', line = 1

;;;;;; test ranges shared and ranges compressed
limit = [12., 15.]

ranges_in1 = [[1., 2.], [7., 8.], [9., 13.]]
ranges_in2 = [[1.5, 2.5], [7.1, 7.9], [8.9, 13.], [14, 15]]

ranges_out = ranges_inter(ranges_in1, ranges_in2)
print, ranges_out
end
