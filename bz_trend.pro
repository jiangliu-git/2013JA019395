function bz_trend, events, method = method, datafolder = datafolder
;;; get statistical peak bzs, output is in seconds.
;;; events (bbf events): | start time| end time| type (s,m)  | probe |
;;; method: specify the method to get it
;;;		'x': interpolate to several defined Bz depending on X to get a Bz.
;;; datafolder: the folder of the data needed
;;;	  if method is 'x': the folder of stat_pos
if strcmp(method, 'x') then begin
	;;; the marks of x and bz
	;; for 0809 events
	x_marks = [-20.9, -16.5, -14.1, -12.05, -10.3, -8.75, -7.] ;;; in RE for all events
	bz_marks = [8.31028, 8.27705, 7.98683, 8.61190, 8.91469, 9.03745, 9.54003] ;;; in nT
	;;; get x of the events
	pos_list = event_status_instant(events, 'pos', time_length=3, datafolder=datafolder)
	x = pos_list(2,*)/6371.
	bz_out = interpol(bz_marks, x_marks, x) 
endif
return, bz_out
end
