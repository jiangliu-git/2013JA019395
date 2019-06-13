function dfb_duration, events, method = method, datafolder = datafolder
;;; get dfb durations, output is in seconds.
;;; events: | dbdt time| bz time| type (s,m)  | probe |
;;; method: specify the method to get it
;;;		'x': interpolate to several defined duration depending on X to get a duration.
;;;		'bz': get the average bz of t0+2min to t0+3min and select the first meet point after DF as the duration
;;; datafolder: the folder of the data needed
;;;	  if method is 'x': the folder of stat_pos
;;;	  if method is 'bz': the folder of fgs
if strcmp(method, 'x') then begin
	;;; the marks of x and duration
	;; for all events
	;x_marks = [-17.722568, -14.013473, -10.260011, -7.9005638] ;;; in RE
	;dur_marks = [97.32, 47.81, 36.444336, 36.148856] ;;; in seconds
	;; for 0809 events
	x_marks = [-17.722568, -14.013473, -10.338759, -7.9557382] ;;; in RE for all events
	dur_marks = [97.32, 47.81, 34.555398, 27.695121] ;;; in seconds
	;;; get x of the events
	pos_list = event_status_instant(events, 'pos', time_length=3, datafolder=datafolder)
	x = pos_list(2,*)/6371.
	durations = interpol(dur_marks, x_marks, x) 
endif

return, durations
end
