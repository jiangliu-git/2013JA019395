pro quantity_loc_plot, x, y, circ_q = circ_q, color_q = color_q, xtitle = xtitle, ytitle = ytitle, circ_title = circ_title, color_title = color_title, circ_legend = circ_legend, loc_legend = loc_legend, scale_circ_x = scale_circ_x, scale_circ_y = scale_circ_y, linear_circ = linear_circ, color_range = color_range, title = title, xrange = xrange, xstyle = xstyle, yrange = yrange, ystyle = ystyle, isotropic = isotropic
;; make circle-color-location plot, similar to Fu et al., 2012(last) used.
;; x, y: specify location, titles are x/ytitle
;; circ_q: quantity expressed as circle, title is circ_title, if positive it's a circ; if negative it's a triangle
;; circ_legend: can be an array, many circles to express the values circles stand for
;; loc_legend: location of the legend, in the ratio of the axis
;; scale_circ: the maximum radius of the circle quantities (in ratio to the smaller of x or y range)
;; linear_circ: if not set, circ value is expressed by circ area, if set, expressed by circ radius
;; color_q: the quantity expressed by color, title is color_title
;; color_range: the range of color for this
;; others are standard plot keywords
;; Liu Jiang 2013-2-4

if ~ keyword_set(cq_range) then cq_range = [min(color_q, /nan), max(color_q, /nan)]
if ~ keyword_set(color_range) then color_range = [7, 254]
if ~ keyword_set(scale_circ_x) then scale_circ_x = 0.02
if ~ keyword_set(scale_circ_y) then scale_circ_y = scale_circ_x
if ~ keyword_set(loc_legend) then loc_legend = [0.7, 0.95]

;; unit circle
theta=findgen(20)/19*2*!pi
unit_circ_x = cos(theta)
unit_circ_y = sin(theta)
unit_tri_x = [-1,1,0,-1]*sqrt(!pi/sqrt(3.))
unit_tri_y = [-1/sqrt(3.),-1/sqrt(3.),2/sqrt(3.),-1/sqrt(3.)]*sqrt(!pi/sqrt(3.))

;; get the color scale for color_q
color_data = color_range(0)+(color_range(1)-color_range(0))*(color_q-cq_range(0))/(cq_range(1)-cq_range(0))
i_large = where(color_data gt color_range(1), j_large)
if j_large gt 0 then color_data(i_large) = color_range(1)
i_small = where(color_data lt color_range(0), j_small)
if j_small gt 0 then color_data(i_small) = color_range(0)

;; draw the plot
plot, x, y, xtitle = xtitle, ytitle = ytitle, xrange = xrange, xstyle = xstyle, yrange = yrange, ystyle = ystyle, color = 0, psym = 3, title = title, isotropic = isotropic

;; get the scale for circ_q
xp_range = !x.crange(1)-!x.crange(0)
yp_range = !y.crange(1)-!y.crange(0)
circ_max = max(abs(circ_q))
div_circ_x = scale_circ_x*xp_range
if ~keyword_set(isotropic) then div_circ_y = scale_circ_y*yp_range else div_circ_y = div_circ_x

for i = 0, n_elements(x)-1 do begin
	if circ_q(i) gt 0 then begin
		unit_shape_x = unit_circ_x
		unit_shape_y = unit_circ_y
	endif else begin
		unit_shape_x = unit_tri_x
		unit_shape_y = unit_tri_y
	endelse
	if keyword_set(linear_circ) then factor = abs(circ_q(i)/circ_max) else factor = sqrt(abs(circ_q(i)/circ_max))
	oplot, unit_shape_x*factor*div_circ_x+x(i), unit_shape_y*factor*div_circ_y+y(i), color = color_data(i)
endfor

;;; draw legend
x_legend = !x.crange(0)+loc_legend(0)*xp_range
y_legend = !y.crange(0)+loc_legend(1)*yp_range
inc = 0.1*xp_range
for i = 0, n_elements(circ_legend)-1 do begin
	if circ_legend(i) gt 0 then begin
		unit_shape_x = unit_circ_x
		unit_shape_y = unit_circ_y
	endif else begin
		unit_shape_x = unit_tri_x
		unit_shape_y = unit_tri_y
	endelse
	if keyword_set(linear_circ) then factor = abs(circ_legend(i)/circ_max) else factor = sqrt(abs(circ_legend(i)/circ_max))
	polyfill, unit_shape_x*factor*div_circ_x+x_legend+i*inc, unit_shape_y*factor*div_circ_y+y_legend
	xyouts, x_legend+i*inc, y_legend-0.05*yp_range, strcompress(string(circ_legend(i), format='(f0)'), /remove)+'!Uo!n', alignment=0.5
endfor
xyouts, x_legend+0.05*xp_range, y_legend+0.05*yp_range, circ_title, /data

;;; draw the color bar
draw_color_scale, range = cq_range, brange = color_range, title = color_title

end
