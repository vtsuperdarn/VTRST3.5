pro	rad_fit_calc_velocity, force=force

common rad_data_blk
common rt_data_blk

; Find data index
data_index = rad_fit_get_data_index()

; Find run date and time
caldat, (*rad_fit_info[data_index]).sjul, mm1, dd1, yy1, hh1, mn1
caldat, (*rad_fit_info[data_index]).sjul, mm2, dd2, yy2, hh2, mn2
date = yy1*10000L + mm1*100L + dd1
time = [hh1*100L + mn1, hh2*100L + mn2]

; run ray tracing
rt_run, date, (*rad_fit_info[data_index]).code, time=time, force=force

; Apply correction to velocity (mask)
njuls = n_elements((*rad_fit_data[data_index]).juls)
novelinds = where((*rad_fit_data[data_index]).velocity eq 10000.)
for it=0,n_elements(rt_data.juls[*,0])-2 do begin
	for ib=0,n_elements(rt_data.beam[0,*])-1 do begin
		for j=0,rt_info.ngates-2 do begin
			julinds = where((*rad_fit_data[data_index]).juls ge rt_data.juls[it,0] and $
							(*rad_fit_data[data_index]).juls lt rt_data.juls[it+1,0] and $
							(*rad_fit_data[data_index]).beam eq rt_data.beam[it,ib], cc)
			; If there is a correction to apply, do it
			if cc gt 0 and rt_data.nr[it,ib,j] gt 0 then begin
				(*rad_fit_data[data_index]).velocity[julinds + j*njuls] = (*rad_fit_data[data_index]).velocity[julinds + j*njuls] * 1./rt_data.nr[it,ib,j]
			endif
		endfor
	endfor
endfor
; Enforce no velocity value (the correction overwrote this, which is why it needs to be enforced here)
(*rad_fit_data[data_index]).velocity[novelinds] = 10000.

rad_fit_data2ascii,date,rt_info.name



end