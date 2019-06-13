;;; declare things
;; the folder to put pics
dataroot_folder = '/Work/data'
;plasma_folder = '/plasma_old'
plasma_folder = '/plasma_new'

if keyword_set(computer) then begin
	fgs_folder=computer + '/Work/data/fgs'
	fgh_folder=computer + '/Work/data/fgh'
	fgl_folder=computer + '/Work/data/fgl'
	pos_folder=computer + '/Work/data/pos'
	efs_folder=computer + '/Work/data/efs_dsl'
	kyoto_ae_folder=computer + '/Work/data/ground/kyoto_ae'
	pseudo_ae_folder=computer + '/Work/data/ground/pseudo_ae'
	pseudo_al_folder=computer + '/Work/data/ground/pseudo_al'
	omni_b_folder=computer + '/Work/data/sw/omni_b_gsm'
	omni_Pdyn_folder=computer + '/Work/data/sw/omni_Pdyn'
	omni_vxb_folder=computer + '/Work/data/sw/omni_vxb'
	Pall_folder=computer + dataroot_folder + plasma_folder + '/Pall'
	Pth_folder=computer + dataroot_folder + plasma_folder + '/Pth'
	beta_folder=computer + dataroot_folder + plasma_folder + '/beta'
	viperp_folder = computer + dataroot_folder + plasma_folder + '/viperp'
	vixy_folder = computer + dataroot_folder + plasma_folder + '/vixy_infer'
	vi_folder = computer + dataroot_folder + plasma_folder + '/vi'
	vixb_folder = computer + dataroot_folder + plasma_folder + '/vixb'
	vexb_folder = computer + dataroot_folder + plasma_folder + '/vexb'
	vexb_dsl_folder = computer + dataroot_folder + plasma_folder + '/vexb_dsl'
	ni_folder = computer + dataroot_folder + plasma_folder + '/ni'
endif

listfolder='../../lists'
pic_folder = '../../pics/publication/bbf_flux'
pic_folder_events = '../../pics_test/events'
save_folder = 'variables'

RE = 6371.
mu0 = 4*!pi*1e-7
nTesla2_to_nPa=0.01/25.132741

;;;line thickness for plot
l_thick = 2.4

;;; constants and signs
Re = 6371.
mu0 = 4*!pi*1e-7
nTesla2_to_nPa = 0.01/25.132741
perp_sign = '!9'+string("136B)+'!X'
cross = '!9'+string("264B)+'!X'
theta_letter = '!9'+string("161B)+'!X'
delta_letter = '!9'+string("144B)+'!X'
phi_letter = '!9'+string("146B)+'!X'
intg_sign = '!9'+string("362B)+'!X'
jiao_l = '!9'+string("341B)+'!X'
jiao_r = '!9'+string("361B)+'!X'

;;; turn off the timestamp of tplot
time_stamp,/off
