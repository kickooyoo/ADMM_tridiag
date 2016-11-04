clear;
tridiag_setup;
orient = 'sagittal'; 
exp_ADMM_ALP2;
display('done with sagittal slice 69')

clear;
tridiag_setup;
orient = 'sim';
exp_ADMM_ALP2;
display('done with sim');

send_mai_text('done with all tests on ir72')
