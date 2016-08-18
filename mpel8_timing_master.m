%clear;
%tridiag_setup;
%orient = 'sagittal'; truncate = 2; niters = 1000;
%exp_ADMM_ALP2;
%display('done with slice 38')

clear;
tridiag_setup;
orient = 'coronal'; truncate = 2; niters = 1000;
exp_ADMM_ALP2;
display('done with slice 38')

send_mai_text('done with sagittal tests on iv1')
return

clear;
tridiag_setup;
orient = 'axial'; slice = 38; truncate = 2;
exp_ADMM_ALP2;
display('done with slice 38')

clear;
tridiag_setup;
orient = 'axial'; slice = 90; truncate = 2;
exp_ADMM_ALP2;
display('done with slice 90');

%clear;
%tridiag_setup;
%orient = 'sim';
%exp_ADMM_ALP2;
%display('done with sim');

send_mai_text('done with all tests on mpel8')
