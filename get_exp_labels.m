function [slice_str, curr_folder] = get_exp_labels(orient, slice, reduction)
%function [slice_str, curr_folder] = get_exp_labels(orient, slice, reduction)
% for file saving and loading

base = './reviv/';
switch orient
case 'sim'
	slice_str = 'sim';
	curr_folder = sprintf('%ssim/sim_r%d', base, reduction);
case 'axial'
	slice_str = sprintf('slice%d', slice);
	curr_folder = sprintf('%saxial/slice%d/axial_slice%d_r%d', base, slice, slice, reduction);
case 'sagittal'
	slice_str = sprintf('sagittal%d', slice);
	curr_folder = sprintf('%ssagittal/slice%d/sagittal_slice%d_r%d', base, slice, slice, reduction);
case 'coronal'
	curr_folder = sprintf('%scoronal/slice%d/coronal_slice%d_r%d', base, slice, slice, reduction);
	slice_str = sprintf('coronal%d', slice);

otherwise
	display(sprintf('unknwon orientation %s', orient))
	keyboard
end
if ~exist(curr_folder, 'dir')
	fsep = strfind(curr_folder ,'/');
	parent = curr_folder(1:fsep(end)-1);
	child = curr_folder(fsep(end)- (-1):end);
	display(sprintf('mkdir %s in %s?', child, parent));
	keyboard;
	[status, msg] = mkdir(parent, child);
	if ~status || ~isempty(msg)
		display('msg');
		display('directory not created');
		keyboard;
	end
end

