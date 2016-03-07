function trunc_sense_map = truncate_sense_maps(sense_maps, mask)
%function trunc_sense_map = truncate_sense_maps(sense_maps, mask)

Nc = size(sense_maps, 3);
mc_mask = repmat(mask, [1 1 Nc]);
assert(numel(sense_maps) == numel(mc_mask), 'size of mask does not match sense_maps');
trunc_sense_map = zeros(size(sense_maps)); 
for coil_ndx = 1:Nc
	curr_map = sense_maps(:,:,coil_ndx);
	max_in_mask = max(col(abs(curr_map(mask))));
	truncate_mask = abs(curr_map) > max_in_mask;
	trunc_sense_map(:,:,coil_ndx) = curr_map.*(~truncate_mask) + max_in_mask*exp(1i*angle(curr_map)).*truncate_mask;
end
