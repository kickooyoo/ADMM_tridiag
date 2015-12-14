function S = construct_sensefat(sense_maps)

% input sense_maps has dims [nx ny nc]

% maybe wrong dims??
dims = [size(sense_maps,1) size(sense_maps,2)];
%dims = [size(sense_maps,2) size(sense_maps,1)];

nc = size(sense_maps,3);

for ii = 1:nc
    curr_smap = sense_maps(:,:,ii);
    %curr_smap = sense_maps(:,:,ii).';
    S{ii} = Gdiag(curr_smap(:));
end
S = vertcat(S{:});
