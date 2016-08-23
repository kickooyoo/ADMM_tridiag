img = imread('textbook.jpg');
maxval = 256;
minval = 0;

tmp = img(:, 1:3827, :);
[Nx, Ny, Ncolor] = size(tmp);

level = 0.6;
bw = im2bw(tmp, level);
stats = regionprops(bw, 'PixelIdxList', 'PixelList', 'Area');

for ii = 1:length(stats)
	areas(ii) = stats(ii).Area;
end

bg_ndx = find(areas == max(areas));
bg_mask = stats(bg_ndx).PixelIdxList;
tmp_wh = reshape(tmp, Nx*Ny, Ncolor);
tmp_wh(bg_mask, :) = min(maxval - (maxval - tmp_wh(bg_mask,:))/2, maxval);
tmp_wh = reshape(tmp_wh, size(tmp));

level = 0.4;
bw = im2bw(tmp, level);
tmp_blk = reshape(tmp_wh, Nx*Ny, Ncolor);
tmp_blk(~bw, :) = max(tmp_blk(~bw,:)/2, minval);
tmp_blk = reshape(tmp_blk, size(tmp));
img = tmp_blk;
save('textbook_contrast.mat', 'img')
imwrite(img, 'textbook_contrast', 'jpg')

figure; imagesc(tmp)
figure; imagesc(tmp_wh)
