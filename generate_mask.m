function mask = generate_mask(im_type,loose,nx,ny);

% im_type describes whether I am working with the simulated brainweb image,
% or the real data. The position and ratios are hard-coded in...

switch im_type
    case {'sim','brainweb'}
        mask = gen_sim_mask(loose,nx,ny);
    case {'slice38'}
        mask = gen_slice38_mask(loose,nx,ny);
    case {'slice67'}
        mask = gen_slice67_mask(loose,nx,ny);
    otherwise
        error('invalid im_type');
end
        
end

function mask = gen_sim_mask(loose,nx,ny)
    assert(nx <= 258,'nx too large');
    assert(ny <= 258,'ny too large');
    fov = 25;
    ig = image_geom( 'nx', 258, 'ny', 258, 'fov', 25 );
    mask = logical(ig.circ( loose*fov / 2.4, loose*fov / 2.85, -5/fov, -3/fov));
    % crop to nx x ny
    extra_x = 258-nx;
    extra_y = 258-ny;
    if (mod(extra_x,2) == 0)
        mask = mask(1+extra_x/2:end-extra_x/2,:);
    else
        mask = mask(1+(extra_x+1)/2:end-(extra_x-1)/2,:);
    end
    if (mod(extra_x,2) == 0)
        mask = mask(:,1+extra_y/2:end-extra_y/2);
    else
        mask = mask(:,1+(extra_y+1)/2:end-(extra_y-1)/2);
    end
end

% to do: let crop to nx x ny
function mask = gen_slice38_mask(loose,nx,ny)
    fov = 25;
    ig = image_geom( 'nx', nx, 'ny', ny, 'fov', 25 );
    mask = logical(ig.circ( loose*fov/2.3, loose*fov/4, -4/fov, 2/fov));
end

% to do: let crop to nx x ny
function mask = gen_slice67_mask(loose,nx,ny)
    fov = 25;
    ig = image_geom( 'nx', nx, 'ny', ny, 'fov', 25 );
    mask = logical(ig.circ( loose*fov/2.5, loose*fov/4, -4/fov, 2/fov));
end
