function mask = generate_mask(orient, slice, Nx, Ny, varargin);
%function mask = generate_mask(orient, slice, Nx, Ny, varargin)
% varragin:
% 	loose (scalar), if > 1, makes mask roomier
%		default: 1

arg.loose = 1;
arg = vararg_pair(arg, varargin);

switch orient
case 'axial'
	mask = gen_axial_mask(slice, arg.loose, Nx, Ny);
case 'sim'
	mask = gen_sim_mask(arg.loose, Nx, Ny);
case 'coronal'
	mask = cat(1, zeros(5, 128), ones(144-13,128), zeros(8,128));
	display('warning: poor man coronal mask');
otherwise
	display(sprintf('no mask determined for %s orientation, using all true', orient));
	mask = true(Nx, Ny);
end
        
end

function mask = gen_sim_mask(loose, Nx, Ny)
    assert(Nx <= 258,'Nx too large');
    assert(Ny <= 258,'Ny too large');
    fov = 25;
    ig = image_geom( 'Nx', 258, 'Ny', 258, 'fov', 25 );
    mask = logical(ig.circ( loose*fov / 2.4, loose*fov / 2.85, -5/fov, -3/fov));
    % crop to Nx x Ny
    extra_x = 258-Nx;
    extra_y = 258-Ny;
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

function mask = gen_axial_mask(slice, loose, Nx, Ny)
fov = 25;
ig = image_geom('nx', Nx, 'ny', Ny, 'fov', fov);

% to do: let crop to Nx x Ny
switch slice
case 38
	mask = logical(ig.circ( loose*fov/2.3, loose*fov/4, -4/fov, 2/fov));

case 67
	mask = logical(ig.circ( loose*fov/2.5, loose*fov/4, -4/fov, 2/fov));

case 90
	mask = logical(ig.circ( loose*fov/2.7, loose*fov/4.2, 0/fov, 4/fov));

otherwise
	display(sprintf('no mask determined for axial slice %d, using all true', slice));
	mask = true(Nx, Ny);
end

end

