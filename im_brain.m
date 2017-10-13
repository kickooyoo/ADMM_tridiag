function im_brain(img, varargin)
arg.clim = [];
arg.orient = 'axial';
arg.abs = true;
arg = vararg_pair(arg, varargin);

switch arg.orient 
        case 'axial'
                perm = [2 1 3];
        case 'coronal'
                perm = [1 2 3];
                img = flipdim(img, 2);
        case 'sagittal'
                perm = [1 2 3];
                img = flipdim(img, 2);
        case 'sim'
                perm = [1 2 3];
        otherwise
                display(sprintf('%s is unknown orientation for im_brain()'));
                keyboard
end

        
if arg.abs 
        img = abs(img);
end
if isempty(arg.clim)
        im(permute(img, perm))
else
        im(permute(img, perm), arg.clim)
end

switch arg.orient
        case 'axial'
                set(gca, 'DataAspectRatio', [1 1.35 1])
        case 'coronal'
                set(gca, 'DataAspectRatio', [1 1.35 1])
        case 'sagittal'
                % no stretching needed
        case 'sim'
                % no stretching needed
        otherwise
                display(sprintf('%s is unknown orientation for im_brain()'));
                keyboard
end