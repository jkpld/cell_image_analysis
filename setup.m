function setup
% SETUP Add required files to the path and try to compile .c functions for
% increased speed.
%
% setup()

% James Kapaldo

% Check dependancies -----------------------------------------------------
if verLessThan('matlab','9.3')
    error('cell_image_analysis:setup','The cell_image_analysis code requires at least Matlab 2017b. (Dependence on function str2sym.)')
end
verImg = ver('images');
if isempty(verImg)
    error('cell_image_analysis:setup','The Image Processing Toolbox is required but not found.')
end
verStats = ver('stats');
if isempty(verStats)
    error('cell_image_analysis:setup','The Statistics and Machine Learning Toolbox is required but not found.')
end

fprintf('Setting up...\n')

% Get the path to the current folder -------------------------------------
fileLocation = mfilename('fullpath');
path = fileparts(fileLocation);

fprintf('...Found path\n')

% Compile the needed C files into .mex files -----------------------------
have_compiler = true;
try 
    evalc('mex(''-setup'',''c'')');
catch % ME
    have_compiler = false;
    warning('cell_image_analysis:setup','There is no supported C compiler installed.\nThe code will still run; however, it could be slower for 2D data without the compiled mex functions.')
%     rethrow(ME)
end

if have_compiler
    try

        mex(fullfile(path,'utilities','interp2mex.c'),'-outdir',fullfile(path,'utilities'),'-silent')
        fprintf('...Compiled interp2mex.c\n')
        
        mex(fullfile(path,'utilities','interp2mex_wExpand.c'),'-outdir',fullfile(path,'utilities'),'-silent')
        fprintf('...Compiled interp2mex_wExpand.c\n')

        mex(fullfile(path,'utilities','nakeinterp1.c'),'-outdir',fullfile(path,'utilities'),'-silent')
        fprintf('...Compiled nakeinterp1.c\n')
    catch % ME
    %     rethrow(ME)
        warning('cell_image_analysis:setup','There was an error compiling the required C functions ''interp2mex.c'', ''interp2mex_wExpand.c'' and ''nakeinterp1.c''. Make sure that the function ''mex'' is coorectly setup to compile C code.\nThe code will still run; however, it could be slower without the compiled mex functions.')
    end
end

% Copy ca_tifflib into the utilities folder --------------------------------

% if ~any(strncmp('salr_pdistmex',names,11))
    folder = fullfile(matlabroot,'toolbox','matlab','imagesci','private');
    d = dir(folder);
    names = {d.name};
    nameIdx = strncmp('tifflib',names,7);

    if ~any(nameIdx)
        error('cell_image_analysis:setup','File ''tifflib'' was not found in\n %s\nThis file is required. Try manually locating it; if found copy into the utilities folder and rename it to ''ca_tifflib.%s''.',folder,mexext)
    end

    copyfile(fullfile(folder,names{nameIdx}),fullfile(path,'utilities',['ca_' names{nameIdx}]))
% end
fprintf('...Added tifflib\n')

% K:\MATLAB\R2017a\toolbox\matlab\imagesci\private
% Add subfolders of current location to path -----------------------------
% (but do not include any .git repositories
pths = split(string(genpath(path)),';');
toIgnore = {'.git','docs','@'};
pths = pths(~pths.contains(toIgnore)).join(';');
addpath(pths.char());
fprintf('...Added subfolders to path\n')

fprintf('Setup finished!\n')
    
end
