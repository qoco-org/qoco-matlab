% Copyright (c) 2025 Govind M. Chari
% Copyright (c) 2017 Bartolomeo Stellato
%
% This source code is licensed under the BSD 3-Clause License

function make_qoco(varargin)
% Matlab MEX makefile for QOCO.
%

%    MAKE_QOCO(VARARGIN) is a make file for QOCO solver. It
%    builds QOCO and its components from source.
%
%    WHAT is the last element of VARARGIN and cell array of strings,
%    with the following options:
%
%    {}, '' (empty string) or 'all': build all components and link.
%
%    'qoco': builds the QOCO solver using CMake
%
%    'qoco_mex': builds the QOCO mex interface and links it to the QOCO
%    library
%
%    VARARGIN{1:NARGIN-1} specifies the optional flags passed to the compiler
%
%    Additional commands:
%
%    'clean': delete all object files (.o and .obj)
%    'purge' : same as above, and also delete the mex files.


if( nargin == 0 )
    what = {'all'};
    verbose = false;
elseif ( nargin == 1 && ismember('-verbose', varargin) )
    what = {'all'};
    verbose = true;
else
    what = varargin{nargin};
    if(isempty(strfind(what, 'all'))        && ...
            isempty(strfind(what, 'qoco'))      && ...
            isempty(strfind(what, 'qoco_mex'))  && ...
            isempty(strfind(what, 'clean'))     && ...
            isempty(strfind(what, 'purge')))
        fprintf('No rule to make target "%s", exiting.\n', what);
    end
    
    verbose = ismember('-verbose', varargin);
end

%% Try to unlock any pre-existing version of qoco_mex

% this prevents compile errors if a user builds, runs qoco
% and then tries to recompile
if(mislocked('qoco_mex'))
    munlock('qoco_mex');
end



%% Basic compile commands

% Get make and mex commands
make_cmd = 'cmake --build .';
mex_cmd = sprintf('mex -O -silent');
mex_libs = '';


% Add arguments to cmake and mex compiler
cmake_args = '-DMATLAB=ON';
mexoptflags = '-DMATLAB';

% % Add specific generators for windows linux or mac
% if (ispc)
%     cmake_args = sprintf('%s %s', cmake_args, '-G "MinGW Makefiles"');
% else
%     cmake_args = sprintf('%s %s', cmake_args, '-G "Unix Makefiles"');
% end

% Pass Matlab root to cmake
Matlab_ROOT = strrep(matlabroot, '\', '/');
cmake_args = sprintf('%s %s%s%s', cmake_args, ...
    '-DMatlab_ROOT_DIR="', Matlab_ROOT, '"');

% Add parameters options to mex and cmake
% CTRLC
% if (ispc)
%    ut = fullfile(matlabroot, 'extern', 'lib', computer('arch'), ...
%                  'mingw64', 'libut.lib');
%    mex_libs = sprintf('%s "%s"', mex_libs, ut);
% else
%    mex_libs = sprintf('%s %s', mex_libs, '-lut');
% end
% Shared library loading
if (isunix && ~ismac)
    mex_libs = sprintf('%s %s', mex_libs, '-ldl');
end

% Set OS flag for timer.
if (isunix && ~ismac)
    mex_libs = sprintf('%s %s', mex_libs, '-DIS_LINUX');
end

if (ismac)
    mex_libs = sprintf('%s %s', mex_libs, '-DIS_MACOS');
end


% Add large arrays support if computer is 64 bit and a pre-2018 version
% Release R2018a corresponds to Matlab version 9.4
if (~isempty(strfind(computer, '64')) && verLessThan('matlab', '9.4'))
    mexoptflags = sprintf('%s %s', mexoptflags, '-largeArrayDims');
end

%Force Matlab to respect old-style usage of mxGetPr in releases after 2018a,
%which use interleaved complex data.   Note that the -R2017b flag is badly
%named since it indicates that non-interleaved complex data model is being used;
%it is not really specific to the release year
if ~verLessThan('matlab', '9.4')
    mexoptflags = sprintf('%s %s', mexoptflags, '-R2017b');
end


% Set optimizer flag
if (~ispc)
    mexoptflags = sprintf('%s %s', mexoptflags, 'COPTIMFLAGS=''-O3''');
end

% Set library extension
lib_ext = '.a';
lib_name = sprintf('libqocostatic%s', lib_ext);
qdldl_lib_name = 'libqdldl.a';


% Set qoco directory and qoco_build directory
current_dir = pwd;
[makefile_path,~,~] = fileparts(which('make_qoco.m'));
qoco_dir = fullfile(makefile_path, 'qoco');
qoco_build_dir = fullfile(qoco_dir, 'build');
qdldl_dir = fullfile(qoco_dir, 'lib', 'qdldl');
qdldl_types_dir = fullfile(qoco_build_dir, 'lib', 'qdldl');

% Include directory
inc_dir = [
    fullfile(sprintf(' -I%s', qoco_dir), 'include'), ...
    sprintf(' -I%s', qdldl_dir), ...
    fullfile(sprintf(' -I%s', qdldl_dir), 'include'), ...
    fullfile(sprintf(' -I%s', qdldl_types_dir), 'include'),
    ];


%% QOCO Solver
if( any(strcmpi(what,'qoco')) || any(strcmpi(what,'all')) )
    fprintf('Compiling QOCO solver...');
    
    % Create build directory and go inside
    if exist(qoco_build_dir, 'dir')
        rmdir(qoco_build_dir, 's');
    end
    mkdir(qoco_build_dir);
    cd(qoco_build_dir);
    
    % Extend path for CMake mac (via Homebrew)
    PATH = getenv('PATH');
    if ((ismac) && (isempty(strfind(PATH, 'opt/homebrew/bin'))))
        setenv('PATH', [PATH ':/opt/homebrew/bin']);
    end
    
    % Compile static library with CMake
    [status, output] = system(sprintf('%s %s ..', 'cmake', cmake_args));
    if(status)
        fprintf('\n');
        disp(output);
        error('Error configuring CMake environment');
    elseif(verbose)
        fprintf('\n');
        disp(output);
    end
    
    [status, output] = system(sprintf('%s %s', make_cmd, '--target qocostatic'));
    if (status)
        fprintf('\n');
        disp(output);
        error('Error compiling QOCO');
    elseif(verbose)
        fprintf('\n');
        disp(output);
    end
    
    
    % Change directory back to matlab interface
    cd(makefile_path);
    
    % Copy static library to current folder
    lib_origin = fullfile(qoco_build_dir, lib_name);
    copyfile(lib_origin, lib_name);
    
    % Copy libqdldl.a to current folder
    lib_origin = fullfile(qoco_build_dir, 'lib', 'qdldl', 'out', qdldl_lib_name);
    copyfile(lib_origin, qdldl_lib_name);
    
    fprintf('\t\t\t\t\t\t[done]\n');
    
end

%% qocomex
if( any(strcmpi(what,'qoco_mex')) || any(strcmpi(what,'all')) )
    % Compile interface
    fprintf('Compiling and linking qocomex...');
    
    % Compile command
    
    cmd = sprintf('%s %s %s %s %s qoco_mex.cpp %s', ...
        mex_cmd, mexoptflags, inc_dir, lib_name, qdldl_lib_name, mex_libs);
    
    % Compile
    eval(cmd);
    fprintf('\t\t\t\t\t[done]\n');
    
end


%% clean
if( any(strcmpi(what,'clean')) || any(strcmpi(what,'purge')) )
    fprintf('Cleaning mex files and library...');
    
    % Delete mex file
    mexfiles = dir(['*.', mexext]);
    for i = 1 : length(mexfiles)
        delete(mexfiles(i).name);
    end
    
    % Delete static libraries
    lib_full_path = fullfile(makefile_path, lib_name);
    if( exist(lib_full_path,'file') )
        delete(lib_full_path);
    end
    qdldl_lib_full_path = fullfile(makefile_path, qdldl_lib_name);
    if( exist(qdldl_lib_full_path,'file') )
        delete(qdldl_lib_full_path);
    end
    
    
    fprintf('\t\t\t[done]\n');
end


%% purge
if( any(strcmpi(what,'purge')) )
    fprintf('Cleaning QOCO build directory...');
    
    % Delete QOCO build directory
    if exist(qoco_build_dir, 'dir')
        rmdir(qoco_build_dir, 's');
    end
    
    fprintf('\t\t[done]\n');
end


%% Go back to the original directory
cd(current_dir);

end