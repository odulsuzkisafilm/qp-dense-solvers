function root = suite_init()

    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root, 'src'));
end
