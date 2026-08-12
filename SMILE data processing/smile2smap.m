clc; clear; close all;
smile_path = '..\Template data for tilt SMILE\';
smap_path = '.\smap_sum_sml.mat';    % Change the file path
load([smile_path 'frame_all.mat']);

% 1: xy
% 2: 6N
% 3: init
mode = 2;

if mode == 1
    load([smile_path, 'result_xy_all']);
    load(smap_path);
    
    loc_num = length(frame_all);
    pre_num = length(saveloc.loc.frame);
    
    field_name = fieldnames(saveloc.loc);
    for ii = 1:length(field_name)
        temp = saveloc.loc.(field_name{ii});
        saveloc.loc.(field_name{ii}) = [];
    end
    
    saveloc.loc.frame = frame_all;
    saveloc.loc.xnm = result_xy_all(:, 5);
    saveloc.loc.ynm = result_xy_all(:, 6);
    saveloc.loc.znm = result_xy_all(:, 7);
    saveloc.loc.phot = result_xy_all(:, 1) + result_xy_all(:, 2);
    saveloc.loc.channel = zeros(length(frame_all), 1);
    saveloc.loc.filenumber = ones(length(frame_all), 1);
    
    save([smile_path, 'smile_xy_smap_sml.mat'], "saveloc", "fileformat", "parameters");

elseif mode == 2
    load([smile_path, 'result_6N_all']);
    load(smap_path);
    
    loc_num = length(frame_all);
    pre_num = length(saveloc.loc.frame);
    
    field_name = fieldnames(saveloc.loc);
    for ii = 1:length(field_name)
        temp = saveloc.loc.(field_name{ii});
        saveloc.loc.(field_name{ii}) = [];
    end
    
    saveloc.loc.frame = frame_all;
    saveloc.loc.xnm = result_6N_all(:, 7);
    saveloc.loc.ynm = result_6N_all(:, 8);
    saveloc.loc.znm = result_6N_all(:, 9);
    saveloc.loc.channel = zeros(length(frame_all), 1);
    saveloc.loc.filenumber = ones(length(frame_all), 1);
    saveloc.loc.phot = sum(result_6N_all(:, 1:6), 2);
    
    save([smile_path, 'smile_6N_smap_sml.mat'], "saveloc", "fileformat", "parameters");

elseif mode == 3
    load([smile_path, 'result_init_all']);
    load(smap_path);
    
    loc_num = length(frame_all);
    pre_num = length(saveloc.loc.frame);
    
    field_name = fieldnames(saveloc.loc);
    for ii = 1:length(field_name)
        temp = saveloc.loc.(field_name{ii});
        saveloc.loc.(field_name{ii}) = [];
    end
    
    saveloc.loc.frame = frame_all;
    saveloc.loc.xnm = result_init_all(:, 2);
    saveloc.loc.ynm = result_init_all(:, 3);
    saveloc.loc.znm = result_init_all(:, 4);
    saveloc.loc.bg = result_init_all(:, 5);
    saveloc.loc.channel = zeros(length(frame_all), 1);
    saveloc.loc.filenumber = ones(length(frame_all), 1);
    saveloc.loc.phot = result_init_all(:, 1);
    
    save([smile_path, 'smile_init_smap_sml.mat'], "saveloc", "fileformat", "parameters");
    

end

fprintf("convert done\n");



