%% DEMO

% This script was created by Umut Güçlü and last modified on 2 May 2014 at
% 21:15
% Please visit http://www.ccnlab.net/ for more information

%% Cleanup

close all;
clear all;
clc;

%% Seed random number generator

rng(1);

%% Load stimuli

stimuli = load('Stimuli.mat');

train_im = permute(double(stimuli.stimTrn), [2 3 1]);
val_im   = permute(double(stimuli.stimVal), [2 3 1]);

%% Train feature model

var.im_patch_num  = 50000;
var.im_patch_off  = [23 23];
var.im_patch_size = [32 32];
var.it_num        = 2000;
var.lay_size      = [25 25];
var.neigh_size    = [5 5];

fea_mo = train_feature_model(train_im, var);

%% Simulate feature model

train_im_fea = simulate_feature_model(fea_mo, train_im);
val_im_fea   = simulate_feature_model(fea_mo, val_im);

%% Load responses

var.ROI = 1;
var.sub = 1;

estimated_responses = load('EstimatedResponses.mat');

if var.sub == 1
    
    train_vox = estimated_responses.dataTrnS1(var.ROI == estimated_responses.roiS1 & all(isfinite(estimated_responses.dataTrnS1), 2) & all(isfinite(estimated_responses.dataValS1), 2), :)';
    val_vox   = estimated_responses.dataValS1(var.ROI == estimated_responses.roiS1 & all(isfinite(estimated_responses.dataTrnS1), 2) & all(isfinite(estimated_responses.dataValS1), 2), :)';
    
elseif var.sub == 2
    
    train_vox = estimated_responses.dataTrnS2(var.ROI == estimated_responses.roiS2 & all(isfinite(estimated_responses.dataTrnS2), 2) & all(isfinite(estimated_responses.dataValS2), 2), :)';
    val_vox   = estimated_responses.dataValS2(var.ROI == estimated_responses.roiS2 & all(isfinite(estimated_responses.dataTrnS2), 2) & all(isfinite(estimated_responses.dataValS2), 2), :)';
    
end

%% Train voxel model

var.df_num = 100;
var.df_tol = 1e-3;

vox_mo = train_voxel_model(train_im_fea, var, train_vox);

%% Simulate voxel model

val_vox_hat = simulate_voxel_model(val_im_fea, vox_mo);

rho = sort(diag(corr(val_vox, val_vox_hat)), 'descend');

figure(1);

semilogx(rho);

title('encoding performance');
xlabel('voxel');
ylabel('rho');

%% Simulate decoding model

vox_num = 500;

im_hat = simulate_decoding_model(fea_mo, val_im, val_vox, vox_mo, vox_num);

val_im_num = size(val_im, 3);

ac = round(100 * sum(im_hat == 1 : val_im_num) / val_im_num);

disp(['decoding performance: ' num2str(ac) '%']);