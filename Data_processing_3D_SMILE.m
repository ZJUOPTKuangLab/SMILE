%% This is code is for processing the experiment data （3D-SMILE）



close all;
clear;
clc

addpath('shared/')
addpath('code/')
addpath('helpfunction/')
support_path = pwd ; %check support path
addpath([support_path '\Segmentation']);
addpath([support_path '\Helpers']);
single_maxpoint = 100000;   % The maximum candidates number for the GPU to process (to avoid exceeding the GPU memory)
% addpath('Visualization\')
addpath('bigtiff\')


%% 1. Read PSF, Generate cubic spline coefficients (Remember to revise)
load('Template data for 3D SMILE\psf\coeff_VISPR_NPC_3D_SMILE.mat')    % load cubic-spline parameters ========================================================================================
PSF = insitu_PSF;

for i = 1:size(PSF,3) % slice by slice normalization for energy conservation
    tmp = PSF(:,:,i);
    tmp = tmp/sum(tmp(:));
    PSF(:,:,i) = tmp;
end
coeff = Spline3D_interp(PSF);

x_pixelsize = 83.6;   
z_pixelsize = 10;

[x, ~, z , ~] = size(coeff);

% % for test only +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% figure;orthosliceViewer(PSF(:,:,1:10:end)); title('PSF');

% other parameters (wavelength, NA, phase interval)
lambda = 560;
NA = 1.45;
phase_difference = (120)/180*pi;

% Select a region of the rendering images for pattern estimation, no need to change
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------%
xhalf = (7500)/x_pixelsize;  % half FOV of xy dimension, unit: nm
xhalf_z =  (5500)/x_pixelsize;  % half size for z dimension, unit nm
yhalf = xhalf;

n = 1.3384;   % The refractive index of the sample medium 
swcycle = lambda/2/n;   % minimum period of the non-evanescent waves (unit: nm)
swcycle_pixel = swcycle/x_pixelsize;   % minimum period of the non-evanescent waves (unit: pixel)
k_pixel = 2*pi/swcycle_pixel;   % maximum pattern k value (unit: pixel^-1)
Tmin = swcycle;






%% 2. Read the raw image stack (batch 1 round 1 to determine the parameters)，and apply sliding window
% ===================================================================== Remember to revise these parameters ==================================================================================================================
total_batch = 12;
batch_frame = 30000;  % how many frames in each batch
round_size = 6000;  % for each batch, round_size frames are used to estimate local parameters
total_round = floor(batch_frame/round_size)  % number of round for each batch
file_title = 'Template data for 3D SMILE\raw\';
Last_batch = 12;
Last_round = 5;   % This means the reconstruction will end at the Last_round of Last_batch

save_path = 'Template data for 3D SMILE\'; % Remember to revise the save path ====================================================================================================================================================================================
% ==================================================================== END：Remember to revise these parameters ==================================================================================================================

filepath = [file_title, '1.tif'];    % ====================================================================================================================================================================================
tmp_image = loadtiff(filepath);
tmp_image = double(tmp_image);
tmp_image = tmp_image(:,:,1+0*round_size:1*round_size);   % round 1  1~6000
% tmp_image = tmp_image(:,:,1+1*round_size:2*round_size);   % round 2  6001~12000 
% tmp_image = tmp_image(:,:,1+2*round_size:3*round_size);   % round 3  12001~18000
% tmp_image = tmp_image(:,:,1+3*round_size:4*round_size);   % round 4  18001~24000
% tmp_image = tmp_image(:,:,1+4*round_size:5*round_size);   % round 5  24001~30000

% tmp_image = tmp_image(:,:,1:600);   % for test ============================================================================================================================================
  

% Pay attention to the positive directions of x and y: 
% the positive direction of x is along the rows, and the positive direction of y is along the columns.
% x1 with rolling
tmp1_rolling0 = tmp_image(:,:,1:6:end);     % 500 frames
tmp1_rolling1 = tmp_image(:,:,1+6:6:end);   % 499 frames
tmp1_rolling2 = tmp_image(:,:,1+6:6:end);   % 499 frames
tmp1_rolling3 = tmp_image(:,:,1+6:6:end);   % 499 frames
tmp1_rolling4 = tmp_image(:,:,1+6:6:end);   % 499 frames
tmp1_rolling5 = tmp_image(:,:,1+6:6:end);   % 499 frames
tmp1 = cat(3,tmp1_rolling0, tmp1_rolling1, tmp1_rolling2, tmp1_rolling3, tmp1_rolling4, tmp1_rolling5);
clear tmp1_rolling1 tmp1_rolling2 tmp1_rolling3 tmp1_rolling4 tmp1_rolling5

% x2 with rolling
tmp2_rolling0 = tmp_image(:,:,2:6:end);     % 500 frames
tmp2_rolling1 = tmp_image(:,:,2:6:end-6);
tmp2_rolling2 = tmp_image(:,:,2+6:6:end);
tmp2_rolling3 = tmp_image(:,:,2+6:6:end);
tmp2_rolling4 = tmp_image(:,:,2+6:6:end);
tmp2_rolling5 = tmp_image(:,:,2+6:6:end);
tmp2 = cat(3,tmp2_rolling0, tmp2_rolling1, tmp2_rolling2, tmp2_rolling3, tmp2_rolling4, tmp2_rolling5);
clear tmp2_rolling1 tmp2_rolling2 tmp2_rolling3 tmp2_rolling4 tmp2_rolling5

% x3 with rolling
tmp3_rolling0 = tmp_image(:,:,3:6:end);
tmp3_rolling1 = tmp_image(:,:,3:6:end-6);
tmp3_rolling2 = tmp_image(:,:,3:6:end-6);
tmp3_rolling3 = tmp_image(:,:,3+6:6:end);
tmp3_rolling4 = tmp_image(:,:,3+6:6:end);
tmp3_rolling5 = tmp_image(:,:,3+6:6:end);
tmp3 = cat(3,tmp3_rolling0, tmp3_rolling1, tmp3_rolling2, tmp3_rolling3, tmp3_rolling4, tmp3_rolling5);
clear tmp3_rolling1 tmp3_rolling2 tmp3_rolling3 tmp3_rolling4 tmp3_rolling5

% y1 with rolling
tmp4_rolling0 = tmp_image(:,:,4:6:end);
tmp4_rolling1 = tmp_image(:,:,4:6:end-6);
tmp4_rolling2 = tmp_image(:,:,4:6:end-6);
tmp4_rolling3 = tmp_image(:,:,4:6:end-6);
tmp4_rolling4 = tmp_image(:,:,4+6:6:end);
tmp4_rolling5 = tmp_image(:,:,4+6:6:end);
tmp4 = cat(3,tmp4_rolling0, tmp4_rolling1, tmp4_rolling2, tmp4_rolling3, tmp4_rolling4, tmp4_rolling5);
clear tmp4_rolling1 tmp4_rolling2 tmp4_rolling3 tmp4_rolling4 tmp4_rolling5

% y2 with rolling
tmp5_rolling0 = tmp_image(:,:,5:6:end);
tmp5_rolling1 = tmp_image(:,:,5:6:end-6);
tmp5_rolling2 = tmp_image(:,:,5:6:end-6);
tmp5_rolling3 = tmp_image(:,:,5:6:end-6);
tmp5_rolling4 = tmp_image(:,:,5:6:end-6);
tmp5_rolling5 = tmp_image(:,:,5+6:6:end);
tmp5 = cat(3,tmp5_rolling0, tmp5_rolling1, tmp5_rolling2, tmp5_rolling3, tmp5_rolling4, tmp5_rolling5);
clear tmp5_rolling1 tmp5_rolling2 tmp5_rolling3 tmp5_rolling4 tmp5_rolling5

% y3 with rolling
tmp6_rolling0 = tmp_image(:,:,6:6:end);
tmp6_rolling1 = tmp_image(:,:,6:6:end-6);
tmp6_rolling2 = tmp_image(:,:,6:6:end-6);
tmp6_rolling3 = tmp_image(:,:,6:6:end-6);
tmp6_rolling4 = tmp_image(:,:,6:6:end-6);
tmp6_rolling5 = tmp_image(:,:,6:6:end-6);
tmp6 = cat(3,tmp6_rolling0, tmp6_rolling1, tmp6_rolling2, tmp6_rolling3, tmp6_rolling4, tmp6_rolling5);
clear tmp6_rolling1 tmp6_rolling2 tmp6_rolling3 tmp6_rolling4 tmp6_rolling5

image_stack = zeros(size(tmp1,1), size(tmp1,2), size(tmp1,3), 6);
image_stack(:,:,:,1) = tmp1;
image_stack(:,:,:,2) = tmp2;
image_stack(:,:,:,3) = tmp3;
image_stack(:,:,:,4) = tmp4;
image_stack(:,:,:,5) = tmp5;
image_stack(:,:,:,6) = tmp6;
clear tmp_image tmp1 tmp2 tmp3 tmp4 tmp5 tmp6

image_all = image_stack;
clear image_stack



frame_num_round = size(image_all,3)


%% 3. Segmentation, select candidate points that meet the requirements
total_frame = size(image_all,3);
max_frame = 20000;   % Process 20000 frames at a time
if mod(total_frame,max_frame)==0
    iter = total_frame /max_frame;
else
    iter =  floor(total_frame /max_frame)+1;
end

setup.is_imgsz = 1 ;
setup.is_sCMOS = 0;
setup.offset = 104;  % the offset of the sCMOS
setup.gain = 1;
% image_stack = image_stack1+image_stack2+image_stack3+image_stack4+image_stack5+image_stack6;

% Start segmentation
% boxsz = 15;
% thresh_dist = 12;
boxsz = 17;  % The box size of each sub-image
thresh_dist = 13.5;  % The minimum dist of adjacent candidates, to avoid cross-talk
thresh_low = 105;   % The parameter to adjust segmentation threshold, increasing this value can segment out weaker candidates.
thresh_high = 105;
thresh = [thresh_low thresh_high];
r_boxsz = (boxsz-1)/2;

off = (size(PSF,1)-boxsz)/2

x_list = zeros(0,1);
y_list = zeros(0,1);
frame = zeros(0,1);
image_stack = zeros(boxsz,boxsz,0,6);
for i = 1:iter
    if i~=iter
        fanwei = (i-1)*max_frame+1 : i*max_frame;
    else
        fanwei = (i-1)*max_frame+1 : total_frame;
       
    end 
    tic
    [t_image_stack, t_seg_display] = crop_subregion_ast6N(image_all(:,:,fanwei,:),boxsz,thresh,thresh_dist, setup);
    toc
    image_stack = cat(3,image_stack,t_image_stack);
    x_list = cat(1,x_list,t_seg_display.allcds_mask(:,2));
    y_list = cat(1,y_list,t_seg_display.allcds_mask(:,1));
%     t_seg_display.allcds_mask(:,3) = t_seg_display.allcds_mask(:,3)+(i-1)*max_frame;
    frame = cat(1,frame,t_seg_display.allcds_mask(:,3)+(i-1)*max_frame);
    
end
x_list = x_list+5; % -1? =========================================================================================
y_list = y_list+5; % -1? =========================================================================================

frame = frame+1;

% Show segemntation result
num_display = 1+(0)*6;  % frame index  
if num_display < 1 || num_display > size(t_seg_display.ims_ch1,3)
    msgbox('Please input the correct number!');
end
disp('Show segmentation results');
raw = t_seg_display.ims_ch1(:,:,num_display)/max(max(t_seg_display.ims_ch1(:,:,num_display)));
index_rec = find(t_seg_display.allcds_mask(:,3) == num_display-1);
rec_vector = cat(2,t_seg_display.t1(index_rec),t_seg_display.l1(index_rec),repmat(boxsz,length(index_rec),2));
img_select = insertShape(raw,'Rectangle',rec_vector,'LineWidth',1, 'Color', 'green');
figure; imshow(img_select);
axis tight
title('Segmentation results');


a = sum(image_stack(:,:,1:100,:),4);
imageslicer(single(a));

total_candidates = size(x_list,1)






%% 4. Initial estimation
locs = total_candidates;
batchsize = 500000;
it_all = floor(total_candidates/batchsize)+1;
image_stack_sum = sum(image_stack,4);
sCMOSvarmap = 0;
Pcspline = zeros(locs, 6);
for seq = 1:it_all
    if seq == it_all
        range = batchsize*(seq-1)+1:locs;
    else
        range = batchsize*(seq-1)+1:batchsize*seq;
    end
    tic
    [t_Pcspline,CRLB]=mleFit_LM( image_stack_sum(:, :, range), 6, 50,single(coeff),sCMOSvarmap,1);
    toc
    Pcspline(range, :) =  t_Pcspline;
end






%% 5. Ratiometric estimation
result_init = zeros(locs,5); % N, x, y, z, bg
result_init(:,1) = Pcspline(:,3);
result_init(:,2) = Pcspline(:,1);
result_init(:,3) = Pcspline(:,2);
result_init(:,4) = Pcspline(:,5);
result_init(:,5) = Pcspline(:,4);
result_init = single(result_init);
data_x1 = single(image_stack(:, :, :, 1));
data_x2 = single(image_stack(:, :, :, 2));
data_x3 = single(image_stack(:, :, :, 3));
data_x4 = single(image_stack(:, :, :, 4));
data_x5 = single(image_stack(:, :, :, 5));
data_x6 = single(image_stack(:, :, :, 6));
coeff = single(coeff);
x_start = single(x_list- r_boxsz);
y_start = single(y_list- r_boxsz);

P2 = zeros(locs, 16);
for seq = 1:it_all
    if seq == it_all
        range = batchsize*(seq-1)+1:locs;
    else
        range = batchsize*(seq-1)+1:batchsize*seq;
    end
    tic
    [t_P2,Like2] = csplineMexCuda_6N(data_x1(:, :, range),data_x2(:, :, range),data_x3(:, :, range),...
        data_x4(:, :, range),data_x5(:, :, range),data_x6(:, :, range), result_init(range, :), (coeff), 50);
    toc
    P2(range, :) = t_P2;
end






%% 6. Preprocessing for estimating pattern parameters
result_6N = double(P2);
tx = result_6N(:,7)+ x_start;   % tx is the x coordinate and is already in pixel coordinates (global pixel coordinates)
ty = result_6N(:,8)+ y_start;   % 
tz = result_6N(:,9)- size(coeff,3)/2;   % 

tx_6bg = result_init(:,2)+ x_start;   % 
ty_6bg = result_init(:,3)+ y_start;   % 
tz_6bg = result_init(:,4)- size(coeff,3)/2;   % 





%% 7.1 iteration 1: Generate XY image
close all
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------%
mask_factor = 0.75;  % mask to determine the minimum estimate pattern period 

pixel = 5;  % The pixelsize of the rendered image
sigma = pixel;   % The sigma of the gaussian blob (unit: nm)
sigma_p = sigma/pixel;  % The sigma of the gaussian blob (unit: pixel)

r_ROI = floor(3*sigma_p)+1;
if mod(r_ROI,2)==0
    r_ROI = r_ROI+1;
end
gap = r_ROI*pixel+5;   % Leave blank area to rendered image，unit: nm
x1 = tx*x_pixelsize;   
y1 = ty*x_pixelsize;
x1 = x1 - min(x1) + gap;
y1 = y1 - min(y1) + gap;   
d = abs(max(x1)-max(y1))/2;
if max(x1) < max(y1)
    x1 = x1+d;
    xnum = floor( (max(x1)+ gap + d)/pixel ) + 2;
    ynum = floor( (max(y1)+ gap)/pixel ) + 2;
else
    y1 = y1+d;
    xnum = floor( (max(x1)+ gap )/pixel ) + 2;
    ynum = floor( (max(y1)+ gap + d)/pixel ) + 2;
end

image_x1 = zeros(xnum, ynum, 6);

N = result_6N(:,1:6);

[cy, cx] = meshgrid(-r_ROI:r_ROI, -r_ROI:r_ROI);
tmp_gauss = exp(-(cx.^2+cy.^2)/2/sigma_p.^2);
% figure;imshow(tmp_gauss,[]);title('Gaussian blob')

x11 = x1/pixel;  
y11 = y1/pixel;  
x_all = floor(x11);
y_all = floor(y11);   
x2 = x11 - x_all;
y2 = y11 - y_all;
tmp_stack = zeros(2*r_ROI+1, 2*r_ROI+1, size(x1,1), 6);

for j = 1:6
    parfor i = 1:size(x1,1)
        tmp_stack(:,:,i,j) = N(i,j)*exp(-( (cx-x2(i)).^2 + (cy-y2(i)).^2 )/2/sigma_p.^2);
    end
end

for j = 1:6
    for i = 1:size(x1,1)    
        image_x1(x_all(i)-r_ROI+1:x_all(i)+r_ROI+1, y_all(i)-r_ROI+1:y_all(i)+r_ROI+1, j) = ...
            image_x1(x_all(i)-r_ROI+1:x_all(i)+r_ROI+1, y_all(i)-r_ROI+1:y_all(i)+r_ROI+1, j)+tmp_stack(:,:,i,j);
    end
    tmp = image_x1(:,:,j);
    image_x1(:,:,j) = tmp./max(tmp(:));
end

% figure;
% subplot(231);imshow(image_x1(:,:,1),[]);title('Iteration1:image(x1)')
% subplot(232);imshow(image_x1(:,:,2),[]);title('Iteration1:image(x2)')
% subplot(233);imshow(image_x1(:,:,3),[]);title('Iteration1:image(x3)')
% subplot(234);imshow(image_x1(:,:,4),[]);title('Iteration1:image(y1)')
% subplot(235);imshow(image_x1(:,:,5),[]);title('Iteration1:image(y2)')
% subplot(236);imshow(image_x1(:,:,6),[]);title('Iteration1:image(y3)')


para_it1_xy = est_kvector(image_x1, 2,3, pixel, mask_factor, 1, lambda, NA);
% function para = est_kvector(raw_image, a_num, p_num, psize, mask_factor, flag, lambda, NA)
% flag = 1 represents x pattern only  or  xz/yz  or  x and y patterns estimation
% flag = 2 represents y pattern only estimation

% END iteration 1: Generate XY image
%====================================================================================================================================================================%





%% 7.2 iteration 1: Generate XZ image
close all
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------%
sitax_est = para_it1_xy(1,1)/180*pi;
vector = [cos(sitax_est), sin(sitax_est)];
x = tx*x_pixelsize; 
y = ty*x_pixelsize;
midx = (max(x)+min(x))/2;
midy = (max(y)+min(y))/2;
xlow = midx - xhalf_z*x_pixelsize;
xup = midx + xhalf_z*x_pixelsize;
ylow = midy - xhalf_z*x_pixelsize;
yup = midy + xhalf_z*x_pixelsize;
mask_spatial = x>xlow & x<xup & y>ylow & y<yup;
x = x(mask_spatial);
y = y(mask_spatial);
x1 = x*vector(1) + y*vector(2);

gap = r_ROI*pixel+5;   
gap_z = gap; 
 
z1 = tz(mask_spatial)*z_pixelsize;
x1 = x1 - min(x1) + gap;
z1 = z1 - min(z1) + gap_z;   
xnum = floor( (max(x1)+ gap)/pixel ) + 2;
znum = floor( (max(z1)+ gap_z)/pixel ) + 2;
if znum<xnum
   gap_z = floor(gap_z+ (xnum-znum)/2*pixel)+1;
end
x1 = x1 - min(x1) + gap;
z1 = z1 - min(z1) + gap_z;   
xnum = floor( (max(x1)+ gap)/pixel ) + 2;
znum = floor( (max(z1)+ gap_z)/pixel ) + 2;

image_xz = zeros(xnum, znum, 3);

N = result_6N(mask_spatial,1:3);

[cz, cx] = meshgrid(-r_ROI:r_ROI, -r_ROI:r_ROI);

x11 = x1/pixel;  
z11 = z1/pixel;  

x_all = floor(x11);
z_all = floor(z11);   
x2 = x11 - x_all;
z2 = z11 - z_all;
tmp_stack = zeros(2*r_ROI+1, 2*r_ROI+1, size(x1,1), 3);

for j = 1:3
    parfor i = 1:size(x1,1)
        tmp_stack(:,:,i,j) = N(i,j)*exp(-( (cx-x2(i)).^2 + (cz-z2(i)).^2 )/2/sigma_p.^2);
    end
end

for j = 1:3
    for i = 1:size(x1,1)    
        image_xz(x_all(i)-r_ROI+1:x_all(i)+r_ROI+1, z_all(i)-r_ROI+1:z_all(i)+r_ROI+1, j) = ...
            image_xz(x_all(i)-r_ROI+1:x_all(i)+r_ROI+1, z_all(i)-r_ROI+1:z_all(i)+r_ROI+1, j)+tmp_stack(:,:,i,j);
    end
    tmp = image_xz(:,:,j);
    image_xz(:,:,j) = tmp./max(tmp(:));
end

tmp_sum = sum(image_xz,3);
tmp_sum(tmp_sum==0) = 0.01;
tmp_xz = zeros(size(image_xz));
tmp_xz(:,:,1) = image_xz(:,:,1)./tmp_sum;
tmp_xz(:,:,2) = image_xz(:,:,2)./tmp_sum;
tmp_xz(:,:,3) = image_xz(:,:,3)./tmp_sum;
figure;
% subplot(231);imshow(image_xz(:,:,1),[]);title('Iteration1:image(xz1)')
% subplot(232);imshow(image_xz(:,:,2),[]);title('Iteration1:image(xz2)')
% subplot(233);imshow(image_xz(:,:,3),[]);title('Iteration1:image(xz3)')
subplot(234);imshow(tmp_xz(:,:,1),[]);
subplot(235);imshow(tmp_xz(:,:,2),[]);
subplot(236);imshow(tmp_xz(:,:,3),[]);

para_it1_xz = fdomain_est(tmp_xz, pixel, Tmin, mask_factor);


% END iteration 1: Generate XZ image
%====================================================================================================================================================================%






%% 7.3 iteration 1: Generate YZ image
close all
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------%
sitay_est = para_it1_xy(2,1)/180*pi;
vector = [cos(sitay_est), sin(sitay_est)];
x = tx*x_pixelsize; 
y = ty*x_pixelsize;
midx = (max(x)+min(x))/2;
midy = (max(y)+min(y))/2;
xlow = midx - xhalf_z*x_pixelsize;
xup = midx + xhalf_z*x_pixelsize;
ylow = midy - xhalf_z*x_pixelsize;
yup = midy + xhalf_z*x_pixelsize;
mask_spatial = x>xlow & x<xup & y>ylow & y<yup;
x = x(mask_spatial);
y = y(mask_spatial);
y1 = x*vector(1) + y*vector(2);


gap = r_ROI*pixel+5;  
gap_z = gap; 
z1 = tz(mask_spatial)*z_pixelsize;
y1 = y1 - min(y1) + gap;
z1 = z1 - min(z1) + gap_z;   
ynum = floor( (max(y1)+ gap)/pixel ) + 2;
znum = floor( (max(z1)+ gap_z)/pixel ) + 2;
if znum<ynum
   gap_z = floor(gap_z+ (ynum-znum)/2*pixel)+1;
end
y1 = y1 - min(y1) + gap;
z1 = z1 - min(z1) + gap_z;  
ynum = floor( (max(y1)+ gap)/pixel ) + 2;
znum = floor( (max(z1)+ gap_z)/pixel ) + 2;

image_yz = zeros(ynum, znum, 3);

N = result_6N(mask_spatial,4:6);

[cz, cy] = meshgrid(-r_ROI:r_ROI, -r_ROI:r_ROI);

y11 = y1/pixel;  
z11 = z1/pixel; 

y_all = floor(y11);
z_all = floor(z11);   
y2 = y11 - y_all;
z2 = z11 - z_all;
tmp_stack = zeros(2*r_ROI+1, 2*r_ROI+1, size(y1,1), 3);

for j = 1:3
    parfor i = 1:size(y1,1)
        tmp_stack(:,:,i,j) = N(i,j)*exp(-( (cy-y2(i)).^2 + (cz-z2(i)).^2 )/2/sigma_p.^2);
    end
end

for j = 1:3
    for i = 1:size(y1,1)    
        image_yz(y_all(i)-r_ROI+1:y_all(i)+r_ROI+1, z_all(i)-r_ROI+1:z_all(i)+r_ROI+1, j) = ...
            image_yz(y_all(i)-r_ROI+1:y_all(i)+r_ROI+1, z_all(i)-r_ROI+1:z_all(i)+r_ROI+1, j)+tmp_stack(:,:,i,j);
    end
    tmp = image_yz(:,:,j);
    image_yz(:,:,j) = tmp./max(tmp(:));
end


tmp_sum = image_yz(:,:,1)+image_yz(:,:,2)+image_yz(:,:,3);
tmp_sum(tmp_sum==0) = 0.01;
tmp_yz = zeros(size(image_yz));
tmp_yz(:,:,1) = image_yz(:,:,1)./tmp_sum;
tmp_yz(:,:,2) = image_yz(:,:,2)./tmp_sum;
tmp_yz(:,:,3) = image_yz(:,:,3)./tmp_sum;

figure;
% subplot(231);imshow(image_yz(:,:,1),[]);title('Iteration1:image(yz1)')
% subplot(232);imshow(image_yz(:,:,2),[]);title('Iteration1:image(yz2)')
% subplot(233);imshow(image_yz(:,:,3),[]);title('Iteration1:image(yz3)')
subplot(234);imshow(tmp_yz(:,:,1),[]);
subplot(235);imshow(tmp_yz(:,:,2),[]);
subplot(236);imshow(tmp_yz(:,:,3),[]);

para_it1_yz = fdomain_est(tmp_yz, pixel, Tmin, mask_factor);
% para = fdomain_est(image_stack, psize, Tmin, mask_factor)

% END iteration 1: Generate YZ image
%====================================================================================================================================================================%






%% 7.4 Iteration1 End
close all
sita_z1 = para_it1_xz(1,1)/180*pi;
sita_z2 = para_it1_yz(1,1)/180*pi;
vector_z1 = [cos(sita_z1), sin(sita_z1)];
vector_z2 = [cos(sita_z2), sin(sita_z2)];

para_all_it1 = [para_it1_xy(1,1),sita_z1/pi*180, para_it1_xz(1,2); para_it1_xy(2,1),sita_z2/pi*180, para_it1_yz(1,2)];






%% 8.1 iteration 2: Generate XY image
close all
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------%
% part1: First calculate the direction and period of pattern x
gap = r_ROI*pixel+5;   
tmp_x1 = tx*x_pixelsize;
tmp_z1 = tz*z_pixelsize;
x1 = tmp_x1*vector_z1(1) + tmp_z1*vector_z1(2);

y1 = ty*x_pixelsize;
x1 = x1 - min(x1) + gap;
y1 = y1 - min(y1) + gap;   
d = abs(max(x1)-max(y1))/2;
if max(x1) < max(y1)
    x1 = x1+d;
    xnum = floor( (max(x1)+ gap + d)/pixel ) + 2;
    ynum = floor( (max(y1)+ gap)/pixel ) + 2;
else
    y1 = y1+d;
    xnum = floor( (max(x1)+ gap )/pixel ) + 2;
    ynum = floor( (max(y1)+ gap + d)/pixel ) + 2;
end

image_x = zeros(xnum, ynum, 3);
N = result_6N(:,1:3);
[cy, cx] = meshgrid(-r_ROI:r_ROI, -r_ROI:r_ROI);

x11 = x1/pixel;  
y11 = y1/pixel;  

x_all = floor(x11);
y_all = floor(y11);   
x2 = x11 - x_all;
y2 = y11 - y_all;
tmp_stack = zeros(2*r_ROI+1, 2*r_ROI+1, size(x1,1), 6);

for j = 1:3
    parfor i = 1:size(x1,1)
        tmp_stack(:,:,i,j) = N(i,j)*exp(-( (cx-x2(i)).^2 + (cy-y2(i)).^2 )/2/sigma_p.^2);
    end
end

for j = 1:3
    for i = 1:size(x1,1)    
        image_x(x_all(i)-r_ROI+1:x_all(i)+r_ROI+1, y_all(i)-r_ROI+1:y_all(i)+r_ROI+1, j) = ...
            image_x(x_all(i)-r_ROI+1:x_all(i)+r_ROI+1, y_all(i)-r_ROI+1:y_all(i)+r_ROI+1, j)+tmp_stack(:,:,i,j);
    end
    tmp = image_x(:,:,j);
    image_x(:,:,j) = tmp./max(tmp(:));
end

% figure;
% subplot(131);imshow(image_x(:,:,1),[]);title('Iteration2:image(x1)')
% subplot(132);imshow(image_x(:,:,2),[]);title('Iteration2:image(x2)')
% subplot(133);imshow(image_x(:,:,3),[]);title('Iteration2:image(x3)')

para_it2_x = est_kvector(image_x, 1,3, pixel, mask_factor, 1, lambda, NA);


% part2:  Then calculate the direction and period of pattern y
gap = r_ROI*pixel+5;   
tmp_y1 = ty*x_pixelsize;
tmp_z1 = tz*z_pixelsize;
y1 = tmp_y1*vector_z2(1) + tmp_z1*vector_z2(2);

x1 = tx*x_pixelsize;
x1 = x1 - min(x1) + gap;
y1 = y1 - min(y1) + gap;   
d = abs(max(x1)-max(y1))/2;
if max(x1) < max(y1)
    x1 = x1+d;
    xnum = floor( (max(x1)+ gap + d)/pixel ) + 2;
    ynum = floor( (max(y1)+ gap)/pixel ) + 2;
else
    y1 = y1+d;
    xnum = floor( (max(x1)+ gap )/pixel ) + 2;
    ynum = floor( (max(y1)+ gap + d)/pixel ) + 2;
end

image_y = zeros(xnum, ynum, 3);

N = result_6N(:,4:6);
[cy, cx] = meshgrid(-r_ROI:r_ROI, -r_ROI:r_ROI);

x11 = x1/pixel;  
y11 = y1/pixel;  

x_all = floor(x11);
y_all = floor(y11);  
x2 = x11 - x_all;
y2 = y11 - y_all;
tmp_stack = zeros(2*r_ROI+1, 2*r_ROI+1, size(x1,1), 6);

for j = 1:3
    parfor i = 1:size(x1,1)
        tmp_stack(:,:,i,j) = N(i,j)*exp(-( (cx-x2(i)).^2 + (cy-y2(i)).^2 )/2/sigma_p.^2);
    end
end

for j = 1:3
    for i = 1:size(x1,1)    
        image_y(x_all(i)-r_ROI+1:x_all(i)+r_ROI+1, y_all(i)-r_ROI+1:y_all(i)+r_ROI+1, j) = ...
            image_y(x_all(i)-r_ROI+1:x_all(i)+r_ROI+1, y_all(i)-r_ROI+1:y_all(i)+r_ROI+1, j)+tmp_stack(:,:,i,j);
    end
    tmp = image_y(:,:,j);
    image_y(:,:,j) = tmp./max(tmp(:));
end

% figure;
% subplot(131);imshow(image_y(:,:,1),[]);title('Iteration2:image(y1)')
% subplot(132);imshow(image_y(:,:,2),[]);title('Iteration2:image(y2)')
% subplot(133);imshow(image_y(:,:,3),[]);title('Iteration2:image(y3)')

para_it2_y = est_kvector(image_y, 1,3, pixel, mask_factor, 2, lambda, NA);
% END iteration 2: Generate XY image
%====================================================================================================================================================================%




%% 8.2 iteration 2: Generate XZ image
close all
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------%
sitax_est = para_it2_x(1,1)/180*pi;
vector = [cos(sitax_est), sin(sitax_est)];
x = tx*x_pixelsize; 
y = ty*x_pixelsize;
midx = (max(x)+min(x))/2;
midy = (max(y)+min(y))/2;
xlow = midx - xhalf_z*x_pixelsize;
xup = midx + xhalf_z*x_pixelsize;
ylow = midy - xhalf_z*x_pixelsize;
yup = midy + xhalf_z*x_pixelsize;
mask_spatial = x>xlow & x<xup & y>ylow & y<yup;
x = x(mask_spatial);
y = y(mask_spatial);
x1 = x*vector(1) + y*vector(2);

gap = r_ROI*pixel+5;   
gap_z = gap; 
 
z1 = tz(mask_spatial)*z_pixelsize;
x1 = x1 - min(x1) + gap;
z1 = z1 - min(z1) + gap_z;   
xnum = floor( (max(x1)+ gap)/pixel ) + 2;
znum = floor( (max(z1)+ gap_z)/pixel ) + 2;
if znum<xnum
   gap_z = floor(gap_z+ (xnum-znum)/2*pixel)+1;
end
x1 = x1 - min(x1) + gap;
z1 = z1 - min(z1) + gap_z;  
xnum = floor( (max(x1)+ gap)/pixel ) + 2;
znum = floor( (max(z1)+ gap_z)/pixel ) + 2;

image_xz = zeros(xnum, znum, 3);

N = result_6N(mask_spatial,1:3);

[cz, cx] = meshgrid(-r_ROI:r_ROI, -r_ROI:r_ROI);


x11 = x1/pixel;  
z11 = z1/pixel;  

x_all = floor(x11);
z_all = floor(z11);   
x2 = x11 - x_all;
z2 = z11 - z_all;
tmp_stack = zeros(2*r_ROI+1, 2*r_ROI+1, size(x1,1), 3);

for j = 1:3
    parfor i = 1:size(x1,1)
        tmp_stack(:,:,i,j) = N(i,j)*exp(-( (cx-x2(i)).^2 + (cz-z2(i)).^2 )/2/sigma_p.^2);
    end
end

for j = 1:3
    for i = 1:size(x1,1)    
        image_xz(x_all(i)-r_ROI+1:x_all(i)+r_ROI+1, z_all(i)-r_ROI+1:z_all(i)+r_ROI+1, j) = ...
            image_xz(x_all(i)-r_ROI+1:x_all(i)+r_ROI+1, z_all(i)-r_ROI+1:z_all(i)+r_ROI+1, j)+tmp_stack(:,:,i,j);
    end
    tmp = image_xz(:,:,j);
    image_xz(:,:,j) = tmp./max(tmp(:));
end

tmp_sum = sum(image_xz,3);
tmp_sum(tmp_sum==0) = 0.01;
tmp_xz = zeros(size(image_xz));
tmp_xz(:,:,1) = image_xz(:,:,1)./tmp_sum;
tmp_xz(:,:,2) = image_xz(:,:,2)./tmp_sum;
tmp_xz(:,:,3) = image_xz(:,:,3)./tmp_sum;
figure;
% subplot(231);imshow(image_xz(:,:,1),[]);title('Iteration2:image(xz1)')
% subplot(232);imshow(image_xz(:,:,2),[]);title('Iteration2:image(xz2)')
% subplot(233);imshow(image_xz(:,:,3),[]);title('Iteration2:image(xz3)')
subplot(234);imshow(tmp_xz(:,:,1),[]);
subplot(235);imshow(tmp_xz(:,:,2),[]);
subplot(236);imshow(tmp_xz(:,:,3),[]);

para_it2_xz = fdomain_est(tmp_xz, pixel, Tmin, mask_factor);

% END iteration 2: Generate XZ image
%====================================================================================================================================================================%







% close all
%% 8.2_2 manually refine angle x
close all
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------%
Nbound = 1000;  % Remove all points whose total photon count are less than Nbound
sitax_est = para_it2_x(1,1)/180*pi;
vector = [cos(sitax_est), sin(sitax_est)];
x = tx*x_pixelsize; 
y = ty*x_pixelsize;
x1 = x*vector(1) + y*vector(2);
z1 = tz*z_pixelsize;

sitaz_est = para_it2_xz(1,1)/180*pi;
angle_estx = [sitax_est/pi*180, sitaz_est/pi*180];
vector = [cos(sitaz_est), sin(sitaz_est)];
x11 = x1*vector(1) + z1*vector(2);
N = result_6N(:,1:3);

N1 = N(:,1);
N2 = N(:,2);
N3 = N(:,3);
%--------------------------------------------------- Apply filter on photon count ------------------------------------------------------%
maskx = (N1+N2+N3)>Nbound;
N1 = N1(maskx);
N2 = N2(maskx);
N3 = N3(maskx);
tmp_x = x11(maskx);
% N_total = N1 + N2 + N3;
%------------------------------------------------- END: Apply filter on photon count ---------------------------------------------------%

N11 = N1./(N1+N2+N3);
N22 = N2./(N1+N2+N3);
N33 = N3./(N1+N2+N3);
figure;
subplot(311)
scatter(tmp_x, N11,1.5,'filled');
title('N1(x) case1')
subplot(312)
scatter(tmp_x, N22,1.5,'filled');
title('N2(x) case1')
subplot(313)
scatter(tmp_x, N33,1.5,'filled');
title('N3(x) case1')

promptx1 = "select the value of uplim/lowlim for x1&3, x2: ";
tmplimit = input(promptx1)/100;
uplimx13 = tmplimit(1);
lowlimx13 = tmplimit(2);
uplimx2 = tmplimit(3);
lowlimx2 = tmplimit(4);

mask_x1 = (N11<uplimx13 & N11>lowlimx13);    
mask_x2 = (N22<uplimx2 & N22>lowlimx2);
mask_x3 = (N33<uplimx13 & N33>lowlimx13);
mask_xall = mask_x1 & mask_x2 & mask_x3;
N11 = N11(mask_xall);
tmp_x = tmp_x(mask_xall);




%% 8.2_3 Adjust the azimuthal angle of x pattern
close all
angle_list = -2:0.05:2;
tmp_data = cell(1,numel(angle_list));
for ii = 1:numel(angle_list)
    sitax_est = (para_it2_x(1,1)+angle_list(ii))/180*pi;

    vector = [cos(sitax_est), sin(sitax_est)];
    x = tx*x_pixelsize;
    y = ty*x_pixelsize;
    x1 = x*vector(1) + y*vector(2);
    z1 = tz*z_pixelsize;

    sitaz_est = para_it2_xz(1,1)/180*pi;

    vector = [cos(sitaz_est), sin(sitaz_est)];
    x11 = x1*vector(1) + z1*vector(2);
    N = result_6N(:,1:3);

    N1 = N(:,1);
    N2 = N(:,2);
    N3 = N(:,3);
    tmp_maskx = (N1+N2+N3)>Nbound;
    N1 = N1(tmp_maskx);
    N2 = N2(tmp_maskx);
    N3 = N3(tmp_maskx);
    tmp_x = x11(tmp_maskx);

    N11 = N1./(N1+N2+N3);
    N22 = N2./(N1+N2+N3);
    N33 = N3./(N1+N2+N3);


    tmp_mask_x1 = (N11<uplimx13 & N11>lowlimx13);    
    tmp_mask_x2 = (N22<uplimx2 & N22>lowlimx2);
    tmp_mask_x3 = (N33<uplimx13 & N33>lowlimx13);
    tmp_mask_xall = tmp_mask_x1 & tmp_mask_x2 & tmp_mask_x3;
    N11 = N11(tmp_mask_xall);
%     N22 = N22(tmp_mask_xall);
%     N33 = N33(tmp_mask_xall);
    tmp_x = tmp_x(tmp_mask_xall);
    
    tmp_data{1,ii} = [tmp_x, N11];

end

slider_and_button_demo1(tmp_data, angle_list)

prompt = "select optimal angle(x) change: ";
delta_phase = input(prompt);
para_it2_x(1,1) = para_it2_x(1,1) + delta_phase;









%% 8.2_4 manually refine Tx
close all
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------%
T_list = -8:0.1:8;
tmp_data = cell(1,numel(T_list));
for ii = 1:numel(T_list)
    sitax_est = para_it2_x(1,1)/180*pi;
    T = para_it2_x(1,2) + T_list(ii);
    vector = [cos(sitax_est), sin(sitax_est)];
    x = tx*x_pixelsize;
    y = ty*x_pixelsize;
    x1 = x*vector(1) + y*vector(2);
    z1 = tz*z_pixelsize;

    sitaz_est = (para_it2_xz(1,1))/180*pi;

    vector = [cos(sitaz_est), sin(sitaz_est)];
    x11 = x1*vector(1) + z1*vector(2);
    N = result_6N(:,1:3);

    N1 = N(:,1);
    N2 = N(:,2);
    N3 = N(:,3);
    tmp_maskx = (N1+N2+N3)>Nbound;
    N1 = N1(tmp_maskx);
    N2 = N2(tmp_maskx);
    N3 = N3(tmp_maskx);
    tmp_x = x11(tmp_maskx);

    N11 = N1./(N1+N2+N3);
    N22 = N2./(N1+N2+N3);
    N33 = N3./(N1+N2+N3);


    tmp_mask_x1 = (N11<uplimx13 & N11>lowlimx13);   
    tmp_mask_x2 = (N22<uplimx2 & N22>lowlimx2);
    tmp_mask_x3 = (N33<uplimx13 & N33>lowlimx13);
    tmp_mask_xall = tmp_mask_x1 & tmp_mask_x2 & tmp_mask_x3;
    N11 = N11(tmp_mask_xall);
%     N22 = N22(tmp_mask_xall);
%     N33 = N33(tmp_mask_xall);
    tmp_x = tmp_x(tmp_mask_xall);
    
    tmp_data{1,ii} = [mod(tmp_x,T), N11];

end
slider_and_button_demo1(tmp_data, T_list)
prompt = "select optimal T change: ";
delta_phase = input(prompt);
para_it2_x(1,2) = para_it2_x(1,2) + delta_phase;

% END refine Tx manually
%====================================================================================================================================================================%










%% 8.3 iteration 2: Generate YZ image
close all
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------%
sitay_est = para_it2_y(1,1)/180*pi;
vector = [cos(sitay_est), sin(sitay_est)];
x = tx*x_pixelsize; 
y = ty*x_pixelsize;
midx = (max(x)+min(x))/2;
midy = (max(y)+min(y))/2;
xlow = midx - xhalf_z*x_pixelsize;
xup = midx + xhalf_z*x_pixelsize;
ylow = midy - xhalf_z*x_pixelsize;
yup = midy + xhalf_z*x_pixelsize;
mask_spatial = x>xlow & x<xup & y>ylow & y<yup;
x = x(mask_spatial);
y = y(mask_spatial);
y1 = x*vector(1) + y*vector(2);

gap = r_ROI*pixel+5;   
gap_z = gap; %******************************************************************************************************************************************** 
z1 = tz(mask_spatial)*z_pixelsize;
y1 = y1 - min(y1) + gap;
z1 = z1 - min(z1) + gap_z;   
ynum = floor( (max(y1)+ gap)/pixel ) + 2;
znum = floor( (max(z1)+ gap_z)/pixel ) + 2;
if znum<ynum
   gap_z = floor(gap_z+ (ynum-znum)/2*pixel)+1;
end
y1 = y1 - min(y1) + gap;
z1 = z1 - min(z1) + gap_z;   
ynum = floor( (max(y1)+ gap)/pixel ) + 2;
znum = floor( (max(z1)+ gap_z)/pixel ) + 2;

image_yz = zeros(ynum, znum, 3);

N = result_6N(mask_spatial,4:6);
[cz, cy] = meshgrid(-r_ROI:r_ROI, -r_ROI:r_ROI);
y11 = y1/pixel; 
z11 = z1/pixel;  

y_all = floor(y11);
z_all = floor(z11);  
y2 = y11 - y_all;
z2 = z11 - z_all;
tmp_stack = zeros(2*r_ROI+1, 2*r_ROI+1, size(y1,1), 3);

for j = 1:3
    parfor i = 1:size(y1,1)
        tmp_stack(:,:,i,j) = N(i,j)*exp(-( (cy-y2(i)).^2 + (cz-z2(i)).^2 )/2/sigma_p.^2);
    end
end

for j = 1:3
    for i = 1:size(y1,1)    
        image_yz(y_all(i)-r_ROI+1:y_all(i)+r_ROI+1, z_all(i)-r_ROI+1:z_all(i)+r_ROI+1, j) = ...
            image_yz(y_all(i)-r_ROI+1:y_all(i)+r_ROI+1, z_all(i)-r_ROI+1:z_all(i)+r_ROI+1, j)+tmp_stack(:,:,i,j);
    end
    tmp = image_yz(:,:,j);
    image_yz(:,:,j) = tmp./max(tmp(:));
end

% 20230714, revised xz/yz parameter estimation 
tmp_sum = image_yz(:,:,1)+image_yz(:,:,2)+image_yz(:,:,3);
tmp_sum(tmp_sum==0) = 0.01;
tmp_yz = zeros(size(image_yz));
tmp_yz(:,:,1) = image_yz(:,:,1)./tmp_sum;
tmp_yz(:,:,2) = image_yz(:,:,2)./tmp_sum;
tmp_yz(:,:,3) = image_yz(:,:,3)./tmp_sum;
figure;
% subplot(231);imshow(image_yz(:,:,1),[]);title('Iteration2:image(yz1)')
% subplot(232);imshow(image_yz(:,:,2),[]);title('Iteration2:image(yz2)')
% subplot(233);imshow(image_yz(:,:,3),[]);title('Iteration2:image(yz3)')
subplot(234);imshow(tmp_yz(:,:,1),[]);
subplot(235);imshow(tmp_yz(:,:,2),[]);
subplot(236);imshow(tmp_yz(:,:,3),[]);

para_it2_yz = fdomain_est(tmp_yz, pixel, Tmin, mask_factor);
% END iteration 2: Generate YZ image
%====================================================================================================================================================================%



%% 8.3-2 refine angle y manually
close all
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------%
sitay_est = para_it2_y(1,1)/180*pi;
vector = [cos(sitay_est), sin(sitay_est)];
x = tx*x_pixelsize; 
y = ty*x_pixelsize;
y1 = x*vector(1) + y*vector(2);
z1 = tz*z_pixelsize;

sitaz_est = para_it2_yz(1,1)/180*pi;
angle_esty = [sitay_est/pi*180, sitaz_est/pi*180];
vector = [cos(sitaz_est), sin(sitaz_est)];
y11 = y1*vector(1) + z1*vector(2);
N = result_6N(:,4:6);

N4 = N(:,1);
N5 = N(:,2);
N6 = N(:,3);
%--------------------------------------------------- Apply filter on photon count ------------------------------------------------------%
masky = (N4+N5+N6)>Nbound;
N4 = N4(masky);
N5 = N5(masky);
N6 = N6(masky);
tmp_y = y11(masky);
% N_total = N1 + N2 + N3;
%------------------------------------------------- END: Apply filter on photon count ---------------------------------------------------%

N44 = N4./(N4+N5+N6);
N55 = N5./(N4+N5+N6);
N66 = N6./(N4+N5+N6);

% figure;
% scatter(tmp_y, N44,1.5,'filled');
% title('N4(y) case1')

figure;
subplot(311)
scatter(tmp_y, N44,1.5,'filled');
title('N4(y) case1')
subplot(312)
scatter(tmp_y, N55,1.5,'filled');
title('N5(y) case1')
subplot(313)
scatter(tmp_y, N66,1.5,'filled');
title('N6(y) case1')


prompty1 = "select the value of uplim/lowlim for y1&3, y2: ";
tmplimity = input(prompty1)/100;
uplimy13 = tmplimity(1);
lowlimy13 = tmplimity(2);
uplimy2 = tmplimity(3);
lowlimy2 = tmplimity(4);

% uplim = 0.67;
% lowlim = 0.00;
mask_y1 = (N44<uplimy13 & N44>lowlimy13);   
mask_y2 = (N55<uplimy2 & N55>lowlimy2);
mask_y3 = (N66<uplimy13 & N66>lowlimy13);
mask_yall = mask_y1 & mask_y2 & mask_y3;
N44 = N44(mask_yall);
N55 = N55(mask_yall);
N66 = N66(mask_yall);
tmp_y = tmp_y(mask_yall);
% axiscoor = min(tmp_x):1:max(tmp_x);
N4 = N4(mask_yall);
N5 = N5(mask_yall);
N6 = N6(mask_yall);







%% 8.3-3 y角度调整
close all
tmp_data = cell(1,numel(angle_list));
for ii = 1:numel(angle_list)
    sitay_est = (para_it2_y(1,1)+angle_list(ii))/180*pi;

    vector = [cos(sitay_est), sin(sitay_est)];
    x = tx*x_pixelsize;
    y = ty*x_pixelsize;
    y1 = x*vector(1) + y*vector(2);
    z1 = tz*z_pixelsize;

    sitaz_est = para_it2_yz(1,1)/180*pi;

    vector = [cos(sitaz_est), sin(sitaz_est)];
    y11 = y1*vector(1) + z1*vector(2);
    N = result_6N(:,4:6);
    
    N4 = N(:,1);
    N5 = N(:,2);
    N6 = N(:,3);
    
    tmp_masky = (N4+N5+N6)>Nbound;
    N4 = N4(tmp_masky);
    N5 = N5(tmp_masky);
    N6 = N6(tmp_masky);
    tmp_y = y11(tmp_masky);
    N_total = N4 + N5 + N6;
    
    N44 = N4./(N4+N5+N6);
    N55 = N5./(N4+N5+N6);
    N66 = N6./(N4+N5+N6);


    tmp_mask_y1 = (N44<uplimy13 & N44>lowlimy13);    
    tmp_mask_y2 = (N55<uplimy2 & N55>lowlimy2);
    tmp_mask_y3 = (N66<uplimy13 & N66>lowlimy13);
    tmp_mask_yall = tmp_mask_y1 & tmp_mask_y2 & tmp_mask_y3;
    N44 = N44(tmp_mask_yall);
    N55 = N55(tmp_mask_yall);
    N66 = N66(tmp_mask_yall);
    tmp_y = tmp_y(tmp_mask_yall);
    
    tmp_data{1,ii} = [tmp_y, N44];

end

slider_and_button_demo1(tmp_data, angle_list)
prompt = "select optimal angle(y) change: ";
delta_phase = input(prompt);
para_it2_y(1,1) = para_it2_y(1,1) + delta_phase;
% END refine y angle manually
%====================================================================================================================================================================%




%% 8.3-4 refine Ty manually
close all
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------%
tmp_data = cell(1,numel(T_list));
for ii = 1:numel(T_list)
    sitay_est = (para_it2_y(1,1))/180*pi;
    T = T_list(ii)+para_it2_y(1,2);
    vector = [cos(sitay_est), sin(sitay_est)];
    x = tx*x_pixelsize;
    y = ty*x_pixelsize;
    y1 = x*vector(1) + y*vector(2);
    z1 = tz*z_pixelsize;

    sitaz_est = para_it2_yz(1,1)/180*pi;

    vector = [cos(sitaz_est), sin(sitaz_est)];
    y11 = y1*vector(1) + z1*vector(2);
    N = result_6N(:,4:6);
    
    N4 = N(:,1);
    N5 = N(:,2);
    N6 = N(:,3);
    
    tmp_masky = (N4+N5+N6)>Nbound;
    N4 = N4(tmp_masky);
    N5 = N5(tmp_masky);
    N6 = N6(tmp_masky);
    tmp_y = y11(tmp_masky);
    N_total = N4 + N5 + N6;
    
    N44 = N4./(N4+N5+N6);
    N55 = N5./(N4+N5+N6);
    N66 = N6./(N4+N5+N6);


    tmp_mask_y1 = (N44<uplimy13 & N44>lowlimy13);    
    tmp_mask_y2 = (N55<uplimy2 & N55>lowlimy2);
    tmp_mask_y3 = (N66<uplimy13 & N66>lowlimy13);
    tmp_mask_yall = tmp_mask_y1 & tmp_mask_y2 & tmp_mask_y3;
    N44 = N44(tmp_mask_yall);
    N55 = N55(tmp_mask_yall);
    N66 = N66(tmp_mask_yall);
    tmp_y = tmp_y(tmp_mask_yall);
    
    tmp_data{1,ii} = [mod(tmp_y,T), N44];

end

slider_and_button_demo1(tmp_data, T_list)
prompt = "select optimal T change: ";
delta_phase = input(prompt);
para_it2_y(1,2) = para_it2_y(1,2) + delta_phase;

% END refine yT angle manually
%====================================================================================================================================================================%






%% Iteration2 End
close all
para_all_it2 = [para_it2_x(1,1),para_it2_xz(1,1), para_it2_x(1,2); para_it2_y(1,1), para_it2_yz(1,1), para_it2_y(1,2)];
% Filter out the points that meet the criteria (points where both xy subimages are within on-state)
index_x = (1:locs)';
index_y = (1:locs)';
logic_x = zeros(locs,1);
logic_y = zeros(locs,1);
index_x = index_x(maskx);
index_x = index_x(mask_xall);
index_y = index_y(masky);
index_y = index_y(mask_yall);

logic_x(index_x,1)=1;
logic_y(index_y,1)=1;
logic_all = logic_y & logic_x;  




%% 9. Preprocessing for estimating global phase 
tx = double(tx(logic_all,1));
ty = double(ty(logic_all,1));
tz = double(tz(logic_all,1));
frame = frame(logic_all,1);  
tx_6bg = double(tx_6bg(logic_all,1));
ty_6bg = double(ty_6bg(logic_all,1));
tz_6bg = double(tz_6bg(logic_all,1));

result_6N = double(result_6N(logic_all,:));
result_init = double(result_init(logic_all,:));

% phase0 = 120/180*pi;
phase0 = phase_difference;
phasedata = zeros(size(result_6N,1), 8); %[phase1 phase2 md1 md2 amp1 amp2 offset1 offset2]
parfor ii =1:size(result_6N,1)
    intlist = result_6N(ii,1:6);
    [phase1, amp1, offset1] = MyCalPhase0828(intlist(1:3), phase0);
    [phase2, amp2, offset2] = MyCalPhase0828(intlist(4:6), phase0);
    phasedata(ii,:) = [phase1, phase2, amp1/offset1, amp2/offset2, amp1, amp2, offset1, offset2];
end
mean_mx = mean(phasedata(:,3))
mean_my = mean(phasedata(:,4))
std_mx = std(phasedata(:,3))
std_my = std(phasedata(:,4))
p1 = phasedata(:,1);   % p1 is phasex，ranging from -pi to pi
p2 = phasedata(:,2);   % p2 is phasey，ranging from -pi to pi

mean_intx = mean(sum(result_6N(:,1:3),2))*0.24
mean_inty = mean(sum(result_6N(:,4:6),2))*0.24




%% 10. Refine the pattern parameters using Least square fitting
rangez = 300;
mask_z = tz> -rangez/z_pixelsize & tz< rangez/z_pixelsize;
tx1 = tx(mask_z);
ty1 = ty(mask_z);
tz1 = tz(mask_z);

tic
init_parax = [para_all_it2(1,1)/180*pi, 2*pi/(para_all_it2(1,3)/x_pixelsize),  para_all_it2(1,2)/180*pi];   
% azimuthal angle，period，tilt angle

phase_list = (0)/180*pi;
z_list = (-8:0.5:8)/180*pi;
num_phase = numel(phase_list);
num_anlgez = numel(z_list);
error_list = zeros(num_anlgez,num_phase);

[phase, angle] = meshgrid(phase_list, z_list);
parfor ii = 1:num_anlgez
    for jj = 1:num_phase
        [~, error_list(ii,jj)] = FitPhasePlane3D0828_v2(tx1, ty1, tz1/x_pixelsize*z_pixelsize, p1(mask_z), ...
            [init_parax(1) init_parax(2) init_parax(3)+z_list(ii) phase_list(jj)], 1);
    end
end
tmp_m = (error_list == min(error_list));
[lpx, ~] = FitPhasePlane3D0828_v2(tx1, ty1, tz1/x_pixelsize*z_pixelsize, p1(mask_z), ...
        [init_parax(1) init_parax(2) init_parax(3)+mean(angle(tmp_m)) mean(phase(tmp_m))], 0); 

% angleZ_list = (-1:0.01:1)/180*pi;
% num_anlgeZ = numel(angleZ_list);
% error_list = zeros(num_anlgeZ,1);
% tic
% parfor ii = 1:num_anlgeZ
%     [~, error_list(ii)] = FitPhasePlane3D0828_v2(tx, ty, tz/x_pixelsize*z_pixelsize, p1, ...
%         [tmp_lpx(1) tmp_lpx(2) tmp_lpx(3)+angleZ_list(ii) tmp_lpx(4)], 1); 
% end
% toc
% select_anlgeZ = angleZ_list(error_list == min(error_list));
% [lpx, ~] = FitPhasePlane3D0828_v2(tx, ty, tz/x_pixelsize*z_pixelsize, p1, ...
%         [tmp_lpx(1) tmp_lpx(2) tmp_lpx(3)+select_anlgeZ tmp_lpx(4)], 1);
% [lpx(1)/pi*180, lpx(3)/pi*180, 2*pi/lpx(2)*x_pixelsize, lpx(4)/pi*180]


init_paray = [para_all_it2(2,1)/180*pi, 2*pi/(para_all_it2(2,3)/x_pixelsize),  para_all_it2(2,2)/180*pi];
% init_paray = [para_all_it2(2,1)/180*pi, 2*pi/(219.4299/x_pixelsize),  para_all_it2(2,2)/180*pi];
error_list = zeros(num_anlgez,num_phase);
parfor ii = 1:num_anlgez
    for jj = 1:num_phase
        [~, error_list(ii,jj)] = FitPhasePlane3D0828_v2(tx1, ty1, tz1/x_pixelsize*z_pixelsize, p2(mask_z), ...
            [init_paray(1) init_paray(2) init_paray(3)+z_list(ii) phase_list(jj)], 1); 
    end
end
tmp_m = (error_list == min(error_list));
[lpy, ~] = FitPhasePlane3D0828_v2(tx1, ty1, tz1/x_pixelsize*z_pixelsize, p2(mask_z), ...
        [init_paray(1) init_paray(2) init_paray(3)+mean(angle(tmp_m)) mean(phase(tmp_m))], 0); 

% [lpy, ~] = FitPhasePlane3D0828_v2(tx1, ty1, tz1/x_pixelsize*z_pixelsize, p2(mask_z), ...
%         [init_paray(1) init_paray(2) init_paray(3)+mean(angle(tmp_m)) -170/180*pi], 0); 

para_est = [lpx(1)/pi*180, lpx(3)/pi*180, 2*pi/lpx(2)*x_pixelsize, lpx(4)/pi*180;
lpy(1)/pi*180,  lpy(3)/pi*180, 2*pi/lpy(2)*x_pixelsize,lpy(4)/pi*180]
toc







%% 11. Reconstruction
result_init6N = zeros(size(result_6N,1),5); % N, x, y, z, bg
result_init6N(:,1) = result_6N(:,1)+result_6N(:,2)+result_6N(:,3)+result_6N(:,4)+result_6N(:,5)+result_6N(:,6);
result_init6N(:,2) = result_6N(:,7);
result_init6N(:,3) = result_6N(:,8);
result_init6N(:,4) = result_6N(:,9);
result_init6N(:,5) = result_6N(:,10)+result_6N(:,11)+result_6N(:,12)+result_6N(:,13)+result_6N(:,14)+result_6N(:,15);
result_init6N = single(result_init6N);

data_1 = single(data_x1(:, :, logic_all));
data_2 = single(data_x2(:, :, logic_all));
data_3 = single(data_x3(:, :, logic_all));
data_4 = single(data_x4(:, :, logic_all));
data_5 = single(data_x5(:, :, logic_all));
data_6 = single(data_x6(:, :, logic_all));
x_start1 = single(x_start(logic_all,:));
y_start1 = single(y_start(logic_all,:));

para = [lpx(2); lpy(2); lpx(1); lpy(1); lpx(3); lpy(3); lpx(4); lpy(4);];
para = single(para);
% para = [k_pixel; k_pixel;sita_x;sita_y;sita_z;sita_z;phase_x;phase_y];
% init_para = [para_all_it2(2,1)/180*pi, 2*pi/(para_all_it2(2,3)/x_pixelsize),  para_all_it2(2,2)/180*pi];
% [result_xy,~] = csplineMexCuda_xyz0311_v2(data_1,data_2,data_3,data_4,data_5,data_6, ...
%     result_init6N, para, single(coeff), x_start1, y_start1, 50, single(z_pixelsize), single(x_pixelsize));


result_xy = zeros(size(result_6N,1), 14);
it_all = floor(size(result_6N,1)/batchsize)+1;
for seq = 1:it_all
    if seq == it_all
        range = batchsize*(seq-1)+1:size(result_6N,1);
    else
        range = batchsize*(seq-1)+1:batchsize*seq;
    end
    [tresult_xy,~] = csplineMexCuda_xyz0311_v2(data_1(:,:,range),data_2(:,:,range),data_3(:,:,range),data_4(:,:,range),data_5(:,:,range),data_6(:,:,range), ...
    result_init6N(range,:), para, single(coeff), x_start1(range,:), y_start1(range,:), 50, single(z_pixelsize), single(x_pixelsize));
    result_xy(range, :) = tresult_xy;
end






%% 12. Convert to xy coordinates and save the data
final_x = (x_start1 + result_xy(:,5))*x_pixelsize;
final_y = (y_start1 + result_xy(:,6))*x_pixelsize;
final_z = (result_xy(:,7) - size(coeff,3)/2)*z_pixelsize;

N_x = tx*x_pixelsize;
N_y = ty*x_pixelsize;
N_z = tz*z_pixelsize;

bg_x = tx_6bg*x_pixelsize;
bg_y = ty_6bg*x_pixelsize;
bg_z = tz_6bg*z_pixelsize;

para_save_path = [save_path,'para\'];
result_save_path = [save_path,'result\'];


% Check and create the para save folder
if ~exist(para_save_path, 'dir')
    mkdir(para_save_path);
end

% Check and create the result save folder
if ~exist(result_save_path, 'dir')
    mkdir(result_save_path);
end


% save([para_save_path, 'para_b1_r1.mat'], 'lpx', 'lpy', 'para_est');
% save([result_save_path, 'result_b1_r1.mat'], 'frame', 'result_xy', 'result_6N', 'result_init','final_x', 'final_y', 'final_z', 'N_x', 'N_y', 'N_z','bg_x', 'bg_y', 'bg_z');  % 1126更新

























%% 15. Processing the data by each round, and save the corresponding results
% One can change batch_frame parameters if needed
% ===================================================================== Remember to revise these parameters ==================================================================================================================
total_batch = 12;
batch_frame = 30000;  % how many frames in each batch
round_size = 3000;  % for each batch, round_size frames are used to estimate local parameters
total_round = floor(batch_frame/round_size)  % number of round for each batch
Last_batch = 12;
Last_round = 10;   % This means the reconstruction will end at the Last_round of Last_batch

frame_num_round = round_size-5
% ==================================================================== END：Remember to revise these parameters ==================================================================================================================


lpx_last_round = lpx;
lpy_last_round = lpy;
for batch_idx = 1:Last_batch   
    if batch_idx == 1
        round_seq = 1:total_round;     
    else
        if batch_idx == Last_batch
            round_seq = 1:Last_round;
        else
            round_seq = 1:total_round;
        end
    end

    filepath = [file_title, num2str(batch_idx) ,'.tif'];
    tmp_image_all = loadtiff(filepath);
    tmp_image_all = double(tmp_image_all);

    for round_idx =  round_seq
        tmp_image = tmp_image_all(:,:,1+(round_idx-1)*round_size:round_idx*round_size);

        [lpx, lpy, para_est, result_xy, result_6N, result_init, final_x, final_y, final_z, N_x, N_y, N_z, bg_x, bg_y, bg_z, frame] = process_round_data(tmp_image, coeff, phase_difference, x_pixelsize, z_pixelsize, ...
            boxsz, thresh_dist, thresh_low, thresh_high, Nbound, uplimx13, lowlimx13, uplimx2, lowlimx2, uplimy13, lowlimy13, uplimy2, lowlimy2, ...
            rangez, lpx_last_round, lpy_last_round);
        
        frame_offset = ( (batch_idx-1)*total_round+round_idx-1)* frame_num_round;
        frame = frame + frame_offset;
        para_save_name_now = ['para_b',num2str(batch_idx),'_r',num2str(round_idx),'.mat'];
        result_save_name_now = ['result_b',num2str(batch_idx),'_r',num2str(round_idx),'.mat']
        save([para_save_path, para_save_name_now], 'lpx', 'lpy', 'para_est');
        save([result_save_path, result_save_name_now], 'frame', 'result_xy', 'result_6N', 'result_init','final_x', 'final_y', 'final_z', 'N_x', 'N_y', 'N_z','bg_x', 'bg_y', 'bg_z');
        
        % Update the pattern parameters
        lpx_last_round = lpx;
        lpy_last_round = lpy;


    end

end













