% This function is to process the 3D-SMILE data for each round of data
function [lpx, lpy, para_est, result_xy, result_6N, result_init, final_x, final_y, final_z, N_x, N_y, N_z, bg_x, bg_y, bg_z, frame] = process_round_data_modloc_single_V3(tmp_image, coeff, phase_difference, x_pixelsize, z_pixelsize, boxsz, thresh_dist, thresh_low, thresh_high, ...
    Nbound, uplimx13, lowlimx13, uplimx2, lowlimx2,...
    uplimy13, lowlimy13, uplimy2, lowlimy2, rangez,lpx_last_round, lpy_last_round)

    %% 1. Pre-process
    % x1 with rolling
    tmp1_rolling0 = tmp_image(:,:,1:3:end);     % 500 frames
    tmp1_rolling1 = tmp_image(:,:,1+3:3:end);   % 499 frames
    tmp1_rolling2 = tmp_image(:,:,1+3:3:end);   % 499 frames
    tmp1 = cat(3,tmp1_rolling0, tmp1_rolling1, tmp1_rolling2);
    clear tmp1_roliing0 tmp1_rolling1 tmp1_rolling2 
    
    % x2 with rolling
    tmp2_rolling0 = tmp_image(:,:,2:3:end);     % 500 frames
    tmp2_rolling1 = tmp_image(:,:,2:3:end-3);
    tmp2_rolling2 = tmp_image(:,:,2+3:3:end);
    tmp2 = cat(3,tmp2_rolling0, tmp2_rolling1, tmp2_rolling2);
    clear tmp2_rolling1 tmp2_rolling2 tmp2_rolling0 
    
    % x3 with rolling
    tmp3_rolling0 = tmp_image(:,:,3:3:end);
    tmp3_rolling1 = tmp_image(:,:,3:3:end-3);
    tmp3_rolling2 = tmp_image(:,:,3:3:end-3);
    tmp3 = cat(3,tmp3_rolling0, tmp3_rolling1, tmp3_rolling2);
    clear tmp3_rolling1 tmp3_rolling2 tmp3_rolling0
    
    
    image_stack = zeros(size(tmp1,1), size(tmp1,2), size(tmp1,3), 6);
    image_stack(:,:,:,1) = tmp1;
    image_stack(:,:,:,2) = tmp2;
    image_stack(:,:,:,3) = tmp3;
    image_stack(:,:,:,4) = tmp1;   % 1210 revised ------------------------------------------------------------------------------------------------------------------------------------------------
    image_stack(:,:,:,5) = tmp2;   % 1210 revised ------------------------------------------------------------------------------------------------------------------------------------------------
    image_stack(:,:,:,6) = tmp3;   % 1210 revised ------------------------------------------------------------------------------------------------------------------------------------------------
    clear tmp_image tmp1 tmp2 tmp3 
    
    image_all = image_stack;
    clear image_stack
    
    
    
    frame_num_round = size(image_all,3)






    %% 3. 子图分割，选择符合要求的候选点
    total_frame = size(image_all,3);
    max_frame = 20000;   % 单次处理10000frame
    if mod(total_frame,max_frame)==0
        iter = total_frame /max_frame;
    else
        iter =  floor(total_frame /max_frame)+1;
    end

    setup.is_imgsz = 1 ;
    setup.is_sCMOS = 0;
    setup.offset = 104;
    setup.gain = 1;
    % image_stack = image_stack1+image_stack2+image_stack3+image_stack4+image_stack5+image_stack6;

    % Start segmentation
%     boxsz = 21;
%     thresh_dist = 15;
%     thresh_low = 105;   % 提高这个数值能够分割出更多的点
%     thresh_high = 105;
    thresh = [thresh_low thresh_high];
    r_boxsz = (boxsz-1)/2;

%     off = (size(PSF,1)-boxsz)/2

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
    % 注意这里的加5是因为seg code里面有is_imgsz, 并且不用减1
    frame = frame+1;

%     % Show segemntation result
%     num_display = 200;  % frame index
%     if num_display < 1 || num_display > size(t_seg_display.ims_ch1,3)
%         msgbox('Please input the correct number!');
%     end
%     disp('Show segmentation results');
%     raw = t_seg_display.ims_ch1(:,:,num_display)/max(max(t_seg_display.ims_ch1(:,:,num_display)));
%     index_rec = find(t_seg_display.allcds_mask(:,3) == num_display-1);
%     rec_vector = cat(2,t_seg_display.t1(index_rec),t_seg_display.l1(index_rec),repmat(boxsz,length(index_rec),2));
%     img_select = insertShape(raw,'Rectangle',rec_vector,'LineWidth',1, 'Color', 'green');
%     figure; imshow(img_select);
%     axis tight
%     title('Segmentation results');
% 
% 
%     a = sum(image_stack(:,:,1:100,:),4);
%     imageslicer(single(a));

    total_candidates = size(x_list,1)





    
    %% 4. 初始估计
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






    %% 5. 6N 估计
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






    %% 6. 进行预处理和Iteration 1
    result_6N = double(P2);
    tx = result_6N(:,7)+ x_start;   % tx为x的坐标，并且已经是像素坐标了(全局像素坐标)
    ty = result_6N(:,8)+ y_start;   % ty为y的坐标，并且已经是像素坐标了(全局像素坐标)
    tz = result_6N(:,9)- size(coeff,3)/2;   % tz为z的坐标，并且已经是像素坐标了(全局像素坐标)

    tx_6bg = result_init(:,2)+ x_start;   % tx为x的坐标，并且已经是像素坐标了(全局像素坐标)
    ty_6bg = result_init(:,3)+ y_start;   % ty为y的坐标，并且已经是像素坐标了(全局像素坐标)
    tz_6bg = result_init(:,4)- size(coeff,3)/2;   % tz为z的坐标，并且已经是像素坐标了(全局像素坐标)






    %% 7. Generate relevant mask and select locs
    % 7.1 xmasks
    N = result_6N(:,1:3);
    N1 = N(:,1);
    N2 = N(:,2);
    N3 = N(:,3);

    maskx = (N1+N2+N3)>Nbound;
    N1 = N1(maskx);
    N2 = N2(maskx);
    N3 = N3(maskx);

    N11 = N1./(N1+N2+N3);
    N22 = N2./(N1+N2+N3);
    N33 = N3./(N1+N2+N3);

    mask_x1 = (N11<uplimx13 & N11>lowlimx13);    % 这里改成mask_x1，是为了之后的
    mask_x2 = (N22<uplimx2 & N22>lowlimx2);
    mask_x3 = (N33<uplimx13 & N33>lowlimx13);
    mask_xall = mask_x1 & mask_x2 & mask_x3;



    % 7.2 ymasks
    N = result_6N(:,4:6);

    N4 = N(:,1);
    N5 = N(:,2);
    N6 = N(:,3);
    %--------------------------------------------------- 加光子数bound ------------------------------------------------------%
    masky = (N4+N5+N6)>Nbound;
    N4 = N4(masky);
    N5 = N5(masky);
    N6 = N6(masky);
    % N_total = N1 + N2 + N3;
    %------------------------------------------------- END: 加光子数bound ---------------------------------------------------%
    
    N44 = N4./(N4+N5+N6);
    N55 = N5./(N4+N5+N6);
    N66 = N6./(N4+N5+N6);

    mask_y1 = (N44<uplimy13 & N44>lowlimy13);   
    mask_y2 = (N55<uplimy2 & N55>lowlimy2);
    mask_y3 = (N66<uplimy13 & N66>lowlimy13);
    mask_yall = mask_y1 & mask_y2 & mask_y3;


    % 7.3 筛选下符合条件的点(xy子图都亮的点)
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
    logic_all = logic_y & logic_x;  % logic_all 就是筛选出来的候选点







    %% 8. 进一步估算全局phase
    tx = double(tx(logic_all,1));
    ty = double(ty(logic_all,1));
    tz = double(tz(logic_all,1));
    frame = frame(logic_all,1);  % +++++++++++++++++++++++++++++++++ 20251119 新增: 把frame也相应的修改 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
    p1 = phasedata(:,1);   % p1为phasex，在-pi到pi之间
    p2 = phasedata(:,2);   % p2为phasey，在-pi到pi之间





    %% 9. 进一步refine 最后的结果
%     rangez = 300;
    mask_z = tz> -rangez/z_pixelsize & tz< rangez/z_pixelsize;
    tx1 = tx(mask_z);
    ty1 = ty(mask_z);
    tz1 = tz(mask_z);
    
    tic
    [lpx, ~] = FitPhasePlane3D0828_v2(tx1, ty1, tz1/x_pixelsize*z_pixelsize, p1(mask_z), ...
            [lpx_last_round(1) lpx_last_round(2) lpx_last_round(3) lpx_last_round(4)], 0); 
    

    [lpy, ~] = FitPhasePlane3D0828_v2(tx1, ty1, tz1/x_pixelsize*z_pixelsize, p2(mask_z), ...
            [lpy_last_round(1) lpy_last_round(2) lpy_last_round(3) lpy_last_round(4)], 0); 
    

    
    para_est = [lpx(1)/pi*180, lpx(3)/pi*180, 2*pi/lpx(2)*x_pixelsize, lpx(4)/pi*180;
    lpy(1)/pi*180,  lpy(3)/pi*180, 2*pi/lpy(2)*x_pixelsize,lpy(4)/pi*180]
    toc







    %% 10. Reconstruction
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


    

    %% 20260128 利用modloc的方法进行处理(只针对Y方向进行处理) ==================================================================================================================================
    final_x = (x_start1 + result_xy(:,5))*x_pixelsize;
    final_y = (y_start1 + result_xy(:,6))*x_pixelsize;
    % final_z = (result_xy(:,7) - size(coeff,3)/2)*z_pixelsize;  
    
    tz_modify = tz/x_pixelsize*z_pixelsize;
    
    %lp: (theta, cycle, phase)
    theta_1 = lpx(1);
    cycle_1 = lpx(2);
    thetaz_1 = lpx(3);
    phase_1 = lpx(4);
    result_p1 = (tx.*cos(theta_1).*cos(thetaz_1) + ty.*sin(theta_1).*cos(thetaz_1)+tz_modify .*sin(thetaz_1)) * cycle_1 + phase_1;
    
    k_num = round((result_p1- p1)/2/pi);
    
    
    T1 = 2*pi/lpx(2)*x_pixelsize;
    p1_vector = [cos(theta_1),sin(theta_1)];
    r0 = p1_vector(1)*final_x + p1_vector(2)*final_y;
    
    
    fine_z2 = (k_num + p1/2/pi - phase_1/2/pi)*T1/sin(thetaz_1)-r0/tan(thetaz_1);
    
    final_z = fine_z2;
    
    
    
    N_x = tx*x_pixelsize;
    N_y = ty*x_pixelsize;
    N_z = tz*z_pixelsize;
    
    bg_x = tx_6bg*x_pixelsize;
    bg_y = ty_6bg*x_pixelsize;
    bg_z = tz_6bg*z_pixelsize;



end