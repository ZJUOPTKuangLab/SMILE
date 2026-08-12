%% This is for stitching all results together
close all
clear
clc


% ===================================================================== Remember to revise these parameters ==================================================================================================================
total_batch = 10;
batch_frame = 30000;  % how many frames in each batch
round_size = 3000;  % for each batch, round_size frames are used to estimate local parameters
total_round = batch_frame/round_size;  % number of round for each batch
file_title = 'Template data for tilt SMILE\result\result_';
Last_batch = 10;
Last_round = 10;   % This means the reconstruction will end at the Last_round of Last_batch
save_img_path = 'Template data for tilt SMILE\';
% ==================================================================== END：Remember to revise these parameters ==================================================================================================================


zrange = 600;
m = 0.4;
factor = 2.5;
N = (3000)/3/0.24   % photon count filter


batch_list = 1:Last_batch;
start_idx = batch_list(1);

for batch_idx = batch_list  
   
    if batch_idx == Last_batch
        round_seq = 1:Last_round;
  else
        round_seq = 1:total_round;
    end

    for round_idx =  round_seq
        filepath = [file_title, 'b',num2str(batch_idx),'_r', num2str(round_idx)];
        load(filepath)

        result_xy(:,5) = final_x;
        result_xy(:,6) = final_y;
        result_xy(:,7) = final_z;

        result_6N(:,7) = N_x;
        result_6N(:,8) = N_y;
        result_6N(:,9) = N_z;

        result_init(:,2) = bg_x;
        result_init(:,3) = bg_y;
        result_init(:,4) = bg_z;

        mask = final_z<zrange & final_z>-zrange & result_xy(:,3)>m & result_xy(:,4)>m & result_xy(:,1)./result_xy(:,2)<factor & result_xy(:,1)./result_xy(:,2)>1/factor...
            &(result_xy(:,1)+result_xy(:,2))>N;
        reserve_ratio = sum(mask)/size(mask,1)*100
        frame = frame(mask);

        result_xy = result_xy(mask,:);
        result_6N = result_6N(mask,:);
        result_init = result_init(mask,:);

        if batch_idx == start_idx && round_idx ==1
            result_xy_all = result_xy;
            result_6N_all = result_6N;
            result_init_all = result_init;
            frame_all = frame;
        else
            result_xy_all = cat(1,result_xy_all,result_xy);
            result_6N_all = cat(1,result_6N_all,result_6N);
            result_init_all = cat(1,result_init_all,result_init);
            frame_all = cat(1,frame_all,frame);
        end

    end

end

%% analyze 
N = result_init_all(:,1)*0.24;
bg = (result_init_all(:,5)-6*104)*0.24;
avg_N = mean(N)
avg_bg = mean(bg)

mx = result_xy_all(:,3);
my = result_xy_all(:,4);
avg_mx = mean(mx)
avg_my = mean(my)

figure;
hist(N,100);title('Photons')

figure;
hist(mx,100);title('mx')

figure;
hist(my,100);title('my')




%% save 
frame_all = frame_all - min(frame_all) +1;
start_frame = min(frame_all)
save([save_img_path,'result_init_all.mat'], 'result_init_all');
save([save_img_path,'result_6N_all.mat'], 'result_6N_all');
save([save_img_path,'result_xy_all.mat'], 'result_xy_all');
save([save_img_path,'frame_all.mat'], 'frame_all');


