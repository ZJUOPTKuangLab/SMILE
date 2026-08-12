function para = fdomain_est(image_stack, psize, Tmin, mask_factor)
    [xsize, ysize, zsize] = size(image_stack);
    f = zeros(xsize,ysize);
    for i = 1:zsize
        tmp = image_stack(:,:,i);
        tmpf = abs(fftshift(fftn(tmp)));
        f = f+tmpf;
    end
    f = (f - min(f(:)))/(max(f(:)) - min(f(:)));
    shiftvalue = zeros(1,2);

    delta_f = 1/psize/(xsize-1);
    T_f = 1/Tmin;
    cutoff = T_f/delta_f*mask_factor;
    Tlimit = Tmin/mask_factor

    [X,Y]=meshgrid(1:ysize,1:xsize);
    xc=floor(xsize/2+1);% the x-coordinate of the center
    yc=floor(ysize/2+1);% the y-coordinate of the center
    yr=Y-yc;
    xr=X-xc;
    R=sqrt((xr).^2+(yr).^2);% distance between the point (x,y) and center (xc,yc)

    fmask = zeros(xsize, ysize);
    fmask(R>cutoff)=1;
    fmask(Y>yc)=0;
    [tmp1, tmp2]  = find( f.*fmask==max(max(f.*fmask)) );
%     [shiftvalue(1),shiftvalue(2)]= find( f.*fmask==max(max(f.*fmask)) );
    shiftvalue(1) = mean(tmp1);
    shiftvalue(2) = mean(tmp2);



    cutoff=round(cutoff);
    circ_x=-cutoff:cutoff;
    circ_y(1:2*cutoff+1)=sqrt(cutoff^2-circ_x.^2);
    circ_y(2*cutoff+2:4*cutoff+2)=-flip(circ_y(1:2*cutoff+1));
    circ_x=[circ_x,flip(circ_x)];
    circ_x(4*cutoff+3)=circ_x(1);
    circ_y(4*cutoff+3)=circ_y(1);

%     circ_x=circ_x+shiftvalue(2);
%     circ_y=circ_y+shiftvalue(1);

    circ_x=circ_x+xc;
    circ_y=circ_y+yc;

    box_x=[3,3,-3,-3,3];
    box_y=[3,-3,-3,3,3];


    figure;imagesc(f);
    hold on;
    plot(circ_x,circ_y,'--g',...
        box_x+shiftvalue(2),box_y+shiftvalue(1),'-r');
%     colormap('parula')
    colormap('gray')
    % the component within the green dashed circle is excluded when estimating k0
    % the center of the red solid box denotes the estimated k0

    %% 下面是计算更加精准的值，对原始图像进行cubic插值，然后再去找极值点
%     r_ROI = 2;
%     tmp = f(shiftvalue(1)-r_ROI:shiftvalue(1)+r_ROI, shiftvalue(2)-r_ROI:shiftvalue(2)+r_ROI);
%     F = imresize(tmp,[25 25],'cubic');
% %     figure;imshow(F,[]);
%     mid = (25+1)/2;
%     shiftvalue1 = zeros(1,2);
%     [shiftvalue1(1),shiftvalue1(2)]= find( F == max(F(:)) );
%     shiftvalue1 = (shiftvalue1 - mid)/5;   

    
    precise_shift = shiftvalue;
    precise_shift(1) = precise_shift(1) - xc;
    precise_shift(2) = precise_shift(2) - yc;
%     shiftvalue(1) = shiftvalue(1) - xc;
%     shiftvalue(2) = shiftvalue(2) - yc;
    para = zeros(1,2);  % angle, T
    angle = 90 - atan( precise_shift(1)/(precise_shift(2)) )/pi*180;
%     angle1 = 90 - atan( shiftvalue(1)/(shiftvalue(2)) )/pi*180;
    if angle>90
        angle = angle - 180;
    end
    angle = angle
    kmin = 2*pi/psize/xsize;
    k_shift = kmin*( precise_shift(1)*precise_shift(1)+precise_shift(2)*precise_shift(2) )^(0.5);
%     k_shift1 = kmin*( shiftvalue(1)*shiftvalue(1)+shiftvalue(2)*shiftvalue(2) )^(0.5);
    T = 2*pi/k_shift
%     T1 = 2*pi/k_shift1;
    para(1,1) = angle;
    para(1,2) = T;
%     para(2,1) = angle1;
%     para(2,2) = T1;




end