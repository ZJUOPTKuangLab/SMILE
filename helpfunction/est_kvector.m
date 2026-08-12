function para = est_kvector(raw_image, a_num, p_num, psize, mask_factor, flag, lambda, NA)
% flag = 1 代表 x单方向或者xz/yz或者xy
% flag = 2 代表 y单方向
    
    %% parameter of the detection system
%     lambda=560;% fluorescence emission wavelength (emission maximum). unit: nm
% %     psize = 5; % psize=pixel size/magnification power. unit: nm
%     NA = 1.45;
    
    %% parameter for reconstruction
    wiener_factor=0.005;
    
%     mask_factor=0.75;%a high-pass mask (fmask) is utilized to estimate the modulation vector;
    % the cutoff frequency of fmask is mask_factor*(cutoff frequency of the detection OTF)
    % recommended value: 0.6 for conventional SIM, 0.8 for TIRF-SIM
    
    pitch = lambda/NA/mask_factor/2
    
    %% visualization option
    show_initial_result_flag=1;
    % the phase-only correlation result without correction will
    % be displayed when show_initial_result_flag equals 1
    
    show_corrected_result_flag=1;
    % the phase-only correlation result after correction will
    % be displayed when show_corrected_result_flag equals 1
    
    %% saving file
    save_flag=0; % save the results if save_flag equals 1;
    
    for ii=1:a_num
        for jj=1:p_num
            noiseimage(:,:,ii,jj)=raw_image(:,:,(ii-1)*p_num+jj);        
        end
    end


%     figure;imshow(noiseimage(:,:,1,1),[]);title('raw image');
%     tmp = noiseimage(:,:,1,1);
%     tmp_f = fftshift(fft2(tmp));
%     figure;imshow((abs(tmp_f)).^0.15,[]);title('Fourier domain')
    

    [xsize,ysize]=size(noiseimage(:,:,1,1));
    % make square image
    if xsize~=ysize
        if xsize<ysize
            noiseimage = noiseimage(1:xsize, 1:xsize, :,:);
        else
            noiseimage = noiseimage(1:ysize, 1:ysize, :,:);
        end
    end
    [xsize,ysize]=size(noiseimage(:,:,1,1));
    
    
    
    [X,Y]=meshgrid(1:ysize,1:xsize);
    
    PSF_edge = fspecial('gaussian',5,40);
    for ii=1:a_num
        for jj=1:p_num
            noiseimage(:,:,ii,jj)=edgetaper(noiseimage(:,:,ii,jj),PSF_edge);
        end
    end
    
    
    xc=floor(xsize/2+1);% the x-coordinate of the center
    yc=floor(ysize/2+1);% the y-coordinate of the center
    yr=Y-yc;
    xr=X-xc;
    R=sqrt((xr).^2+(yr).^2);% distance between the point (x,y) and center (xc,yc)
    %% Generate the PSF
    pixelnum=xsize;
    rpixel=NA*pixelnum*psize/lambda;
    cutoff=round(2*rpixel);% cutoff frequency
    % ctfde=ones(pixelnum,pixelnum).*(R<=rpixel);
    ctfde=ones(xsize,ysize).*(R<=rpixel);
    ctfdeSignificantPix=numel(find(abs(ctfde)>eps(class(ctfde))));
    ifftscalede=numel(ctfde)/ctfdeSignificantPix;
    apsfde=fftshift(ifft2(ifftshift(ctfde)));
    ipsfde=ifftscalede*abs(apsfde).^2;
    OTFde=real(fftshift(fft2(ifftshift(ipsfde))));
    clear apsfde ctfde temp X Y
    %% filter/deconvolution before using noiseimage
    widefield=sum(sum(noiseimage,4),3);
    widefield=quasi_wnr(OTFde,widefield,wiener_factor^2);
    widefield=widefield.*(widefield>0);
    
    for ii=1:a_num
        for jj=1:p_num
            noiseimage(:,:,ii,jj)=quasi_wnr(OTFde,squeeze(noiseimage(:,:,ii,jj)),wiener_factor^2);
    
            %noiseimage(:,:,ii,jj)=deconvlucy(noiseimage(:,:,ii,jj),ipsfde,3);
            %pre-deconvolution. It can be applied to suppress noises in experiments
    
            noiseimage(:,:,ii,jj)=noiseimage(:,:,ii,jj).*(noiseimage(:,:,ii,jj)>0);
        end
    end
    widefield=widefield./max(widefield(:))*max(noiseimage(:));
    
    
    separated_FT=zeros(xsize,ysize,a_num,3);
    noiseimagef=zeros(size(noiseimage));
    for ii=1:a_num
        re0_temp=zeros(xsize,ysize);
        rep_temp=zeros(xsize,ysize);
        rem_temp=zeros(xsize,ysize);
        modulation_matrix=[1,1/2*exp(-1i*(pi*0)),1/2*exp(1i*(pi*0));...
            1,1/2*exp(-1i*(pi*2/3)),1/2*exp(1i*(pi*2/3));...
            1,1/2*exp(-1i*(pi*4/3)),1/2*exp(1i*(pi*4/3))];
        matrix_inv=inv(modulation_matrix);
    
        for jj=1:p_num
            noiseimagef(:,:,ii,jj)=fftshift(fft2(noiseimage(:,:,ii,jj)));
            re0_temp=matrix_inv(1,jj)*noiseimagef(:,:,ii,jj)+re0_temp;
            rep_temp=matrix_inv(2,jj)*noiseimagef(:,:,ii,jj)+rep_temp;
            rem_temp=matrix_inv(3,jj)*noiseimagef(:,:,ii,jj)+rem_temp;
        end
    
        separated_FT(:,:,ii,1)=re0_temp;
        separated_FT(:,:,ii,2)=rep_temp;
        separated_FT(:,:,ii,3)=rem_temp;
    end
    clear re0_temp rep_temp rem_temp noiseimage
    
    fmask=double(sqrt(xr.^2+yr.^2)>cutoff*mask_factor);
    [shiftvalue,~]=frequency_est_tirf_v2(separated_FT,0.008,fmask,show_initial_result_flag,mask_factor*cutoff);
    clear separated_FT
    
    
    for ii=1:a_num
        shiftvalue(ii,2,:)=shiftvalue(ii,2,:)-shiftvalue(ii,1,:);
        shiftvalue(ii,3,:)=shiftvalue(ii,3,:)-shiftvalue(ii,1,:);
        shiftvalue(ii,1,1)=0;
        shiftvalue(ii,1,2)=0;
    end
    
    %% phase correction with inverse matrix based algorithm
    search_range=0.4;%the max radius in the local search algorithm
    
    %obtain a more precise estimation of the period and the directon of sinusodial pattern
    [ precise_shift,~] = precise_frequency_tirf(noiseimagef,shiftvalue,search_range);
    
    
    
    %% to show you the final result
    para = zeros(a_num,2);  % angle, T
    for i = 1:a_num
        angle = 90 - atan( precise_shift(i,2,1)/(precise_shift(i,2,2)) )/pi*180;
        if flag~=2
            if i == 1
                if angle>90
                    angle = angle - 180;
                end
            end
        end


        angle = angle
        kmin = 2*pi/psize/xsize;
        k_shift = kmin*(precise_shift(i,2,1)*precise_shift(i,2,1)+precise_shift(i,2,2)*precise_shift(i,2,2))^(0.5);
        T = 2*pi/k_shift
        para(i,1) = angle;
        para(i,2) = T;
    end





end

