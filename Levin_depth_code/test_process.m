clear; close; clc;
load('FinalKernelResult.mat');
I=imread('\\surrey.ac.uk\personal\HS211\yl01690\.System\Desktop\YeLing\Combine\29.JPG');
I = I(300:480,:,:);
I_gray_double = im2double(I(:,:,1));
[ry, cy] = size(I_gray_double);
we = 0.001;
max_it = 200;

range = 1:8;
depthMapWindow = [9, 9];
lens_length_correction = -0.0209;% m
depths = (2.35:0.1:3.05) + lens_length_correction; % m
wingr = floor(depthMapWindow(1)/2);
wingc = floor(depthMapWindow(2)/2);
depthMap = zeros(ry, cy);
sharpImage=zeros(ry, cy);
for r = (wingr+1):(ry-wingr)
    for c = (wingc+1):(cy-wingc)
        % re-initialize
        minLocalEnergy = inf;
        minInd = 0;
        fprintf('Row: %i, Column: %i\n', r,c);
        % obtain local window of scene
        I_local = I_gray_double((r-wingr):(r+wingr), (c-wingc):(c+wingc));
        for i = range
            % rotate kernel to match convention of deconSps %%%
            %tempDeconv = deconvSps(I_local, rot90(PSFs_11{i}, 2), we, max_it);
            tempDeconv = deconvSps(I_local,PSFs_11{i} , we, max_it); %imrotate(filt,i*18) experient with rotation % Use the right orientation Kernel!!!!!!
            % tempDeconv = deconvL2(I_local, rot90(PSFs_11{i}, 2), we);

            % calculate reconstruction error for local window
            if size(PSFs_11{i}, 1) > size(tempDeconv, 1)
                reconv = conv2(PSFs_11{i}, tempDeconv, 'valid');
            else
                reconv = conv2(tempDeconv, PSFs_11{i}, 'valid');
            end
            reconv_wing = floor(size(reconv, 1)/2);
            reconError = I_local((wingr+1-reconv_wing):(wingr+1+reconv_wing), (wingc+1-reconv_wing):(wingc+1+reconv_wing)) - reconv;

            % calculate average local energy for this window (Frobenius
            % norm)
            % avgLocalEnergy = sum(sum(reconError.^2));
            avgLocalEnergy = mean(mean(reconError))^2;
            %avgLocalEnergy=I_local;
            % find minimum energy
            if avgLocalEnergy < minLocalEnergy
                minLocalEnergy = avgLocalEnergy; % update
                minInd = i;
            end
        end
        % update depthMap 
        depthMap(r, c) = depths(minInd);
        sharpImage(r,c)=I_gray_double(r, c);
%         disp(depths(minInd));
    end

end
% resize
depthMap = depthMap((wingr+1):(ry-wingr), (wingc+1):(cy-wingc));


% figure;
% imshow(I_cropped);
% figure; 
% imagesc(depthMap);
% axis equal;
% title('Depth Map');
% colorbar;