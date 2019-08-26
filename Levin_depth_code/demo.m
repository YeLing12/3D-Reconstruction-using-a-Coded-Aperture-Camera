function [x]=demo()
filt=0;
load demo_inp
I=im2double(imread('C:\Users\LY\Desktop\test_image\result\binary\22_bluring_far.JPG'));
% I = im2double(rgb2gray(I));
I = I(:,:,1);
figure, imshow([I]);
title('input')
drawnow

% 
imresize(filt/max(filt(:)),20);
figure, imshow(imresize(filt/max(filt(:)),20))
title('kernel')
drawnow


[sdI1]=deconvSps(I,filt,0.001,200);
[sdI2]=deconvSps(I,filt,0.004,200);
x=sdI1;

figure, imshow([sdI1,sdI2])
title('sparse deconv (varying smoothness weight)')
drawnow



