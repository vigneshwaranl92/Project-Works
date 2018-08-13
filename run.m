addpath(genpath('..'))
% addpath('..')
clear all;
% Load optical flow:
uv = readFlowFile('000000.flo');

u = abs(uv(:,:,1));
v = abs(uv(:,:,2));
low_v_E =floor(v(:,:));
high_v_E=ceil(v(:,:));

displayImg = imread('img1.png');
plotFlow(u, v, displayImg, 10, 25);
hold on
% img = imread('C:\Users\MAI\Dropbox\new code\Matlab\KITTI Test\000000_10.png');
img = double(displayImg)./255;
global scale;
scale = 1;
[ux,vy]=create_voting_space(uv,scale);
% seg_r = plan_detection(img,vy,187,uv(:,:,2));
seg_r = parabolaabc(img,vy,187,uv(:,:,2));
figure(2)
subplot(1,2,1)
imshow(vy)
subplot(1,2,2)
imshow(ux)
