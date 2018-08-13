pixel_size = 1e-5;
f = 0.005;
%plan_lateralA=[0.5 0.5 1;0.5 0.5 2;0.5 -0.5 1 ; 0.5 -0.5 2];
H=500; W=500;
image_plan = zeros(500,500);
image_plan2 = zeros(500,500);
image_flow = zeros(500,500,2);
y=0.3;
for z=1:0.001:20
    for x = -0.5:0.001:0.5
        z_n = z-0.5;
        y_n = y+0.182; 
        a= 10*pi/180;
        x_temp = y_n;
	
	%for rotation
        %x_n = cos(a)*x_n+sin(a)*z_n;
        %z_n = -sin(a)*x_temp + cos(a)*z_n;
        
        i = (f/pixel_size)*(x/z);
        j = (f/pixel_size)*(y/z);
        ii = round(H/2-i);
        jj = round(W/2-j);
        
        i_n = (f/pixel_size)*(x/z_n);
        j_n = (f/pixel_size)*(y_n/z_n);
        ii_n = round(H/2-i_n);
        jj_n = round(W/2-j_n);
        if (( ii > 0 ) && ( ii <= H ) && ( jj > 0 ) && (jj <= W))
            image_plan(ii,jj)=1;
            if(image_flow(ii,jj,1) && image_flow(ii,jj,1))
                image_flow(ii,jj,1)=(j-j_n)*0.5+image_flow(ii,jj,1)*0.5;
                image_flow(ii,jj,2)=(i-i_n)*0.5+image_flow(ii,jj,1)*0.5;
            else
                image_flow(ii,jj,1)=(j-j_n);
                image_flow(ii,jj,2)=(i-i_n);
            end
        end
        
        if (( ii_n > 0 ) && ( ii_n <= H ) && ( jj_n > 0 ) && (jj_n <= W))
            image_plan2(ii_n,jj_n)=1;
        end
    end
end
figure(1)
subplot(1,3,1)
imshow(image_plan);
subplot(1,3,2)
imshow(image_plan2);
% figure(3)
% surf(image_flow(:,:,1))
subplot(1,3,3)
imshow(abs(image_plan2 - image_plan));
voting_space = zeros(1+round(max(max(abs(image_flow(:,:,1))))),W);
u=image_flow(:,:,1);
v=image_flow(:,:,2);
rak = round(abs(image_flow(:,:,1)))+1;
for i=1:H
    for j =1:W
        ii = round(abs(image_flow(i,j,1)))+1;
        voting_space(ii,j) = voting_space(ii,j)+1;
    end
end
figure(2)
imshow(voting_space,[]);