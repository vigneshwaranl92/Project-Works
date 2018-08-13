%rows and columns
H=500; W=500;
%Planes and flow 
plane = zeros(500,500);
plane2 = zeros(500,500);
flow = zeros(500,500,2);
%FOE Variables
y=0.3;
pixel_size = 1e-5;
f = 0.005;
%Loop to find the optical flow
for z=1:0.001:20
    for x = -0.5:0.001:0.5
        z_n = z-0.5;
        y_n = y+0.182;         
        i = (f/pixel_size)*(y/z);
        j = (f/pixel_size)*(x/z);
        %Convert back to matlab coordinates
        ii = round(H/2-i);
        jj = round(W/2-j);
        %change y to y_n to add the Ty motion
        i_n = (f/pixel_size)*(y_n/z_n);
        j_n = (f/pixel_size)*(x/z_n);
        
        %Convert back to matlab coordinate
        ii_n = round(H/2-i_n);
        jj_n = round(W/2-j_n);
        if (( ii > 0 ) && ( ii <= H ) && ( jj > 0 ) && (jj <= W))
            plane(ii,jj)=1;
            if(flow(ii,jj,1) && flow(ii,jj,1))
                flow(ii,jj,1)=(j-j_n)*0.5+flow(ii,jj,1)*0.5;
                flow(ii,jj,2)=(i-i_n)*0.5+flow(ii,jj,2)*0.5;
            else
                flow(ii,jj,1)=(j-j_n);
                flow(ii,jj,2)=(i-i_n);
            end
        end
        
        if (( ii_n > 0 ) && ( ii_n <= H ) && ( jj_n > 0 ) && (jj_n <= W))
            plane2(ii_n,jj_n)=1;
        end
    end
end
%optical flow figures
u1=flow(:,:,1);
v1=flow(:,:,2);
figure(1)
subplot(1,2,1)
imshow(plane);
subplot(1,2,2)
imshow(plane2);
%changing the direction of v1 to v(ie, +ve to -ve) also giving the scale value to give a better resolution
% Convert back to camera axis
scale=1;
v=(-v1).* scale;
%voting space variables
rows = size(v,1);
columns =size(v,2);
w=ones(rows,columns);
min_v=floor(min(v(:)));
max_v =ceil(max(v(:))); 
max_v_s=(((min_v).*(-1)) +1);

%creating the voting space 
vy = zeros(H,max_v_s);
vy1 = zeros(H,max_v_s);
%positiveIndexes = v >= 0;
%negativeIndexes = v <  0;
low_v_E=floor(v(:,:));
high_v_E=ceil(v(:,:));
indices_1_E= ((max_v_s + low_v_E));
indices_2_E = ((max_v_s + high_v_E));
index = (max_v_s + v);
 for i=1:H
	for j=1:W
		low_v =floor(v(i,j)) ;
		high_v=ceil(v(i,j));
		indices_1 = (max_v_s + low_v);	
		indices_2 = (max_v_s + high_v);	
       if (indices_1 ~= indices_2)
           vy(i,indices_1)=vy(i,indices_1)+w(i,j);
           vy(i,indices_2)=vy(i,indices_2)+w(i,j);
       else 
           vy(i,index(i,j))=vy(i,index(i,j))+w(i,j);
       end
	end
 end


figure(4)
imshow(vy,[])
figure(5)
imshow(vy)

