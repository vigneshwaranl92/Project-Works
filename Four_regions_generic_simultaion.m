%rows and columns
H=1000; W=1000;
%Planes and flow 
plane = zeros(1000,1000);
plane2 = zeros(1000,1000);
flow = zeros(1000,1000,2);
%FOE Variables
y= 0.05; %varies depending on the Yfoe value
pixel_size = 1e-5;
f = 0.005; %focal length
%Translation Parameters
Tz = -0.005; 
Ty = 0.1;
%Loop to find the optical flow
for z=1:0.001:20 %depth 
    for x = -0.5:0.001:0.5 %4 meters 
        z_n = z+Tz;
        y_n = y+Ty;         
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

figure(1)
subplot(1,3,1)
imshow(plane);
subplot(1,3,2)
imshow(plane2);
%Focus of expansion
Yfoe_P = (f) .* (Ty/Tz);
%Yfoe_N = -((f) .* (0.182/(-0.5)));
Yfoe_N = -(f) .* (Ty/Tz);
%optical flow figures
u1=flow(:,:,1);
v1=flow(:,:,2);
%changing the direction of v1 to v(ie, +ve to -ve) also giving the scale value to give a better resolution
% Convert back to camera axis
scale=1;
v=(-v1).* scale;
%voting space variables
rows = size(v,1);
columns =size(v,2);
w=ones(rows,columns);
min_v=floor(min(v(:)));
min_v_s=(((min_v).*(-1)) +1);
max_v_s = (ceil(max(v(:)))+1);

%creating the voting space 
vy_N = zeros(rows,min_v_s);
vy_P = zeros(rows,max_v_s);
vy_NP = zeros(rows, min_v_s + max_v_s);
%Postive and Negative entries of the v values
positiveIndexes = v > 0;
negativeIndexes = v <= 0;
Negative = any(v(:)<0);
Positive = any(v(:)>0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
low_v_E=floor(v(:,:));
high_v_E=ceil(v(:,:));
indices_1_E= ((min_v_s + low_v_E));
indices_2_E = ((min_v_s + high_v_E));
index1 = (min_v_s + v);
index2 = (min_v_s + v);
if ((Positive ==1) && (Negative==1))
  for i=1:rows
     for j=1:columns
		low_v =floor(v(i,j)) ;
		high_v=ceil(v(i,j));
		indices_1 = (min_v_s + low_v);	
		indices_2 = (min_v_s + high_v);	
	   if (indices_1 ~= indices_2)
			vy_NP(i,indices_1)=vy_NP(i,indices_1)+w(i,j);
			vy_NP(i,indices_2)=vy_NP(i,indices_2)+w(i,j);
       else 
			vy_NP(i,index1(i,j))=vy_NP(i,index1(i,j)) + w(i,j);
       end
	end
  end
figure(2)
imshow(vy_NP,[])

elseif ((Positive ==1) && (Negative==0))
  for i=1:rows
	for j=1:columns
        low_v =floor(v(i,j)) ;
		high_v=ceil(v(i,j));	
       if (low_v ~= high_v)
           vy_P(i,low_v+1)=vy_P(i,low_v+1)+(v(i,j)-low_v)*w(i,j);
           vy_P(i,high_v+1)=vy_P(i,low_v+1)+(high_v-v(i,j))*w(i,j);
       else 
           vy_P(i,v(i,j)+1)=vy_P(i,v(i,j)+1) + w(i,j);
       end
	end
  end
figure(3)
imshow(vy_P)
 
else ((Positive ==0) && (Negative==1))
  for i=1:rows
	for j=1:columns
		low_v =floor(v(i,j)) ;
		high_v=ceil(v(i,j));
		indices_1 = (min_v_s + low_v);	
		indices_2 = (min_v_s + high_v);	
       if (indices_1 ~= indices_2)
           vy_N(i,indices_1)=vy_N(i,indices_1)+w(i,j);
           vy_N(i,indices_2)=vy_N(i,indices_2)+w(i,j);
       else 
           vy_N(i,index1(i,j))=vy_N(i,index1(i,j)) + w(i,j);
       end
	end
  end
figure(4)
imshow(vy_N,[])
figure(4)
imshow(vy_N)
end