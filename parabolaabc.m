function [ segmented_plan ] = plan_detection( img,voting_space,foe,flow,varargin )

if nargin > 4
    segmented_plan = varargin{1};
else
    segmented_plan =0;
end
H = size(flow,1); W =size(flow,2);
H_v = size(voting_space,1);

voting_space(:,1)=0;
voting_space_rethresh=voting_space;

if (H_v == H )
    % vy voting space
    figure(11)
    subplot(1,2,1)
    imshow(voting_space_rethresh,[]);
    % Detect parabola for half dbtm of voting space
	h_down = voting_space_rethresh(foe:H_v,:);
	%%%%%%%%%%%%%%%%%%%%%
	sigma = 0.0002;
	
	a=paraa(h_down,sigma)
	b=parab(h_down,sigma)
	c=parac(h_down,sigma)
    drawed = draw_parabola(voting_space_rethresh,a,b,c,foe,1);
	
    subplot(1,2,2);
    imshow(drawed);
else
    % ux voting space
    figure(21)
    subplot(1,2,1)
    imshow(voting_space_rethresh,[]);
    
    h_down = voting_space_rethresh(foe:H_v,:);
    a = para_detection2(h_down,0.00002);
    varargout{1} = a;
    % drawed = draw_parabola(voting_space_rethresh,a,foe,1);
    % [segmented_plan,~] = re_project(img,flow,foe,a,'parabola',0,1,segmented_plan);
    
    h_up = voting_space_rethresh(1:foe,:);
    a = para_detection2(flip(h_up),0.00002);
    varargout{2} = a;
    % drawed = draw_parabola(voting_space_rethresh,a,foe,0,drawed);    
    % [segmented_plan,projected] = re_project(img,flow,foe,a,'parabola',0,0,segmented_plan);
    
    subplot(1,2,2)
    imshow(drawed);
    figure(24);
    imshow(projected);
end
end


function a = paraa(half_voting_space,sigma)
%Detect the parameter a with weighted voting space

H = size(half_voting_space,1);
W = size(half_voting_space,2);
nb = sum(sum(half_voting_space> 0));
mat_a = zeros(nb,2);

k = 1;

x0 = half_voting_space(42,4);
y0 = half_voting_space(42,5);
x1 = half_voting_space(83,13);
y1 = half_voting_space(83,14);
x2 = half_voting_space(184,58);
y2 = half_voting_space(184,59);

numerator1=(y2 * (x1-x0) + y1 * (x0-x2) + y0 * (x2-x1));
denominator1=((y0-y1) * (y0-y2) * (y1-y2));
denominator2=((y0-y1) * (y0-y2) * (y1-y2));
denominator3=((y0-y1) * (y0-y2) * (y1-y2));
numerator2=((y2^2) * (x0-x1) + (y1^2) * (x2-x0) + (y0^2) * (x1-x2));
numerator3=((y1 *y2 * (y1-y2) * x0) + (y2 * y0* (y2-y0) * x1) + (y0 * y1 * (y0-y1) * x2));
for i =1:H
    for j=1:W
        if (half_voting_space(i,j) ~= 0) 
            mat_a(k,1)= (numerator1/denominator1);
            mat_a(k,2)= half_voting_space(i,j);
            k=k+1;
        end
    end
end

% Sort ascending
mat_a =sortrows(mat_a,1); %the most costly operation

% Finding the alpha where have the most points in a range sigma
%sigma =0.0002;
i=1;
low= mat_a(i,1);
weight = mat_a(i,2);
weight_low=weight;
high = low + sigma;
j=i+1;
while((mat_a(j,1)  < high) && (j<nb))
    weight=weight + mat_a(j,2);
    j=j+1;
end
weight_old=weight;
a = mean(mat_a(i:j,1));
%limit = 0; %observing variable to end iteration early
for i =2:nb
    low = mat_a(i,1);% lower bound of a
    weight = weight-weight_low; % Remove the lowest weight
    weight_low = mat_a(i,2); %Redefine the lowest weight
    high = low + sigma;%higher bound of a
    k=j;
    while( (mat_a(k,1) < high) && (k<nb))
        weight=weight + mat_a(k,2);%add the new weight
        k=k+1;
        
    end
    j=k;
    weight_new = weight;
    if (weight_new > weight_old)
        weight_old = weight_new;
        a = mean(mat_a(i:j,1));
    end
end
% figure(41)
% h=histogram(mat_a,(0:0.0002: 0.005));

end
function b = parab(half_voting_space,sigma)
%Detect the parameter b with weighted voting space

H = size(half_voting_space,1);
W = size(half_voting_space,2);
nb = sum(sum(half_voting_space> 0));
mat_b = zeros(nb,2);

k = 1;

x0 = half_voting_space(20,8);
y0 = half_voting_space(20,11);
x1 = half_voting_space(88,12);
y1 = half_voting_space(89,12);
x2 = half_voting_space(184,58);
y2 = half_voting_space(184,59);

numerator1=(y2 * (x1-x0) + y1 * (x0-x2) + y0 * (x2-x1));
denominator1=((y0-y1) * (y0-y2) * (y1-y2));
denominator2=((y0-y1) * (y0-y2) * (y1-y2));
denominator3=((y0-y1) * (y0-y2) * (y1-y2));
numerator2=((y2^2) * (x0-x1) + (y1^2) * (x2-x0) + (y0^2) * (x1-x2));
numerator3=((y1 *y2 * (y1-y2) * x0) + (y2 * y0* (y2-y0) * x1) + (y0 * y1 * (y0-y1) * x2));
for i =1:H
    for j=1:W
        if (half_voting_space(i,j) ~= 0) 
            mat_b(k,1)= (numerator2/denominator2);
            mat_b(k,2)= half_voting_space(i,j);
            k=k+1;
        end
    end
end

% Sort ascending
mat_b =sortrows(mat_b,1); %the most costly operation

% Finding the alpha where have the most points in a range sigma
%sigma =0.0002;
i=1;
low= mat_b(i,1);
weight = mat_b(i,2);
weight_low=weight;
high = low + sigma;
j=i+1;
while((mat_b(j,1)  < high) && (j<nb))
    weight=weight + mat_b(j,2);
    j=j+1;
end
weight_old=weight;
b = mean(mat_b(i:j,1));
%limit = 0; %observing variable to end iteration early
for i =2:nb
    low = mat_b(i,1);% lower bound of a
    weight = weight-weight_low; % Remove the lowest weight
    weight_low = mat_b(i,2); %Redefine the lowest weight
    high = low + sigma;%higher bound of a
    k=j;
    while( (mat_b(k,1) < high) && (k<nb))
        weight=weight + mat_b(k,2);%add the new weight
        k=k+1;
        
    end
    j=k;
    weight_new = weight;
    if (weight_new > weight_old)
        weight_old = weight_new;
        b = mean(mat_b(i:j,1));
    end
end

% figure(42)
% histogram(mat_b,(0:0.05: 1));
end

function c = parac(half_voting_space,sigma)
%Detect the parameter c with weighted voting space

H = size(half_voting_space,1);
W = size(half_voting_space,2);
nb = sum(sum(half_voting_space> 0));
mat_c = zeros(nb,2);

k = 1;

x0 = half_voting_space(20,8);
y0 = half_voting_space(20,11);
x1 = half_voting_space(83,13);
y1 = half_voting_space(83,14);
x2 = half_voting_space(184,58);
y2 = half_voting_space(184,59);

numerator1=(y2 * (x1-x0) + y1 * (x0-x2) + y0 * (x2-x1));
denominator1=((y0-y1) * (y0-y2) * (y1-y2));
denominator2=((y0-y1) * (y0-y2) * (y1-y2));
denominator3=((y0-y1) * (y0-y2) * (y1-y2));
numerator2=((y2^2) * (x0-x1) + (y1^2) * (x2-x0) + (y0^2) * (x1-x2));
numerator3=((y1 *y2 * (y1-y2) * x0) + (y2 * y0* (y2-y0) * x1) + (y0 * y1 * (y0-y1) * x2));
for i =1:H
    for j=1:W
        if (half_voting_space(i,j) ~= 0) 
            mat_c(k,1)= (numerator3/denominator3);
            mat_c(k,2)= half_voting_space(i,j);
            k=k+1;
        end
    end
end

% Sort ascending
mat_c =sortrows(mat_c,1); %the most costly operation

% Finding the alpha where have the most points in a range sigma
%sigma =0.0002;
i=1;
low= mat_c(i,1);
weight = mat_c(i,2);
weight_low=weight;
high = low + sigma;
j=i+1;
while((mat_c(j,1)  < high) && (j<nb))
    weight=weight + mat_c(j,2);
    j=j+1;
end
weight_old=weight;
c = mean(mat_c(i:j,1));
%limit = 0; %observing variable to end iteration early
for i =2:nb
    low = mat_c(i,1);% lower bound of a
    weight = weight-weight_low; % Remove the lowest weight
    weight_low = mat_c(i,2); %Redefine the lowest weight
    high = low + sigma;%higher bound of a
    k=j;
    while( (mat_c(k,1) < high) && (k<nb))
        weight=weight + mat_c(k,2);%add the new weight
        k=k+1;
        
    end
    j=k;
    weight_new = weight;
    if (weight_new > weight_old)
        weight_old = weight_new;
        c = mean(mat_c(i:j,1));
    end
end

% figure(42)
% histogram(mat_c,(0:0.05: 1));
end





function drawed = draw_parabola(voting_space,a,b,c,foe,down,varargin)

if nargin > 6
    drawed = varargin{1};
else
    drawed = double(repmat(voting_space,1,1,3));
end
H = size(voting_space,1);
W = size(voting_space,2);
R=cat(3,0, 1, 0);
if (down)
    for i=foe:0.5:H
        for j=1:0.5:W
            index = round(abs(a*(i-foe)^2) + abs(b*(i-foe)) + abs(c))+1;
			% index = round(abs(a*(i-foe)^2) + 0.4*(i-foe) + 0.4)+1;
            if ((index > 0) && (index <= W))
            drawed(round(i),index,:)= R;
            end
        end
    end
else
    for i=1:0.5:foe
        for j=1:0.5:W
            index = round(abs(a*(foe-i)^2))+1;
            drawed(round(i),index,:)= R;
        end
    end
end

end
