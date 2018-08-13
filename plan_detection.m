function [ segmented_plan,varargout ] = plan_detection( img,voting_space,foe,flow,varargin )
% INPUT: voting_space
%        foe focus of expansion y or x depend on voting space
%        flow the corresponding optical flow component
% OUTPUT: segmented_plan H x W matrix 1 for the parabola plan, 2 for the
% line plan
%
global scale;
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
   % voting_space_rethresh(voting_space_rethresh < W*0.01)=0;
   % voting_space_rethresh(voting_space_rethresh > 0)=1;
    figure(11)
    subplot(1,2,1)
    imshow(voting_space_rethresh,[]);
    % Detect parabola for half dbtm of voting space
    h_down = voting_space_rethresh(foe:H_v,:);
    a = para_detection2(h_down,0.0002);
   % a= 0.0035;
    varargout{1} = a;
    
    drawed = draw_parabola(voting_space_rethresh,a,foe,1);
    [segmented_plan,projected] = re_project(img,flow,foe,a,'parabola',1,1,segmented_plan);
    
%     h_up = voting_space_rethresh(1:foe,:);
%     a = para_detection(h_up,0.002)
%     drawed = draw_parabola(voting_space_rethresh,a,foe,0,drawed);
%     [segmented_plan,projected] = re_project(img,flow,foe,a,'parabola',1,0,segmented_plan);

    subplot(1,2,2);
    imshow(drawed);
%     figure(13);
%     imshow(projected);
else
    % ux voting space
  %  voting_space_rethresh(voting_space_rethresh < H*0.01)=0;
  %  voting_space_rethresh(voting_space_rethresh > 0)=1;
    figure(21)
    subplot(1,2,1)
    imshow(voting_space_rethresh,[]);
    
    h_down = voting_space_rethresh(foe:H_v,:);
    a = para_detection2(h_down,0.00002);
    varargout{1} = a;
    drawed = draw_parabola(voting_space_rethresh,a,foe,1);
    [segmented_plan,~] = re_project(img,flow,foe,a,'parabola',0,1,segmented_plan);
    
    h_up = voting_space_rethresh(1:foe,:);
    a = para_detection2(flip(h_up),0.00002);
    varargout{2} = a;
    drawed = draw_parabola(voting_space_rethresh,a,foe,0,drawed);    
    [segmented_plan,projected] = re_project(img,flow,foe,a,'parabola',0,0,segmented_plan);
    
    subplot(1,2,2)
    imshow(drawed);
    figure(24);
    imshow(projected);
end
end

function a = para_detection(half_voting_space,sigma)
H = size(half_voting_space,1);
W = size(half_voting_space,2);
nb = sum(sum(half_voting_space> 0));
mat_a = zeros(1,nb);
% Calculate of all the alpha of y = alpha*x^2
k = 1;
for i =1:H
    for j=1:W
        if (half_voting_space(i,j) ~= 0)
            mat_a(k)= (j)/(i^2);
            k=k+1;
        end
    end
end

% Sort ascending
mat_a =sort(mat_a); %the most costly operation

% Finding the alpha where have the most points in a range sigma
%sigma =0.0002;
i=1;
low= mat_a(i);
high = low + sigma;
j=i+1;
while((mat_a(j)  < high) && (j<nb))
    j=j+1;
    
end
weight = (j-i)+1;
a = mean(mat_a(i:j));
limit = 0; %observing variable to end iteration early
for i =2:nb
    low = mat_a(i);
    high = low +sigma;
    k=j;
    while( mat_a(k) < high && (k<nb))
        if (k == nb)
            limit = 1;
            break;
        end
        k=k+1;
        
    end
    j=k;
    new_w = j-i+1;
    if (new_w > weight)
        weight = new_w;
        a = mean(mat_a(i:j));
    end
    if (limit)
        break;
    end
end
% figure(41)
% h=histogram(mat_a,(0:0.0002: 0.005));
% figure(42)
% histogram(mat_b,(0:0.05: 1));
end

function a = para_detection2(half_voting_space,sigma)
%Detect the parameter a with weighted voting space

H = size(half_voting_space,1);
W = size(half_voting_space,2);
nb = sum(sum(half_voting_space> 0));
mat_a = zeros(nb,2);
% Calculate of all the alpha of y = alpha*x^2
k = 1;
for i =1:H
    for j=1:W
        if (half_voting_space(i,j) ~= 0)
            mat_a(k,1)= (j)/(i^2);
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
% figure(42)
% histogram(mat_b,(0:0.05: 1));
end
function drawed = draw_parabola(voting_space,a,foe,down,varargin)
if nargin > 4
    drawed = varargin{1};
else
    drawed = double(repmat(voting_space,1,1,3));
end
H = size(voting_space,1);
W = size(voting_space,2);
R=cat(3,1, 0, 0);
if (down)
    for i=foe:0.5:H
        for j=1:0.5:W
            index = round(abs(a*(i-foe)^2))+1;
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
function [mat_seg,projected]=re_project(img,flow,foe,a,type,vy,down,mat_seg,varargin)
H  = size(img,1);
W = size(img,2);
global scale;
if (mat_seg ==0)
    mat_seg =zeros(H,W,2);
end
low = 0.85;
high = 1.15;
thresh = 0.25;
switch type
    case 'line'
        
    case 'parabola'
        if (vy) % vy voting space, horizontal plane
            if (down) %half bottom
                for i=foe:H
                    v_ref = abs(a*(i-foe)^2-1)/scale;
                    for j =1:W
                        if ( ((v_ref*low<=abs(flow(i,j)) ) && (v_ref*high >= abs(flow(i,j)) ))...
                                || ( (v_ref-thresh<= abs(flow(i,j))) && (v_ref+thresh>=abs(flow(i,j))) )...
                                && (v_ref >0.5))
                            mat_seg(i,j,1)=1;
                        end
                    end
                end
            else %half top
                for i=1:foe
                    v_ref = abs(a*(i-foe+1)^2-1)/scale;
                    for j =1:W
                        if ( ((v_ref*low<=abs(flow(i,j)) ) && (v_ref*high >= abs(flow(i,j)) ))...
                                || ( (v_ref-thresh<= abs(flow(i,j))) && (v_ref+thresh>=abs(flow(i,j))) )...
                                && (v_ref >0.5))
                            mat_seg(i,j,1)=1;
                        end
                    end
                end                
            end
        else %ux voting space, lateral plane
            if (down) %half bottom
                for i=foe:W
                    u_ref = abs(a*(i-foe+1)^2-1)/scale;
                    for j =1:H
                        if ( ((u_ref*low<=abs(flow(j,i)) ) && (u_ref*high >= abs(flow(j,i)) ))...
                                || ( (u_ref-thresh<= abs(flow(j,i))) && (u_ref+thresh>=abs(flow(j,i))) )...
                                && (u_ref >0.5))
                            mat_seg(j,i,2)=1;
                        end
                    end
                end
            else %half top
                for i=1:foe
                    u_ref = abs(a*(foe-i+1)^2-1)/scale;
                    for j =1:H
                        if ( ((u_ref*low<=abs(flow(j,i)) ) && (u_ref*high >= abs(flow(j,i)) ))...
                                || ( (u_ref-thresh<= abs(flow(j,i))) && (u_ref+thresh>=abs(flow(j,i))) )...
                                && (u_ref >0.5))
                            mat_seg(j,i,2)=1;
                        end
                    end
                end
            end
        end
end
projected = color_plan(mat_seg,img);
end
function colored_map = color_plan(mat_seg,img)
colored_map = img;
H= size(img,1);
W= size(img,2);
R=cat(3,1, 0, 0);
G=cat(3,0, 1, 0);
B=cat(3,0, 0, 1);
for i =1:H
    for j=1:W
        coeff = 1/(sum(mat_seg(i,j,:))+1);
        if (j < W/2)
            colored_map(i,j,:)=colored_map(i,j,:)*coeff+mat_seg(i,j,1)*R*coeff+...
                mat_seg(i,j,2)*G*coeff;
        else
            colored_map(i,j,:)=colored_map(i,j,:)*coeff+mat_seg(i,j,1)*R*coeff+...
                mat_seg(i,j,2)*B*coeff;
        end
    end
end

end

