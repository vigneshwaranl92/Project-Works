function [ ux,vy ] = create_voting_space( uv,varargin )
% Create the voting space u-x, v-y.
% INPUT: uv optical flow
% OUTPUT: ux,vy the two voting space

H = size(uv,1);
W = size(uv,2);

% w variable to modify the voting based on its weight
if nargin > 2
    w = varargin{2};
else %default weight
    w=ones(H,W);
end

% Scale variable to improve the visualization of voting space
if nargin > 1
    scale = varargin{1};
else %default value
    scale=1;
end
% Initialize the parameters and variables
u = abs(uv(:,:,1))*scale;
v = abs(uv(:,:,2))*scale;
max_u = ceil(max(u(:)));
max_v = ceil(max(v(:)));
% Output matrix
ux = zeros(W,max_u+1);
vy = zeros(H,max_v+1);
low_v_E =floor(v(:,:));
       high_v_E=ceil(v(:,:));  
% Create voting space:
for i =1:H
    for j=1:W
       low_v =floor(v(i,j)) ;
       high_v=ceil(v(i,j));  
       if (low_v ~= high_v)
           vy(i,low_v+1)=vy(i,low_v+1)+(v(i,j)-low_v)*w(i,j);
           vy(i,high_v+1)=vy(i,low_v+1)+(high_v-v(i,j))*w(i,j);
       else %v component is integer
           vy(i,v(i,j)+1)=vy(i,v(i,j)+1) + w(i,j);
       end
       
       low_u =floor(u(i,j)); high_u=ceil(u(i,j));
       if (low_u ~= high_u)
           ux(j,low_u+1)=ux(j,low_u+1)+(u(i,j)-low_u)*w(i,j);
           ux(j,high_u+1)=ux(j,high_u+1)+(high_u-u(i,j))*w(i,j);
       else %u component is integer
           ux(i,u(i,j)+1)=ux(i,u(i,j)+1) + w(i,j);
       end
    end
end
end

