function [Para1,Para2] = ransac_demo(data,num,threshDist,inlierRatio)

% data = given the h_down space
% num = number of values required while processing randomperm
%Threshold distance 
% According to qiong nie book,
% threshDist = num of inliers in data /num of points in data, hence keeping the rough value
threshDist=200;

inlierRatio=0.1;
num=5;
 
 figure;plot(data(:,:),'o');hold on;
 R_DIM = size(data,1); % Row Dimension
 C_DIM = size(data,2); % Column Dimension 
 bestInNum = 0; % Best fitting line with largest number of inliers
 bestParameter1=0;bestParameter2=0; % parameters for best fitting line
 
 for i=1:100
 %% Selecting Random points
     SRN = randperm(C_DIM,num);
	 sample = data(:,SRN);  
	 
 %% Computing the distances between all points with the fitting line 
     kLine = sample(:,2)-sample(:,1);% two points relative distance
     kLineNorm = kLine/norm(kLine);
     % normVector = [-kLineNorm(2),kLineNorm(1)];%Ax+By+C=0 A=-kLineNorm(2),B=kLineNorm(1)
	 normVector = [kLineNorm']
     distance = ((normVector)*(data - repmat(sample(:,1),1,C_DIM)));
	 
 %% Compute the inliers with distances smaller than the threshold
     inlierIdx = find(abs(distance)<=threshDist);
     inlierNum = length(inlierIdx);
	 
 %% Update the number of inliers and fitting model if better model is found     
     if inlierNum>=round(inlierRatio*C_DIM) && inlierNum>bestInNum
         bestInNum = inlierNum;
         bparameter1 = (sample(2,2)-sample(2,1))/(sample(1,2)-sample(1,1));
         bparameter2 = sample(2,1)-parameter1*sample(1,1);
         Para1=bparameter1; Para2=bparameter2;
     end
 end
 
 %% Plot the fitting line
 xAxis = -C_DIM/2:C_DIM/2; 
 yAxis = bestParameter1*xAxis + bestParameter2;
 plot(xAxis,yAxis,'r-','LineWidth',2);
