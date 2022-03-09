
% ECE4076 Lab 4 
% Written by : Zhiyue Li
% Last Modified Date: 
close all;clc;clear all;

% Task 1 
% Load the images
img1 = imread('im1.jpg');
img2 = imread('im2.jpg');
img3 = imread('im3.jpg');
img4 = imread('im4.jpg');

double_wide12 = [img1,img2];
double_wide13 = [img1,img3];
double_wide14 = [img1,img4];

% Load the keypoints
imKp1 =load('im1.sift');
imKp2 =load('im2.sift');
imKp3 =load('im3.sift');
imKp4 =load('im4.sift');

% figure; imshow(img1);
% figure; imshow(img2);
% figure; imshow(img3);
% figure; imshow(img4);

% Tasks 2 and 3: Read and show SIFT keypoint locations
% Get feature points as a N*2 matrix
featurePts12 = [imKp1(:,1:2); imKp2(:,1:2)+[size(img1,2) 0]];
featurePts13 = [imKp1(:,1:2); imKp3(:,1:2)+[size(img1,2) 0]];
featurePts14 = [imKp1(:,1:2); imKp4(:,1:2)+[size(img1,2) 0]];

% Mark the joined image
imMarkedJoin12 = insertMarker(double_wide12,featurePts12,'x','color','r','size',5);
imMarkedJoin13 = insertMarker(double_wide13,featurePts13,'x','color','r','size',5);
imMarkedJoin14 = insertMarker(double_wide14,featurePts14,'x','color','r','size',5);

% Displays the marked joined image
% figure; imshow(imMarkedJoin12);
% figure; imshow(imMarkedJoin13);
% figure; imshow(imMarkedJoin14);


% Tasks 4 and 5: Match SIFT keypoints and show matches
Kp1_length = length(imKp1);
RightIndex1 = zeros(Kp1_length,1);

% Loop for all keypoints on the left image
for k = 1:Kp1_length

    % Euclidean distance
     dist = sqrt(sum((imKp1(k,5:end) - imKp2(:,5:end)).^2,2));
     
     % Find the nearest and 2nd nearest match (distanc d1 & d2)
     [dist, idx] = sort(dist,'ascend');

     % Calculate the ratio and distance
     if (dist(1)/dist(2) < 0.5)
         RightIndex1(k) = idx(1);
     end
%    RightIndex1(k) = idx(1);

end

% Find out all RightIndex1 >0
Idx = find(RightIndex1>0);

% Get the coordinates that match in both images.
pointsMatchL = imKp1(Idx,1:2);
pointsMatchR = imKp2(RightIndex1(Idx),1:2);

% Mark the points and draw lines
matchMark = insertMarker(imMarkedJoin12,[pointsMatchL; pointsMatchR + [size(img1,2) 0] ],'x','color','g','size',5);
figure;imshow(matchMark);

hold on;
for i = 1:length(Idx)
     line([pointsMatchL(i,1)  pointsMatchR(i,1)+ size(img1,2) ],[pointsMatchL(i,2)  pointsMatchR(i,2) ],'color','g');
end
hold off; 
title("Original image 1 and image 2 joined together");


% Task 6 Stereo reconstruction
flength = 1;
baseline = 1;
disp12 = (pointsMatchL(:,1) - pointsMatchR(:,1)); % disparity
depth12 = flength *baseline./abs(disp12);

m1 = mean(depth12);
d1 = std(depth12);
goodIndex = (depth12 > m1-d1) & (depth12 < m1+d1);

D12 = mean(disp12(goodIndex));

x12 = depth12.*[pointsMatchL ones(length(pointsMatchL),1)]; % 3d- coordinates of the points in image1's frame


% Repeat process for (image 1, image 3) and (image 1, image 4)
% Tasks 4 and 5: Match SIFT keypoints and show matches (image 1 and image
% 3)
Kp1_length = length(imKp1);
RightIndex13 = zeros(Kp1_length,1);

% Loop for all keypoints on the left image
for k = 1:Kp1_length

    % Euclidean distance
     dist = sqrt(sum((imKp1(k,5:end) - imKp3(:,5:end)).^2,2));
     [dist, idx] = sort(dist,'ascend');

     % Compare 2 smallest distances
     if (dist(1)/dist(2) < 0.5)
         RightIndex13(k) = idx(1);
      end

end

% Find out all RightIndex1 >0
Idx13 = find(RightIndex13>0);

% Get the coordinates that match in both images.
pointsMatchL13 = imKp1(Idx13,1:2);
pointsMatchR13 = imKp3(RightIndex13(Idx13),1:2);

% Mark the points and draw lines
matchMark13 = insertMarker(imMarkedJoin13,[pointsMatchL13; pointsMatchR13 + [size(img3,2) 0] ],'x','color','g','size',5);
figure;imshow(matchMark13);

hold on;
for i = 1:length(Idx13)
     line([pointsMatchL13(i,1)  pointsMatchR13(i,1)+ size(img1,2) ],[pointsMatchL13(i,2)  pointsMatchR13(i,2) ],'color','g');
end

hold off; 
title("Original image 1 and image 3 joined together");

disp13 = (pointsMatchL13(:,1) - pointsMatchR13(:,1)); % disparity
depth13 = flength *baseline./disp13;

m2 = mean(depth13);
d2 = std(depth13);
goodIndex = (depth13 > m2-d2) & (depth13 < m2+d2);

D13 = mean(disp13(goodIndex));
x13 = depth13.*[pointsMatchL13 ones(length(pointsMatchL13),1)];

% Repeat process for (image 1, image 3) and (image 1, image 4)
% Tasks 4 and 5: Match SIFT keypoints and show matches (image 1 and image
% 4)
Kp1_length = length(imKp1);
RightIndex14 = zeros(Kp1_length,1);

% Loop for all keypoints on the left image
for k = 1:Kp1_length

    % Euclidean distance
     dist = sqrt(sum((imKp1(k,5:end) - imKp4(:,5:end)).^2,2));
     [dist, idx] = sort(dist,'ascend');

     % Compare 2 smallest distances
     if (dist(1)/dist(2) < 0.5)
         RightIndex14(k) = idx(1);
      end

end

% Find out all RightIndex1 >0
Idx14 = find(RightIndex14>0);

% Get the coordinates that match in both images.
pointsMatchL14 = imKp1(Idx14,1:2);
pointsMatchR14 = imKp4(RightIndex14(Idx14),1:2);

% Mark the points and draw lines
matchMark14 = insertMarker(imMarkedJoin14,[pointsMatchL14; pointsMatchR14 + [size(img4,2) 0] ],'x','color','g','size',5);
figure;imshow(matchMark14);

hold on;
for i = 1:length(Idx14)
     line([pointsMatchL14(i,1)  pointsMatchR14(i,1)+ size(img1,2) ],[pointsMatchL14(i,2)  pointsMatchR14(i,2) ],'color','g');
end
hold off; 
title("Original image 1 and image 4 joined together");

disp14 = (pointsMatchL14(:,1) - pointsMatchR14(:,1));  % disparity
depth14 = flength *baseline./disp14;

m3 = mean(depth14);
d3 = std(depth14);
goodIndex = (depth14 > m3-d3) & (depth14 < m3+d3);

D14 = mean(disp14(goodIndex));

x14 = depth14.*[pointsMatchL14 ones(length(pointsMatchL14),1)]; % 3d- coordinates of the points in image1's frame

x = [x12;x13;x14];
figure;
scatter3(x12(:,1),x12(:,2),x12(:,3),'r');
hold on
scatter3(x13(:,1),x13(:,2),x13(:,3),'b');
hold on
scatter3(x14(:,1),x14(:,2),x14(:,3),'g');
hold off
legend("image 1 and image 2","image 1 and image 3","image 1 and image 4")


% Estimated baselines
D = abs([D12,D13,D14]);
boxwidth = 0.37;
baseLine = boxwidth.*D/479;
fprintf('Image 1 and Image 2 baseline:  %.4f \n',baseLine(1));
fprintf('Image 1 and Image 3 baseline:  %.4f \n',baseLine(2));
fprintf('Image 1 and Image 4 baseline:  %.4f \n',baseLine(3));



 %% Task 7 Reprojection of 3D points
disparity = (baseLine(2)-baseLine(3))/2/baseLine(3)*disp14;         
points = pointsMatchL14 + [disparity zeros(length(pointsMatchL14),1)];
img3Reconstruction = insertMarker(img3,points,'x','color','r','size',5);
img1Marked = insertMarker(img1,pointsMatchL14,'x','color','r','size',5);

figure;
subplot(1,2,1)
imshow(img1Marked);
title('Points from image 1');
subplot(1,2,2);
imshow(img3Reconstruction);
title('Reprojected points on Image 3');



