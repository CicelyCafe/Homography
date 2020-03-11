function practical1B

%the aim of the second part of practical 1 is to use the homography routine
%that you established in the first part of the practical.  We are going to
%make a panorama of several images that are related by a homography.  I
%provide 3 images (one of which is has a large surrounding region) and a
%matching set of points between these images.

%close all open figures
close all;

%load in the required data
load('PracticalDataSm','im1','im2','im3','pts1','pts2','pts3','pts1b');
%im1 is center image with grey background
%im2 is left image 
%pts1 and pts2 are matching points between image1 and image2
%im3 is right image
%pts1b and pts3 are matching points between image 1 and image 3

%show images and points
%figure; set(gcf,'Color',[1 1 1]);image(uint8(im1));axis off;hold on;axis image;
%plot(pts1(1,:),pts1(2,:),'r.'); 
%plot(pts1b(1,:),pts1b(2,:),'m.');
%figure; set(gcf,'Color',[1 1 1]);image(uint8(im2));axis off;hold on;axis image;
%plot(pts2(1,:),pts2(2,:),'r.'); 
%figure; set(gcf,'Color',[1 1 1]);image(uint8(im3));axis off;hold on;axis image;
%plot(pts3(1,:),pts3(2,:),'m.'); 

%****TO DO**** 
%calculate homography from pts1 to pts2
Homo_2 = calcBestHomography(pts1, pts2);
Homo_3 = calcBestHomography(pts1b,pts3);

[m1 n1 z] = size(im1);
[m2 n2 z] = size(im2);
[m3 n3 z] = size(im3);

for p = 1:m1*n1
    v = ceil(p/n1);
    u = p-(v-1)*n1;
    W = [u;v;1];
    X = Homo_2*W;
    Xcart = X(1:2,:)./repmat(X(3,:),2,1);
    if Xcart(2)>=0 && Xcart(2)<= m2 
        if Xcart(1)>=0 && Xcart(1)<=n2
            im1(v,u,:) = im2(ceil(Xcart(2)),ceil(Xcart(1)),:);
        end
    end
    %repeat the above process mapping image 3 to image 1.
    X3 = Homo_3*W;
    Xcart3 = X3(1:2,:)./repmat(X3(3,:),2,1);
    if Xcart3(2)>=0 && Xcart3(2)<= m3 
        if Xcart3(1)>=0 && Xcart3(1)<=n3
            im1(v,u,:) = im3(ceil(Xcart3(2)),ceil(Xcart3(1)),:);
        end
    end
end


figure; set(gcf,'Color',[1 1 1]);image(uint8(im1));axis off;hold on;axis image;

%****TO DO**** 
%for every pixel in image 1


    %transform this pixel position with your homography to find where it 
    %is in the coordinates of image 2
    %if it the transformed position is within the boundary of image 2 then 
        %copy pixel colour from image 2 pixel to current position in image 1 
        %draw new image1 (use drawnow to force it to draw)
    %end
%end;

%****TO DO****
%repeat the above process mapping image 3 to image 1.


function H = calcBestHomography(pts1Cart, pts2Cart)



%should apply direct linear transform (DLT) algorithm to calculate best
%homography that maps the points in pts1Cart to their corresonding matchin in 
%pts2Cart

m = size(pts1Cart,2);
pts1CartHom = [pts1Cart; ones(1,m)];

A = zeros(m*2,9);
for i=1:m
    A((2*i-1),1:3)=0;
    A((2*i),4:6)=0;
    A((2*i-1),4:6)= -pts1CartHom(1:3,i);
    A((2*i),1:3)= pts1CartHom(1:3,i);
    A((2*i-1),7:9)= pts2Cart(2,i)*pts1CartHom(1:3,i);
    A((2*i),7:9)= -pts2Cart(1,i)*pts1CartHom(1:3,i); 
end    
%disp(A);
h = solveAXEqualsZero(A);
H = transpose(reshape(h,[3,3]));

function x = solveAXEqualsZero(A);
[U,S,V] = svd(A);
x = V(:,end);



