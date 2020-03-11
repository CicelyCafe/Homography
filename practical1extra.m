function practical1extra

%The aim of practical 1 is to calculate the homography that best maps two
%sets of points to one another.  We will (eventually) use this for creating
%panoramas, and for calculating the 3d pose of planes.  You should use this
%template for your code and fill in the missing sections marked "TO DO"

%close all open figures
close all;

%set of two dimensional Cartesian points
pts1Cart = [  240.5000   16.8351   33.5890  164.2696  149.1911;...
              248.8770  193.5890   251.3901 168.4581  228.7723];

%turn points to homogeneous representation
pts1Hom = [pts1Cart; ones(1,size(pts1Cart,2))];
      
%define a homography
H = [0.6 0.7 -100; 1.0 0.6 50; 0.001 0.002 1.0]

%apply homography to points
pts2Hom = H*pts1Hom;

%convert back to Cartesian
pts2Cart = pts2Hom(1:2,:)./repmat(pts2Hom(3,:),2,1)

%add a small amount of noise
noiseLevel = 4.0;
pts2Cart = pts2Cart+noiseLevel*randn(size(pts2Cart));

%draw two set of two dimensional points
%opens figure
figure; set(gcf,'Color',[1 1 1]);
%draw lines between each pair of points
nPoint = size(pts1Cart,2)
for (cPoint = 1:nPoint)
    %plot a green line between each pair of points
    plot([pts1Cart(1,cPoint) pts2Cart(1,cPoint)],[pts1Cart(2,cPoint) pts2Cart(2,cPoint)],'g-');
    %make sure we don't replace with next point
    hold on;
end;

%draws first set of points
plot(pts1Cart(1,:),pts1Cart(2,:),'b.','MarkerSize',20);
%remove axis
set(gca,'Box','Off');

%draws second set of points
plot(pts2Cart(1,:),pts2Cart(2,:),'r.','MarkerSize',20);

%now our goal is to transform the first points so that they map to the
%second set of points

%****TO DO****: Fill in the details of this routine 
%At the moment, it just returns and identity matrix (body is below)
%pts2EstCart = affine_transformation(pts1Cart, pts2Cart);
[R,t] = affine_transformation(pts1Cart, pts2Cart)

m = size(pts1Cart,2);
pts2EstCart = zeros(2,m);
for i = 1:m
    pts2EstCart(:,i) = R*pts1Cart(:,i)+t;
end
%now we will see how well the routine works by applying the mapping and
%measuring the square  distance between the desired and actual positions


%calculate mean squared distance from actual points
%sqDiff = mean(sum((pts2Cart-pts2EstCart).^2));

%draw figure with points before and after
%draw two set of two dimensional points
%opens figure
figure; set(gcf,'Color',[1 1 1]);
%draw lines between each pair of points
nPoint = size(pts1Cart,2);
for (cPoint = 1:nPoint)
    %plot a green line pairs of actual and estimated points
    plot([pts2Cart(1,cPoint) pts2EstCart(1,cPoint)],[pts2Cart(2,cPoint) pts2EstCart(2,cPoint)],'g-');
    %make sure we don't replace with next point
    hold on;
end;

%draws second set of points
plot(pts2Cart(1,:),pts2Cart(2,:),'r.','MarkerSize',20);
%remove axis
set(gca,'Box','Off');

%draws estimated positions of second set of points
plot(pts2EstCart(1,:),pts2EstCart(2,:),'m.','MarkerSize',20);





%==========================================================================
function [R,t] = affine_transformation(pts1Cart, pts2Cart)

m = size(pts1Cart,2);
pts1CartHom = [pts1Cart; ones(1,m)];

A = zeros(m*2,6);
for i=1:m
    A((2*i-1),4:6)=0;
    A((2*i),1:3)=0;
    A((2*i-1),1:3)= pts1CartHom(1:3,i);
    A((2*i),4:6)= pts1CartHom(1:3,i);
end    
X = zeros(m*2,1);
for i = 1:m
    X((2*i-1):(2*i))= pts2Cart(:,i)';
end
b = inv(A'*A)*A'*X;

R = zeros(2,2);
R(1,:) = b(1:2);
R(2,:) = b(4:5);
t = zeros(2,1);
t(1) = b(3);
t(2) = b(6);


