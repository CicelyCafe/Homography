function practical2b

%The goal of this part of the practical is to take a real image containing
%a planar black square and figure out the transformation between the square
%and the camera.  We will then draw a wire-frame cube with it's base
%corners at the corner of the square.  You should use this
%template for your code and fill in the missing sections marked "TO DO"

%load in image 
im = imread('test104.jpg');

%define points on image
xImCart = [  140.3464  212.1129  346.3065  298.1344   247.9962;...
             308.9825  236.7646  255.4416  340.7335   281.5895];
         
%define 3D points of plane
XCart = [-50 -50  50  50 0 ;...
          50 -50 -50  50 0;...
           0   0   0   0 0];

%We assume that the intrinsic camera matrix K is known and has values
K = [640  0    320;...
     0    640  240;
     0    0    1];

%draw image and 2d points
figure; set(gcf,'Color',[1 1 1]);
imshow(im); axis off; axis image; hold on;
plot(xImCart(1,:),xImCart(2,:),'r.','MarkerSize',10);
       
%TO DO Use your routine to calculate TEst, the extrinsic matrix relating the
%plane position to the camera position.

T = estimatePlanePose(xImCart,XCart,K)


%define 3D points of plane
XWireFrameCart = [-50 -50  50  50 -50 -50  50  50;...
                   50 -50 -50  50  50 -50 -50  50;...
                    0   0   0   0 -100 -100 -100 -100];
                
WireFrameCart = projectiveCamera(K,T,XWireFrameCart)
plot(WireFrameCart(1,:), WireFrameCart(2,:),'y.','Markersize',10);

for (cPoint = 1:4)
    plot([WireFrameCart(1,cPoint) WireFrameCart(1,cPoint+4)],[WireFrameCart(2,cPoint) WireFrameCart(2,cPoint+4)],'y-');
    %make sure we don't replace with next point
    hold on;
    %end
end;
for cPoint = [1:3 5:7]
    plot([WireFrameCart(1,cPoint) WireFrameCart(1,cPoint+1)],[WireFrameCart(2,cPoint) WireFrameCart(2,cPoint+1)],'y-');
    hold on
end
plot([WireFrameCart(1,5) WireFrameCart(1,8)],[WireFrameCart(2,5) WireFrameCart(2,8)],'y-');
plot([WireFrameCart(1,1) WireFrameCart(1,4)],[WireFrameCart(2,1) WireFrameCart(2,4)],'y-');
%TO DO Draw a wire frame cube, by projecting the vertices of a 3D cube
%through the projective camera and drawing lines betweeen the resulting 2d image
%points



function xImCart = projectiveCamera(K,T,XCart);

%replace this
%xImCart = [];

%TO DO convert Cartesian 3d points XCart to homogeneous coordinates XHom
XHom = [XCart; ones(1,size(XCart,2))];
%TO DO apply extrinsic matrix to XHom to move to frame of reference of
%camera
XexHom = T*XHom
%TO DO project points into normalized camera coordinates xCamHom by (achieved by
%removing fourth row)
xCamHom = XexHom(1:3,:)
%TO DO move points to image coordinates xImHom by applying intrinsic matrix
xImHom = K*xCamHom;
%TO DO convert points back to Cartesian coordinates xImCart
xImCart = xImHom(1:2,:)./repmat(xImHom(3,:),2,1);


%==========================================================================
%==========================================================================

%goal of function is to estimate pose of plane relative to camera
%(extrinsic matrix) given points in image xImCart, points in world XCart
%and intrinsic matrix K.

function T = estimatePlanePose(xImCart,XCart,K)

%replace this
%T = []

%TO DO Convert Cartesian image points xImCart to homogeneous representation
%xImHom
xImHom = [xImCart; ones(1,size(xImCart,2))];
%TO DO Convert image co-ordinates xImHom to normalized camera coordinates
%xCamHom
xCamHom = inv(K)*xImHom;
%TO DO Estimate homography H mapping homogeneous (x,y)
%coordinates of positions in real world to xCamHom.  Use the routine you wrote for
%Practical 1B.
xCamCart = xCamHom(1:2,:)./repmat(xCamHom(3,:),2,1);
H = calcBestHomography(XCart(1:2,:),xCamCart);
%TO DO Estimate first two columns of rotation matrix R from the first two
%columns of H using the SVD
[U,S,V] = svd(H(:,1:2));
R = zeros(3,3);
I = [1 0;0 1;0 0];
R(:,1:2) = U*I*V';
%TO DO Estimate the third column of the rotation matrix by taking the cross
%product of the first two columns
R(:,3) = cross(R(:,1),R(:,2))

%TO DO Check that the determinant of the rotation matrix is positive - if
%not then multiply last column by -1.
if det(R)<0
    R(:,end) = -R(:,end);
end;
%TO DO Estimate the translation t by finding the appropriate scaling factor k
%and applying it to the third colulmn of H
k = sum(sum(H(:,1:2)./R(:,1:2)))/6;
t = H(:,3)/k;
%TO DO Check whether t_z is negative - if it is then multiply t by -1 and
%the first two columns of R by -1.
if t(3)<0
    t = -t;
    R(:,1:2)= -R(:,1:2);
end;
%assemble transformation into matrix form
T  = [R t;0 0 0 1];

%==========================================================================
function H = calcBestHomography(pts1Cart, pts2Cart)

%should apply direct linear transform (DLT) algorithm to calculate best
%homography that maps the points in pts1Cart to their corresonding matchin in 
%pts2Cart

%****TO DO ****: replace this
%H = eye(3);
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
%**** TO DO ****;
%==========================================================================
function x = solveAXEqualsZero(A);
[U,S,V] = svd(A);
x = V(:,end);
%****TO DO **** Write this routine 

%QUESTIONS TO THINK ABOUT...

%Do the results look realistic?
%If not, then what factors do you think might be causing this?

