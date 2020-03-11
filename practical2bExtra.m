function practical2bExtra

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
          50 -50 -50  50 0]

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
XWireFrameCart = [-50 -50  50  50;...
                   50 -50 -50  50]

                
WireFrameCart = projectiveCamera(K,T,XWireFrameCart)
plot(WireFrameCart(1,:), WireFrameCart(2,:),'y.','Markersize',10);

for (cPoint = 1:3)
    plot([WireFrameCart(1,cPoint) WireFrameCart(1,cPoint+1)],[WireFrameCart(2,cPoint) WireFrameCart(2,cPoint+1)],'y-');
    %make sure we don't replace with next point
    hold on;
    %end
end;
plot([WireFrameCart(1,1) WireFrameCart(1,4)],[WireFrameCart(2,1) WireFrameCart(2,4)],'y-');
hold off 
%TO DO Draw a wire frame cube, by projecting the vertices of a 3D cube
%through the projective camera and drawing lines betweeen the resulting 2d image
%points



function xImCart = projectiveCamera(K,T,XCart);

%replace this
%xImCart = [];

%TO DO convert Cartesian 3d points XCart to homogeneous coordinates XHom
XHom = [XCart; ones(1,size(XCart,2))]
%TO DO apply extrinsic matrix to XHom to move to frame of reference of
%camera
XexHom = T*XHom
%TO DO project points into normalized camera coordinates xCamHom by (achieved by
%removing fourth row)
xCamHom = XexHom(1:3,:);
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

xImHom = [xImCart; ones(1,size(xImCart,2))];
%TO DO Convert image co-ordinates xImHom to normalized camera coordinates
%xCamHom
xCamHom = inv(K)*xImHom;
%TO DO Estimate homography H mapping homogeneous (x,y)
%coordinates of positions in real world to xCamHom.  Use the routine you wrote for
%Practical 1B.
xCamCart = xCamHom(1:2,:)./repmat(xCamHom(3,:),2,1);

T = affine_transformation(XCart(1:2,:),xCamCart)



function T = affine_transformation(pts1Cart, pts2Cart)

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

T = [R t;0 0 1];



