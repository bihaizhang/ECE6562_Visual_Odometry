%% Visual odometry example
%   - stereo camera
%   - ICP between frames

close all

%% read images
left = iread('C:\robotic\bridge-l\*.pgm', 'roi', [20 750; 20 440]);
right = iread('C:\robotic\bridge-r\*.pgm', 'roi', [20 750; 20 440]);
% % ###############################################################################
% % # Camera parameter file                                                       #
% % ###############################################################################
% % 
% % [INTERNAL]
% % F        =  985.939 # [pixel] focal length
% % SX       =  1.0     # [pixel] pixel size in X direction
% % SY       =  1.0     # [pixel] pixel size in Y direction
% % X0       =  390.255 # [pixel] X-coordinate of principle
% % Y0       =  242.329 # [pixel] Y-coordinate of principle
% % 
% % [EXTERNAL]
% % B        =  0.20    # [m] width of baseline of stereo camera rig
% % X        = -0.83    # [m] distance of rectified images (virtual camera)
% % Y        =  0.00    # [m] lateral position of rectified images (virtual camera)
% % Z        =  1.28    # [m] height of rectified images (virtual camera)
% % TILT     =  0.0062  # [rad] tilt angle
% % YAW      =  0.0064  # [rad] yaw angle
% % ROLL     =  0.0009  # [rad] roll angle
% % 
% % # Notes:
% % #  In a stereo camera system the internal parameters for both cameras are the
% % #  same.
% % #
% % #  The camera position (X, Y, Z) is given in car coordinates.
% % #  For the definition of the camera and car coordinate system and the rotation 
% % #  angles see the image carcameracoord.png.



%%

f        =  985.939; % [pixel] focal length
u0       =  390.255; % [pixel] X-coordinate of principle
v0       =  242.329; % [pixel] Y-coordinate of principle
b        =  0.20;    % [m] width of baseline of stereo camera rig

cam = CentralCamera('focal', f, 'centre', [u0 v0], 'pixel', [1 1 ]);



randinit

% matching
for i=1:size(left,3)
    % for every frame
    L = left(:,:,i);
    R = right(:,:,i);
    
    % find corner features
    fl = icorner(L, 'nfeat', 400, 'patch', 7,  'suppress', 0);
    fr = icorner(R, 'nfeat', 400, 'patch', 7,  'suppress', 0);
    
    % robustly match left and right corner features
    % - stereo match
    [mstereo,Cs] = fl.match(fr); %mstereo: coordinates of feature points in image plane, Cs: label of those points
    Fstereo = mstereo.ransac(@fmatrix, 1e-4, 'retry', 5); %Fstereo: the fundamental matrix obtained by RANSAC


    % keep the features and match objects for the inliers
    k = Cs(1,mstereo.inlierx);  % index of inlier features
    fl = fl(k);
    
    mstereo = mstereo.inlier;
    
    % triangulate 3D points
    p = mstereo.p;
    p1 = p(1:2,:); p2 = p(3:4,:); %p1:左图中特征点的uv，p2:右图中特征点的uv
    d = p1(1,:) - p2(1,:);
    X = b * (p1(1,:) - u0) ./ d;
    Y = b * (p1(2,:) - v0) ./ d;
    Z = f * b ./ d;
    P = [X; Y; Z];
    count=1;
    mark=1;
    Pf=[];
    p1f=[];
    marker=[];
    for tt=1:size(P,2)
        if norm(P(:,tt))<=20000
            Pf(:,count)=P(:,tt);
            p1f(:,count)=p1(:,tt);
            count=count+1;
        else  marker(mark)=tt;
            mark=mark+1;
        end
    end
    
    if i > 1
        % if we have a previous frame
        
        % display two sequential stereo pairs
%         idisp([L R; Lp Rp], 'nogui');
        
        % show the stereo matching in the current frame
%         mstereo.plot('y', 'offset', [0 numcols(L) 0 0])
        
        % robustly match all the inliers from this frame with the inliers from previous frame
        % - temporal match
        c1w=sum(Pfp,2)/size(Pfp,2);
        A=Pfp'-c1w';
        [V,D]=eig(A'*A);
        c2w=c1w+sqrt(D(1,1))*V(:,1);
        c3w=c1w+sqrt(D(2,2))*V(:,2);
        c4w=c1w+sqrt(D(3,3))*V(:,3);
        alpha=[];
        for jj=1:size(Pfp,2)
            alpha(:,jj)=inv([c1w c2w c3w c4w;1 1 1 1])\[Pfp(:,jj);1];
        end
        M=[];
        for kk=1:2:2*size(Pfp,2)
            M(kk,:)=zeros(1,12);
            for mm=1:4
                M(kk,3*mm-2)=alpha(mm,(kk+1)/2)*f;
                M(kk,3*mm)=alpha(mm,(kk+1)/2)*(u0-p1fp(1,(kk+1)/2));
            end
            M(kk+1,:)=zeros(1,12);
            for nn=1:4
               M(kk+1,3*nn-1)=alpha(nn,(kk+1)/2)*f;
               M(kk+1,3*nn)=alpha(nn,(kk+1)/2)*(v0-p1fp(2,(kk+1)/2));
            end
        end

        [V,S]=eig(M'*M);
        S(1,1)
        X1=V(:,1);
        Alph=alpha';
        Cw=[c1w';c2w';c3w';c4w'];
        Xw=Pfp';
        [Cc,Xc,sc]=compute_norm_sign_scaling_factor(X1,Cw,Alph,Xw);
        Xc=Xc';
        
        for au=1:size(markerp)
            Xc=[Xc(:,1:markerp(au)-1) [30000;0;0] Xc(:,markerp(au):end)];
        end

        [mtemporal,Ct] = fl.match(flp);
        [Ftemporal,eFt] = mtemporal.ransac(@fmatrix, 1e-4, 'retry', 5);
        
        % now create a bundle adjustment problem
        
        ba = BundleAdjust(cam);
        
        c1 = ba.add_camera( SE3(), 'fixed' );  % first camera at origin (prev frame)
        c2 = ba.add_camera( SE3() );  %  second camera initialized at origin (current frame)
        
        for ii=mtemporal.inlierx
            j = Ct(2,ii);   % index to current frame feature and landmark
            if norm(Xc(:,j))<=20000
                lm = ba.add_landmark(Xc(:,j));
                p = mtemporal(ii).p;
                ba.add_projection(c2, lm, p(1:2));  % current camera
                ba.add_projection(c1, lm, p(3:4));  % previous camera
            end
            
        end
        
        % solve bundle adjustment, fix number of iterations
        [baf,e] = ba.optimize('iterations', 15);
        T(:,1)=baf.getcamera(2).n;
        T(:,2)=baf.getcamera(2).o;
        T(:,3)=baf.getcamera(2).a;
        T(:,4)=baf.getcamera(2).t;
        if i==2
            TT(:,:,i)=[T;0 0 0 1];
        else
        TT(:,:,i)=TT(:,:,i-1)*[T;0 0 0 1];
        end
        pose(:,i)=TT(:,:,i)*[0;0;0;1];
        ebundle(i)=e;                      
    end
    
    % keep images and features for next cycle
    flp = fl;
    Lp = L;
    Rp = R;
    Pp=P;
    p1p=p1; 
    Pfp=Pf;
    p1fp=p1f;
    markerp=marker;
    drawnow
end

%% process and display the results
% Eliminating error in frame 227 and 228
HT=TT;
HT(:,:,228)=HT(:,:,227)*(inv(TT(:,:,226))*TT(:,:,227));
for qq=228:250
    HT(:,:,qq+1)=HT(:,:,qq)*(inv(TT(:,:,qq))*TT(:,:,qq+1))
end
for ss=2:251
    fpose(:,ss)=HT(:,:,ss)*[0;0;0;1];
end
% show the result
plot(fpose(1,:),fpose(3,:))
hold on 
axis equal

