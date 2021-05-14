%% Visual odometry example
%   - stereo camera
%   - ICP between frames

close all

%% read images
%  left = iread('C:\robotic\city-l\*.pgm', 'roi', [30 735; 130 410]);
%  right = iread('C:\robotic\city-r\*.pgm', 'roi', [30 735; 130 410]);
left = iread('C:\robotic\bridge-l\*.pgm', 'roi', [20 750; 20 440]);
right = iread('C:\robotic\bridge-r\*.pgm', 'roi', [20 750; 20 440]);
% left=iread('C:\robotic\overhead-l\*.pgm','roi',[20 750; 20 440]);
% right = iread('C:\robotic\overhead-r\*.pgm', 'roi', [20 750; 20 440]);
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


clear T ebundle efund
randinit

% matching
for i=1:size(left,3)
    % for every frame
    L = left(:,:,i);
    R = right(:,:,i);
    
    % find corner features
    fl = icorner(L, 'nfeat', 400, 'patch', 7,  'suppress', 0);
    fr = icorner(R, 'nfeat', 400, 'patch', 7,  'suppress', 0);
%     counter=size(fl,2);
%     jj=1;
%     for de=1:counter
%         if fl(de).u<=350&&(fl(de).v>=100&&fl(de).v<=150)
%             delete(jj)=de;
%             jj=jj+1;
%         end
%     end
%     fl(:,delete)=[];

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
    
     if i > 1
        % if we have a previous frame
        
        % display two sequential stereo pairs
%         idisp([L R; Lp Rp], 'nogui');
        
        % show the stereo matching in the current frame
%         mstereo.plot('y', 'offset', [0 numcols(L) 0 0])
        
        % robustly match all the inliers from this frame with the inliers from previous frame
        % - temporal match
        [mtemporal,Ct] = fl.match(flp);
        [Ftemporal,eFt] = mtemporal.ransac(@fmatrix, 1e-4, 'retry', 5);
        
        % and plot them
%         mtemporal.inlier.plot('y', 'offset', [0 0 0 numrows(L)])
         
%         if i == 10
%             rvcprint(svg')
%         end
        
        % now create a bundle adjustment problem
        
        ba = BundleAdjust(cam);
        
        c1 = ba.add_camera( SE3(), 'fixed' );  % first camera at origin (prev frame)
        c2 = ba.add_camera( SE3() );  %  second camera initialized at origin (current frame)
        
        for ii=mtemporal.inlierx
            j = Ct(1,ii);   % index to current frame feature and landmark
            lm = ba.add_landmark(P(:,j));
            p = mtemporal(ii).p;
            ba.add_projection(c2, lm, p(1:2));  % current camera
            ba.add_projection(c1, lm, p(3:4));  % previous camera
        end
        
        % solve bundle adjustment, fix number of iterations
        [baf,e] = ba.optimize('iterations', 10);
%         Tan=baf.getcamera(2);
%         angle(i)=180/pi*atan(-Tan.n(3)/sqrt(Tan.n(1)^2+Tan.n(2)^2)); %角度制
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
                       
        tz(i) = baf.getcamera(2).t(3); % translation in z axis between the previous frame and current frame
        %利用getcamera去找相邻两帧之间的位置变化从而得到每一帧对应的相机坐标
        ebundle(i) = e;
        efund(i) = eFt;
    end
    
    % keep images and features for next cycle
    flp = fl;
    Lp = L;
    Rp = R;
    
    drawnow
end

%% process the results
 tz2 = tz;
 tztest=tz;
 tz2(ebundle>20) = NaN;
 pose2=pose;
 posetest=pose;
 pose2(ebundle>20)=NaN;

% %中位值平均滤波
% N=10;
% for t=N+1:size(tz,2)
%         pose(:,t)=sum(pose(:,t-N+1:t),2)/N;
% %     tz(t)=(sum(tztest(t-N+1:t))-min(tztest(t-N+1:t))-max(tztest(t-N+1:t)))/(N-2);
% %     tx(t)=(sum(tx(t-N+1:t))-min(tx(t-N+1:t))-max(tx(t-N+1:t)))/(N-2);
% end
% for p=1:size(tz,2)
%     odometryx(p)=sum(tx(1:p));
%     odometryz(p)=sum(tz(1:p));
% end
% for t=N+1:size(tz2,2)
%     tz2(t)=(sum(tz2(t-N+1:t))-min(tz2(t-N+1:t))-max(tz2(t-N+1:t)))/(N-2);
%     tx2(t)=(sum(tx2(t-N+1:t))-min(tx2(t-N+1:t))-max(tx2(t-N+1:t)))/(N-2);
% end

% f=2;
% fpose(:,1)=[0;0;0;1];
% for nb=2:size(tz,2)
%     if ebundle(nb)<=20&&tz(nb)~=0
%         fpose(:,f)=pose(:,nb);
%         f=f+1;
%     end
% end
plot(pose(1,:),pose(3,:))
hold on
axis equal
% figure
% subplot(311); plot(tz2, '.-', 'MarkerSize', 15)
% yaxis(0,1.5)
% xaxis(2,length(tz))
% ylabel('camera displacement (m)');
% grid
% subplot(312); plot(ebundle, '.-', 'MarkerSize', 15)
% set(gca, 'YScale', 'log')
% xaxis(2,length(tz))
% hold on
% plot([2,length(tz)], [20 20], 'r--');
% xlabel('Time step')
% ylabel('total error (pix^2)');
% grid
% subplot(313);plot(angle,'.-', 'MarkerSize', 15)
% yaxis(0,20)
% xaxis(2,length(tz))
% grid
% plot(pose(1,:),pose(3,:))
% hold on
% axis equal
% rvcprint()

% median(tz(ebundle<20))