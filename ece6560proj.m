clc
close all
clear

left = iread('bridge-l/*.png', 'roi', [20 750; 20 440]);
%left = iread('city-l/*.pgm', 'roi', [30 735; 130 410]);
f        =  985.939; % [pixel] focal length
u0       =  390.255; % [pixel] X-coordinate of principle
v0       =  242.329; % [pixel] Y-coordinate of principle
b        =  0.20;    % [m] width of baseline of stereo camera rig

K = [f 0  u0;  
         0  f v0;
        0  0 1];
cam = CentralCamera('focal', f, 'centre', [u0 v0], 'pixel', [1 1 ]);    
path = [];
Rpos = [1 0 0; 0 1 0; 0 0 1;];
tpos = [0 0 0];
checkList = [];
FDlist= [];
count1 = 0;
count2= 0;
eList = [];
for i=2:size(left,3)
    
    left1 = left(:,:,i-1);
    left2 = left(:,:,i);
    %left1 = iread('Bridge-Left/image0001_c0.pgm', 'roi', [20 750; 20 440]);
    %left2 = iread('Bridge-Left/image0002_c0.pgm', 'roi', [20 750; 20 440]);
    %c = icorner(left, 'nfeat', 200, 'patch', 7);
    %l1 = icorner(left1,'nfeat', 200, 'patch', 7, 'suppress', 0);
    %l2 = icorner(left2,'nfeat', 200, 'patch', 7, 'suppress', 0);
    l1 = icorner(left1,'nfeat', 600, 'patch', 7, 'suppress', 0);
    l2 = icorner(left2,'nfeat', 600, 'patch', 7, 'suppress', 0);
    [matched,Cs]  = l1.match(l2);
  
    cd('C:\Users\Yulong Gu\AppData\Roaming\MathWorks\MATLAB Add-Ons\Toolboxes\Machine Vision Toolbox for MATLAB\vision');
%    [Ftemporal,eFt] = matched.ransac(@fmatrix, 1e-4, 'retry', 5);
    k1 = Cs(1,matched.inlierx);
    k2  = Cs(2,matched.inlierx);
    %matched = matched.inlier;
    p = matched.p; %matched points
    p1 = p(1:2,:); %matched points for l1
    p2 = p(3:4,:); %matched points for l2
    L1 = left1(:,:,1);
    L2 = left2(:,:,1);
   % idisp([L1 L2], 'nogui');
  %  matched.plot('y', 'offset', [0 numcols(L1) 0 0]);
    fMatrixPoints = [];
    fMatrixPoints2 = [];
    
    %Calculate fundamental matrix
    Fund = matched.ransac(@fmatrix,1e-4, 'verbose');
    len = length(p1);
    r1 = randi([1 len]);
    X_p1 = [p1(:,r1);1];
    X_n1 = [p2(:,r1);1];
    r2 = randi([1 len]);
    X_p2 = [p1(:,r2);1];
    X_n2 = [p2(:,r2);1];
    r3 = randi([1 len]);
    X_p3 = [p1(:,r3);1];
    X_n3 = [p2(:,r3);1];
    e1 = abs(X_p1' * Fund * X_n1);
    e2 = abs(X_p2' * Fund * X_n2);
    e3 = abs(X_p3' * Fund * X_n3);
    e = e1 + e2 + e3;
    eList = [eList;e];
   % FDlist =[FDlist;FD];
    %check = xprime'*Fund*x; 
    %checkList = [checkList,check];
   
    %Calculate essential matrix
    E = K' * Fund * K;
    sol = cam.invE(E);
    sol1 = double(sol(1)); sol2 = double(sol(2));
    R1= sol1(1:3,1:3); R2 = sol2(1:3,1:3);
    if (det(R1)<0)
        R1 = -R1;
    end
    if (det(R2)<0)
        R2 = -R2;
    end
    t1 = sol1(1:3,4);t2 = sol2(1:3,4);
    P2{1} = [R1 t1;0 0 0 1];
    P2{2} = [R1 t2;0 0 0 1];
    P2{3} = [R2 t1;0 0 0 1];
    P2{4} = [R2 t2;0 0 0 1];
    len = length(p1);
    r = randi([1 len],1,40);
    fMatrixPoints = p1(:,r);
    fMatrixPoints2 = p2(:,r);
    fMatrixPoints = [fMatrixPoints;ones(1,40)];
    fMatrixPoints2 = [fMatrixPoints2;ones(1,40)];
   
    P1 = eye(4);
    N = [];
for k = [1 2 3 4]
    % Triangulating by using all four camera matrices P2
    cd('D:\ECE 6562 Robotics\4560project(demo vedio)');
    X{k} = triangulationCustom(fMatrixPoints,fMatrixPoints2, K,P1,P2{k} );    
    % Projecting the triangulated scene point X on camera P1
    xp{k} = P1*X{k}; 
    % Projecting the triangulated scene point X on camera P2
    x2p{k} = P2{k}*X{k};  
    % Counting all the points in front of camera pairs P1 and P2
    N = [N sum(x2p{k}(3,:)>0)+sum(xp{k}(3,:)>0)]
end
    [value,index]=max(N);  Pgood = P2{index};
    
    rotGood = Pgood(1:3,1:3)
   
    tGood = Pgood(1:3,4)
    if (value < 43)
        %Pgood
       % Pgood(1:3,1:3) = eye(3);
        %Pgood(1:3,4) = 0;
        count1 = count1 + 1
    end
   if(abs(Pgood(3,4)) < 0.001)
        %Pgood(1:3,1:3) = eye(3);
        %Pgood(1:3,4) = 0;
        count2 = count2 + 1;
    end
    
    if (tGood(3) < 0)
        tGood = -tGood;
    end
    if (norm(tGood)~=0)
        tGood = tGood/norm(tGood);
    end
    
    Rpos = rotGood *Rpos;
    tpos = tpos + tGood'*Rpos;
    
    path = [path; tpos];

    
end
plot(path(:,1),path(:,3))
axis equal