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
eList  = [];
count1 = 0;
count2 = 0;
count3 = 0;
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
    
    p = matched.p; %matched points
    p1 = p(1:2,:); %matched points for l1
    p2 = p(3:4,:); %matched points for l2
    L1 = left1(:,:,1);
    L2 = left2(:,:,1);
    fMatrixPoints = [];
    fMatrixPoints2 = [];
    len = length(p1);
    error = 5;
for j=1:1200  
    r = randi([1 len],1,8);
    fMatrixPoints = p1(:,r);
    fMatrixPoints2 = p2(:,r);
      %Points from frame1
    x1 = fMatrixPoints(1,1); y1 = fMatrixPoints(2,1);
    x2 = fMatrixPoints(1,2); y2 = fMatrixPoints(2,2);
    x3 = fMatrixPoints(1,3); y3 = fMatrixPoints(2,3);
    x4 = fMatrixPoints(1,4); y4 = fMatrixPoints(2,4);
    x5 = fMatrixPoints(1,5); y5 = fMatrixPoints(2,5);
    x6 = fMatrixPoints(1,6); y6 = fMatrixPoints(2,6);
    x7 = fMatrixPoints(1,7); y7 = fMatrixPoints(2,7);
    x8 = fMatrixPoints(1,8); y8 = fMatrixPoints(2,8);
        %Points from frame2
    x1P = fMatrixPoints2(1,1); y1P = fMatrixPoints2(2,1);
    x2P = fMatrixPoints2(1,2); y2P = fMatrixPoints2(2,2);
    x3P = fMatrixPoints2(1,3); y3P = fMatrixPoints2(2,3);
    x4P = fMatrixPoints2(1,4); y4P = fMatrixPoints2(2,4);
    x5P = fMatrixPoints2(1,5); y5P = fMatrixPoints2(2,5);
    x6P = fMatrixPoints2(1,6); y6P = fMatrixPoints2(2,6);
    x7P = fMatrixPoints2(1,7); y7P = fMatrixPoints2(2,7);
    x8P = fMatrixPoints2(1,8); y8P = fMatrixPoints2(2,8);
    
  A  =[x1*x1P x1*y1P x1 y1*x1P y1*y1P y1 x1P y1P 1; ...
    x2*x2P x2*y2P x2 y2*x2P y2*y2P y2 x2P y2P 1;
    x3*x3P x3*y3P x3 y3*x3P y3*y3P y3 x3P y3P 1;
    x4*x4P x4*y4P x4 y4*x4P y4*y4P y4 x4P y4P 1;
    x5*x5P x5*y5P x5 y5*x5P y5*y5P y5 x5P y5P 1;
    x6*x6P x6*y6P x6 y6*x6P y6*y6P y6 x6P y6P 1;
    x7*x7P x7*y7P x7 y7*x7P y7*y7P y7 x7P y7P 1;
    x8*x8P x8*y8P x8 y8*x8P y8*y8P y8 x8P y8P 1;];

    A = double(A);
    [U, D, V] = svd(A);

    Ftemp = reshape(V(:, end), 3, 3)';
    [Utemp,Dtemp,Vtemp] = svd(Ftemp);
    Dtemp(3,3)=0;
    Ftemp=Utemp*Dtemp*Vtemp';
    r1 = randi([1 len]);
    X_p1 = [p1(:,r1);1];
    X_n1 = [p2(:,r1);1];
    r2 = randi([1 len]);
    X_p2 = [p1(:,r2);1];
    X_n2 = [p2(:,r2);1];
    r3 = randi([1 len]);
    X_p3 = [p1(:,r3);1];
    X_n3 = [p2(:,r3);1];
    e1 = abs(X_p1' * Ftemp * X_n1);
    e2 = abs(X_p2' * Ftemp * X_n2);
    e3 = abs(X_p3' * Ftemp * X_n3);
    e = e1 + e2 + e3;
    if (e < error)
        error = e;
        F = Ftemp;
    end   
end
    eList = [eList;error];

    E = K' * F * K;
    [Ue,De,Ve]=svd(E);
    De(1,1)=1;De(2,2)=1;De(3,3)=0;
    E = Ue * De* Ve';
    
    [Ue2,De2,Ve2]=svd(E);
   
    P1 = eye(4);
    
    t1 = vex(Ue2*rotz(pi/2)*De2*Ue2');
    t2 = vex(Ue2*rotz(-pi/2)*De2*Ue2');

    W = rotz(pi/2);
    if(t1(3)>0)
        t_temp = t1;
    end
    if(t2(3)>0)
        t_temp = t2;
    end
    R_temp1 = Ue2*W*Ve2';
    R_temp2 = Ue2*W'*Ve2';
    if (det(R_temp1)<0)
        R_temp1 = -R_temp1;
    end
    if (det(R_temp2)<0)
        R_temp2 = -R_temp2;
    end
    P2{1} = [R_temp1 t1];
    P2{2} = [R_temp2 t1];
    P2{3} = [R_temp1 t2];
    P2{4} = [R_temp2 t2];
   % p2{1} = inv(p2{1});
    %p2{2} = inv(p2{2});
    N = [];
    testNum= randi([1 len],1,40);
    testPoints = p1(:,testNum);
    testPoints2 = p2(:,testNum);
    testPoints = [testPoints;ones(1,40)];
    testPoints2 = [testPoints2;ones(1,40)];
for k = [1 2 3 4]
    
    P2{k} = [P2{k};0,0,0,1];
    %P2{k} = inv(P2{k});
    X{k} = triangulationCustom(testPoints,testPoints2, K,P1,P2{k} );
    
    % Projecting the triangulated scene point X on camera P1
    xp{k} = P1*X{k};
    %P2{k} = inv(P2{k});
    % Projecting the triangulated scene point X on camera P2
    x2p{k} = P2{k}*X{k};
    
    % Counting all the points in front of camera pairs P1 and P2
    N = [N sum(x2p{k}(3,:)>0)+sum(xp{k}(3,:)>0)]
end
    [value,index]=max(N);
        

    %X = X{index};
    Pgood = P2{index};
   % recheck = P2{index}*X{index};
   % recheck2 = P1*X{index};
   % if (sum(recheck(3,:)>0)<22 || sum(recheck2(3,:)>0)<22 )
   if (value < 43)
        %Pgood
       % Pgood(1:3,1:3) = eye(3);Pgood(1:3,4) = 0;
        count1 = count1 + 1
    end
 
    if(abs(Pgood(3,4)) < 0.001)
        Pgood(1:3,1:3) = eye(3);
        Pgood(1:3,4) = 0;
        count2 = count2 + 1;
    end
     if (Pgood(1,1) < 0 || Pgood(2,2) < 0 || Pgood(3,3) < 0 )
       % Pgood(1:3,1:3) = eye(3);
        count3 = count3+1;
    end
    rotGood = Pgood(1:3,1:3);
    tGood = Pgood(1:3,4);
    if (tGood(3) < 0)
        tGood = -tGood;
    end
    %if(norm(tGood)!=0)
    if (norm(tGood)~=0)
        tGood = tGood/norm(tGood);
    end
    Rpos = rotGood *Rpos;
    tpos = tpos + tGood'*Rpos;
    
    path = [path; tpos];
end
plot(path(:,1),path(:,3))
axis equal
