function X = triangulationCustom(matched1,matched2,K, P1, P2)
    x11 = inv(K)*matched1;
    x22 = inv(K)*matched2;
    %P2 = inv(P2);
for i=1:size(x11,2)
    sx1 = x11(:,i);
    sx2 = x22(:,i);
    
    
    A1 = sx1(1,1).*P1(3,:) - P1(1,:);
    A2 = sx1(2,1).*P1(3,:) - P1(2,:);
    A3 = sx2(1,1).*P2(3,:) - P2(1,:);
    A4 = sx2(2,1).*P2(3,:) - P2(2,:);
    
    %Set A matric, and find SDV(A)
    A = [A1;A2;A3;A4];
    
     A1n = sqrt(sum(A(1,:).*A(1,:)));
     A2n = sqrt(sum(A(2,:).*A(2,:)));
     A3n = sqrt(sum(A(1,:).*A(1,:)));
     A4n = sqrt(sum(A(1,:).*A(1,:))); 
     Anorm = [A(1,:)/A1n;
                 A(2,:)/A2n;
                 A(3,:)/A3n;
                 A(4,:)/A4n];
    [~,~,V] = svd(Anorm);
    %Point in 3D space is the last column of V
    X_temp = V(:,4);
    X_temp = X_temp./X_temp(4);
    
    X(:,i) = X_temp;
    
end