function [x,y] = ThreeDcoordToTwoD1(QP,Emesh)

% 1. Find electrode triangle which this point belongs to
Etriangles = Emesh.ConnectivityList;


NEtriangles = max(size(Etriangles));
stop = 0;
j= 1;
err = 0;
count = 0;
TRI = [];
while (j <=NEtriangles) && stop==0
	b = cartesianToBarycentric(Emesh, j, QP);
    if all(b >= -0.00001)
        QPbary = b;
        TRI = [TRI j];
        %stop = 1;
        count = count+1;
    end
    j = j+1;
    if j > NEtriangles && count == 0
        disp("ERROR: PS not on shell");
        err = 1;
    end
end

if count > 1
%     %Variables to store distance to each triangle and coordinates of the
%     %closest point in each triangle
     dist_cum=[];
     for j=1:1:length(TRI) %For all the triangles
         [dist] = pointTriangleDistance([Emesh.Points(Etriangles(TRI(j),1),:);Emesh.Points(Etriangles(TRI(j),2),:); Emesh.Points(Etriangles(TRI(j),3),:)],QP); %Measure distance and closest point
         dist_cum=[dist_cum; dist]; %Store all distances         
     end
    [~, Idef]=min(dist_cum); %Find index with shortest distance    
    TRI = TRI(Idef);
end

if err == 1
    x = NaN;y=NaN; % won't show on 2D map
else
    % 2. Figure out spline and electrode number for each edge of the triangle
    TR = Etriangles(TRI,:); % electrode indices 
    splinesI = floor((TR-1)/8)+1;
    electrodesI = TR - ((splinesI-1)*8);

    % 3a. Check if PS does not lie on areas non represented in 2D
    hasSpline1 = ~(sum(ismember(splinesI,1))==0);
    hasSpline8 = ~(sum(ismember(splinesI,8))==0);
    hasBothSpline1and8  = hasSpline1 && hasSpline8;
    

    if sum(electrodesI) == 3 || sum(electrodesI) == 24 || hasBothSpline1and8
        hola=1;
        x = NaN; y = NaN; % won't show on 2D map
    else
    % 3b. Figure out coordinates on 64x64 map

        x = 1+ (splinesI-1)*9;
        y = 1+ (electrodesI-1)*9;

        x = sum(QPbary.*x);
        y = sum(QPbary.*y);
    end
end




