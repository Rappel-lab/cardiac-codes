function [coefficient, ps_registry] = interpolation_method(patientnum, side, ...
                                            NFRAMES_VID, numsub, ...
                                            datapath, patientpath, ...
                                            figpath)


%Interpolation Method - Function to plot and map phase singularities
% Syntax:  coefficient = interpolation_method(patientnum, side, ...
%                                            NFRAMES_VID, numsub, ...
%                                            datapath, patientpath, ...
%                                            figpath)
%
% Inputs:
%    patientnum - Patient Identifier
%    side - 0 for right atrium or 1 for left atrium
%    NFRAMES_VID - Integer indicating number of frames used to generate 
%    phase density maps. If 0, all frames are taken. Defaults to 4000
%    numsub - Integer indicating the number of subdivisions to perform in
%    the electrode mesh. Defaults to 3.
%    datapath - String indicating the directory containing the anatomic
%    files and phase values. Defaults to 'DATA\'
%    patientpath - String indicating the directory containing raw data for
%    the patient. Defaults to 'Patient_Data\'
%    figpath - String indicating the directory where output figures should
%    be saved. Defaults to 'Results\'
%
% Outputs:
%
%   coefficient - coefficient of the correlation between the number of 
%   phase singularities in 2D and in 3D.
%
% Author: Ricardo Abad, Orvil Collart
% Work address
% email: 
% October 2020; Last revision: 12-Oct-2020


if ~exist('datapath', 'var') || isempty(datapath)
    datapath = 'DATA/';
end
if ~exist('patientpath', 'var') || isempty(patientpath)
    patientpath = 'Patient_Data/';
end
if ~exist('figpath', 'var') || isempty(figpath)
    figpath = 'Results/';
end
if ~exist('NFRAMES_VID', 'var') || isempty(NFRAMES_VID)
    NFRAMES_VID = 4000;
end
if ~exist('numsub', 'var') || isempty(numsub)
    numsub = 3;
end


%% INITIALIZATION


%Load termination and control sites
load(strcat(datapath,'termsite_', patientnum, '.mat'));
load(strcat(datapath,'termsite2D_', patientnum, '.mat'));
TS2D=TS2D*2.25-1.25;

%If control site does not exist, assign termination site
try
 load(strcat(datapath,'controlsite_', patientnum, '.mat'));
 catch
     CS=TS;
end

%Electrode locations 
load(strcat(datapath,'shadows_', patientnum, '.mat'));

%% READ ANATOMY

disp('Loading anatomy file')
if not(exist('anatomy','var'));
    filepath = strcat(patientpath, patientnum, ...
                      '/NavX_Export/ModelGroups.xml');
    anatomy = xmlCoordExtract(filepath,side+1); 
end

%% Generate & Downsample Atrial Mesh
if side==1 %Left atrium
vertices = anatomy.Object2;
faces = anatomy.Polygons2;
end

if side==0 %Right atrium
vertices = anatomy.Object1;
faces = anatomy.Polygons1;        
end       
        
[faces,vertices] = generateMesh(faces,vertices,1000,1000,20); 
TR = triangulation(faces,vertices);

%% Manual triangulation of electrodes basket

% Define manually the triangulation of the electrodes basket. This will be 
% the triangulation of the coarse mesh 

% Trianguation in 3D
faces1=[];
for i=1:7
    for j=0:7
    faces1=[faces1; i+j*8, i+j*8+8, i+j*8+9; i+j*8, i+j*8+1, i+j*8+9];
    end
end

%%%Correction to make the basket closed
j=1;
for i=65:72
    faces1(faces1==i)=j;
    j=j+1;
end

%%%Triangulation in the edges at top and bottom
edge1=[1,9,17; 17,25,33; 33,41, 49; 49, 57, 1; 1, 17, 49; 17 , 49, 33];
edge2=edge1+7;
faces1=[faces1; edge1; edge2];

%%%Triangulation in 2D
faces2D=[];
for i=1:7
    for j=0:6
    faces2D=[faces2D; i+j*8, i+j*8+8, i+j*8+9; i+j*8, i+j*8+1, i+j*8+9];
    end
end

%%%Triangulation in 2D when it is interpolated (64x64)
faces2D_ip=[];
for i=1:63
    for j=0:62
    faces2D_ip=[faces2D_ip; i+j*64, i+j*64+64, i+j*64+65; i+j*64, i+j*64+1, i+j*64+65];
    end
end

% faces1 contains the triangulation of the electrodes
 Nfaces = max(size(faces1));  
 FaceVertexAlphaData = ones(Nfaces,1);

%% Subdivide EMesh
% Each subdivision multiplies 4 fold
vx=electrodes;
fx=faces1;
if numsub>0   
    for i=1:numsub
      [vx fx] = subdivide_tri(vx, fx);
    end
end

[X,Y] = meshgrid(1:64,1:64);
plane_points=[X(:), Y(:)]; %%% This variable includes all points

[X,Y] = meshgrid(linspace(1,64,8),linspace(1,64,8));
plane_tri_points=[X(:), Y(:)]; %%%This variables includes only the vertices of the big triangles

%% Triangulation of meshes
Emesh = triangulation(faces1, electrodes); %Electrodes mesh
Emesh_subdiv = triangulation(fx, vx); %Electrodes mesh subdivided num_sub times
plane_subdiv=triangulation(faces2D_ip, [plane_points, zeros(64*64,1)]); % Two-dimensional grid subdivided (64x64)
plane_triang=triangulation(faces2D, [plane_tri_points, zeros(64,1)]); %Two-dimensional grid not subdivided (8x8)

%% Matching of points between meshes

% Match each point on the subdivided electrode mesh with a big triangle
ti=[];              % This variable will contain, for each point in the subdivided electrode mesh,
                    % the index of the triangle that contains it 
coordinates=[];     % This variable will contain the points in the subdivided electorde mesh

disp('Assigning subdivided vertices to original triangles in 3D')

for i=1:size(vx,1)
[row, col]=find(fx==i, 1);
ti(i)=mod(row, 124);
end

ti=ti';
ti(ti==0)=124;
coordinates=vx;

disp('Assigning subdivided vertices to original triangles in 2D')
%%% Same as previous loop for 2D
ti_plane=[];
for i=1:size(plane_points,1) %For all the vertices in the refined mesh
    %Variables to store distance to each triangle and coordinates of the
    %closest point in each triangle
    dist_cum=[];
    PPO_cum=[];
    for j=1:size(faces2D,1) %For all the triangles
        [dist,PPO] = pointTriangleDistance([[plane_triang.Points(faces2D(j,1),:)];[plane_triang.Points(faces2D(j,2),:)]; plane_triang.Points(faces2D(j,3),:)],[plane_points(i,:), 0]); %Measure distance and closes point
        dist_cum=[dist_cum; dist]; %Store all distances
        PPO_cum=[PPO_cum; PPO];  %Store all the closest points in each triangle
    end
    [~, sidx]=min(dist_cum); %Find index with shortest distance
    ti_plane=[ti_plane; sidx]; %Store index of corresponding triangle
    coordinates=vx; %Store 
end

%% Calculation of barycentric coordinates for each point in the coarse mesh
%Calculation of the area of each triangle

disp('Calculate area of each triangle to obtain barycentric coordinates')
areas_ti=[];
for i = 1:size(faces1,1)
    X=[Emesh.Points(faces1(i,1), 1), Emesh.Points(faces1(i,2), 1), ...
       Emesh.Points(faces1(i,3), 1)];
    Y=[Emesh.Points(faces1(i,1), 2), Emesh.Points(faces1(i,2), 2), ... 
        Emesh.Points(faces1(i,3), 2)];
    Z=[Emesh.Points(faces1(i,1), 3), Emesh.Points(faces1(i,2), 3), ...
        Emesh.Points(faces1(i,3), 3)];
    plarea = area3D(X,Y,Z);
    areas_ti=[areas_ti; plarea];
end

%Calculation of the barycentric coordinates
disp('Obtaining barycentric coordinates for all subdivided vertices')
barycentric=[];
for i = 1:length(ti)
    %Coordinate xlam1
    X1=[coordinates(i,1), Emesh.Points(faces1(ti(i),2), 1), Emesh.Points(faces1(ti(i),3), 1)];
    Y1=[coordinates(i,2), Emesh.Points(faces1(ti(i),2), 2), Emesh.Points(faces1(ti(i),3), 2)];
    Z1=[coordinates(i,3), Emesh.Points(faces1(ti(i),2), 3), Emesh.Points(faces1(ti(i),3), 3)];
    plarea1 = area3D(X1,Y1,Z1);
    xlam1=plarea1/areas_ti(ti(i));
    
    %Coordinate xlam2
    X2=[Emesh.Points(faces1(ti(i),1), 1), coordinates(i,1), Emesh.Points(faces1(ti(i),3), 1)];
    Y2=[Emesh.Points(faces1(ti(i),1), 2), coordinates(i,2), Emesh.Points(faces1(ti(i),3), 2)];
    Z2=[Emesh.Points(faces1(ti(i),1), 3), coordinates(i,3), Emesh.Points(faces1(ti(i),3), 3)];
    plarea2 = area3D(X2,Y2,Z2);
    xlam2=plarea2/areas_ti(ti(i));
    
    %Coordinate xlam3
    X3=[Emesh.Points(faces1(ti(i),1), 1), Emesh.Points(faces1(ti(i),2), 1), coordinates(i,1)];
    Y3=[Emesh.Points(faces1(ti(i),1), 2), Emesh.Points(faces1(ti(i),2), 2), coordinates(i,2)];
    Z3=[Emesh.Points(faces1(ti(i),1), 3), Emesh.Points(faces1(ti(i),2), 3), coordinates(i,3)];
    plarea3 = area3D(X3,Y3,Z3);
    xlam3=plarea3/areas_ti(ti(i));
    
    barycentric=[barycentric; xlam1, xlam2, xlam3];
    
end    

%Replace NaN with 0
barycentric(isnan(barycentric))=0;

%% Extract phase/activation data

disp('Reading phase data')
phase=load(strcat(datapath,'phase_', patientnum, '.mat')); %%% This is the signal after Kuklik transform
phase=double(phase.rsig);
maxval=2*pi;


Nvertices=max(size(vx)); %Number of vertices in the subdivided mesh

NFRAMES_VID=min([NFRAMES_VID, size(phase, 2)]); %Number of frames in case 
% the introduced number is bigger than the actual number of frames

if NFRAMES_VID == 0 
    NFRAMES_VID = max(size(phase));  
end 
 
C = zeros(1,Nvertices,NFRAMES_VID); 
%% Calculation of phase for the grid (2D)

disp('Calculating 2D interpolation')
%%% Interpolation of the signal after Kuklik into the 2D mesh
rsig=phase;
rsig_int = nan(64,64,size(rsig,2)); % This will contain the interpolated signal
[X,Y] = meshgrid(linspace(1,64,8),linspace(1,64,8)); 
[Xq,Yq] = meshgrid(1:64,1:64); 

for t = 1:size(rsig,2)
    %interpolate full 8x8 Kuklik signal to 64x64 grid
    sig8x8 = reshape(rsig(:,t),[8,8]);
    rsig_int(:,:,t)=interp2(X,Y,sig8x8,Xq,Yq,'linear');
end

siglen = size(phase,2);

%compute phase via Hilbert Transform
rsig_int = reshape(rsig_int,[64*64,siglen]); % Reshape so that each pixel corresponds to one row
hV = hilbert(rsig_int'); %% Hilbert transform
rphi = nan(64*64,size(hV,1)); 
%%% Phase calculation
for i = 1:64*64
rphi(i,:) = mod(angle(hV(:,i))+pi/2,2*pi); %shift pi/2 so phase inverts at activation
end
rphi2=rphi;
rphi = reshape(rphi,[64,64,size(rphi,2)]);
sig64=rphi;

N_FRAMES_VID = min(NFRAMES_VID, size(phase,2));


%% Calculation of phase for each point in 3D
disp('Calculating 3D interpolation')
raw=[];
if numsub>0 %%% Interpolation of the Kuklik signal
    for frame = 1:NFRAMES_VID 
        phaseFrame = phase(:,frame); 
        for vertex = 1: Nvertices    
            raw(vertex, frame) = sum(barycentric(vertex,:).*[phaseFrame(faces1(ti(vertex),1)), phaseFrame(faces1(ti(vertex),2)), phaseFrame(faces1(ti(vertex),3))]);
        end
    end
end

% Hilbert transform and phase calculation of the interpolated Kuklik signal
hV = hilbert(raw');
for vertex=1:Nvertices
    C(1,vertex,:) = mod(angle(hV(:,vertex))+pi/2,2*pi)';
end

%% Calculate camera position to point at TS

basketC = mean(electrodes); 
unitV = TS - basketC; 
unitV = unitV/norm(unitV); 
camP = TS + 35.*unitV; 


%% Cumulative phase calculation

disp('Calculating phase singularities in 3D')
%%% Calculation of the phase singularities in 3D
Nfaces = max(size(fx));
phi = zeros(Nfaces,NFRAMES_VID);
for frame = 1:NFRAMES_VID
    for triangle = 1:Nfaces
        TRIplet = fx(triangle,:)';
        v1 = C(1,TRIplet(1),frame);
        v2 = C(1,TRIplet(2),frame);
        v3 = C(1,TRIplet(3),frame);
        diff_phi = NaN(3,1);
        diff_phi(1) = v3 - v2;
        diff_phi(2) = v2 - v1;
        diff_phi(3) = v1 - v3;
        for i = 1:3
            if abs(diff_phi(i)) <= pi
                phi(triangle,frame) = phi(triangle,frame) + diff_phi(i);
            elseif diff_phi(i) < -pi
                phi(triangle,frame) = phi(triangle,frame) + diff_phi(i) + 2*pi;
                if diff_phi(i) < - 3*pi
                    disp('ERROR');
                    diff_phi(i,1)
                end
            elseif diff_phi(i) > pi
                phi(triangle,frame) = phi(triangle,frame) + diff_phi(i) - 2*pi;
                if diff_phi(i) > 3*pi
                    disp('ERROR');
                    diff_phi(i)
                end
            end
        end
    end
end

disp('Calculating phase singularities in 2D')
%%% Calculation for 2D
Nfaces2=size(faces2D_ip,1);
phi2D = zeros(Nfaces2,NFRAMES_VID);
for frame = 1:NFRAMES_VID
    for triangle = 1:Nfaces2
        TRIplet = faces2D_ip(triangle,:)';
        v1 = rphi2(TRIplet(1), frame);
        v2 = rphi2(TRIplet(2), frame);
        v3 = rphi2(TRIplet(3), frame);
        diff_phi = NaN(3,1);
        diff_phi(1) = v3 - v2;
        diff_phi(2) = v2 - v1;
        diff_phi(3) = v1 - v3;
        for i = 1:3
            if abs(diff_phi(i)) <= pi
                phi2D(triangle,frame) = phi2D(triangle,frame) + diff_phi(i);
            elseif diff_phi(i) < -pi
                phi2D(triangle,frame) = phi2D(triangle,frame) + diff_phi(i) + 2*pi;
                if diff_phi(i) < - 3*pi
                    disp('ERROR');
                    diff_phi(i,1)
                end
            elseif diff_phi(i) > pi
                phi2D(triangle,frame) = phi2D(triangle,frame) + diff_phi(i) - 2*pi;
                if diff_phi(i) > 3*pi
                    disp('ERROR');
                    diff_phi(i)
                end
            end
        end
    end
end

ps_yesno=zeros(Nfaces, NFRAMES_VID); %%%Phase singularities in 3D
ps_yesno2D=zeros(Nfaces2, NFRAMES_VID); %%%Phase singularities in 2D


%% Generating maps of phase singularities per triangle
% ps_yesno format: each column represents a frame, and each row a triangle. 1 if PS, 0 otherwise
for frame = 1:NFRAMES_VID
    for triangle = 1:Nfaces
        if abs(phi(triangle,frame))==2*pi
            ps_yesno(triangle, frame)=1;
        end
    end
end

for frame = 1:NFRAMES_VID
    for triangle = 1:Nfaces2
        if abs(phi2D(triangle,frame))==2*pi
            ps_yesno2D(triangle, frame)=1;
        end
    end
end

%% Summarize PS results
ps_yesno_final=ps_yesno;

ps_yesno_triangles=ti(fx(:,1)).*ps_yesno_final; %%% Multiply each row by the ID of the triangle in which it is contained.
ps_yesno_triangles_2D=ti_plane(faces2D_ip(:,1)).*ps_yesno2D; %%% Multiply each row by the ID of the triangle in which it is contained.

% Counting of phase singularities for each triangle. Each row represents a big triangle
ps_cum=[];
for i=1:size(faces1,1) %For each triangle
    check=0;
    for frame=1:NFRAMES_VID
        if any(ps_yesno_triangles(:,frame)==i) %%% If there is one or more PS in that triangle at that time point, count 1.
            check=check+1;
        end
    end
    ps_cum(i)=check;
end

%%%Same for 2D
ps_cum2D=[];
for i=1:size(faces2D,1)
    check=0;
    for frame=1:NFRAMES_VID
        if any(ps_yesno_triangles_2D(:,frame)==i)
            check=check+1;
        end
    end
    ps_cum2D(i)=check;
end

%%% Create table comparing number of PS for each triangle in 2D and 3D
comparison=[];
for i=1:size(faces2D,1)
    for j=1:size(faces1,1)
        if isequal(faces2D(i,:),faces1(j,:)) %%% Check that we're comparing same triangle
            comparison(i, :)=[ps_cum2D(i), ps_cum(j)];    
            break
        end
    end
end

%% Save locations
disp('Calculating distance between 2D phase singularities and projected 3D')
ps_registry=[];
for frame = 1:1:NFRAMES_VID
        circumcenters = incenter(plane_subdiv, find(ps_yesno2D(:,frame)==1));

        Ntriangles = size(ps_yesno,1);
        ps_projected=[];
        for triangle = 1:Ntriangles
            if ps_yesno(triangle,frame) == 1
                vertex = incenter(Emesh_subdiv,triangle);
                [x2D,y2D] = ThreeDcoordToTwoD1(vertex,Emesh);
                plot(x2D,y2D, 'w.', 'MarkerSize', 25);
                ps_projected=[ps_projected; vertex, x2D, y2D];
            end
        end
        
        if ~isempty(ps_projected)
            ps_projected=ps_projected(~isnan(ps_projected(:,4)),:);
            if ~isempty(ps_projected)
                if ~isempty(circumcenters)
                    circumcenters = circumcenters(:, 1:2);
                    circumcenters=circumcenters(~isnan(circumcenters(:,1)),:);
                    dbstop if error
                    [k,d] = dsearchn(ps_projected(:,4:5), circumcenters);
                    for i=1:length(k)
                        ps_registry=[ps_registry; frame, circumcenters(i,:), ps_projected(k(i), 4:5), d(i)];
                    end
                else
                    for i=1:size(ps_projected,1)
                    ps_registry=[ps_registry; frame, NaN, NaN, ps_projected(i,4:5), NaN];
                    end
                end
            else
                if ~isempty(circumcenters)
                    circumcenters = circumcenters(:, 1:2);
                    for i=1:size(circumcenters,1)
                        ps_registry=[ps_registry; frame, circumcenters(i,:),NaN, NaN, NaN];
                    end
                end
            end
        end
end 



%% Plottings

%Plot correlation between triangles
correlation_fig = figure;
plot(comparison(:, 1),comparison(:,2), 'x')
title('Number of phase singularities')
xlabel('2D')
ylabel('3D')
savefig(correlation_fig,strcat(figpath,'phase_',patientnum,'_correlation.fig'))

%Normalize
ps_cum=ps_cum./max(ps_cum);
ps_cum2D=ps_cum2D./max(ps_cum2D);

%%Plot using ps_cum as color
phase_density_figure = figure;
subplot(121)

image(reshape(ps_cum2D(ti_plane), [64,64]), 'CDataMapping', 'scaled')
colormap jet
axis([1 64 1 64])
axis equal tight
set(gca,'XTick',linspace(1,64,8));
set(gca,'XTickLabel',circshift({'A','B','C','D','E','F','G','H'},[0,0]));
set(gca,'YTick',linspace(1,64,8));
set(gca,'YTickLabel',{'1','2','3','4','5','6','7','8'});
hold on;
plot(TS2D(1), TS2D(2), 'r.', 'MarkerSize', 35)
hold on;            

[d,TSp] = point2trimesh('Faces',faces1,...
                        'Vertices',electrodes,...
                        'QueryPoints',TS);
subplot(122)
trisurf(Emesh, ps_cum)
colormap jet
hold on
plot3(TSp(1), TSp(2), TSp(3), 'r.', 'MarkerSize', 35)
splineLabels = str2mat("ABCDEFGH");


for spline = 1:8
    for elec = 1:8
        elecI = ((spline-1)*8)+ elec;
        pos = electrodes(elecI, :);
        unitV = pos - basketC; 
        unitV = unitV/norm(unitV); 
        plot3(pos(1,1),pos(1,2),pos(1,3),'k. ','MarkerSize',20);
        pos = pos + 2*unitV;
        elecLabel = strcat(splineLabels(spline),num2str(elec));
        text(pos(1,1),pos(1,2),pos(1,3),elecLabel,'FontSize',12,'Color',[0.2 0.2 0.2]);
    end
end
savefig(phase_density_figure,strcat(figpath,'phase_',patientnum,'_comparison.fig'))

coefficient = corrcoef(comparison)    

ps_filtered = ps_registry;
ps_filtered(any(isnan(ps_filtered), 2), :) = [];
%%% Save data
save(strcat(datapath,'ps_num_comparison_', patientnum, '.mat'), 'comparison')

save(strcat(datapath,'raw_registry_', patientnum, '.mat'), 'ps_registry')
save(strcat(datapath,'filtered_registry_', patientnum, '.mat'), 'ps_filtered')



