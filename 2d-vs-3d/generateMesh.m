function [faces,vertices] = generateMesh(faces,vertices,DesiredNvertices,minSizeFactor,minFormFactor)

%% Generate Downsampled Mesh

BinFigure = figure('visible','off');
P = patch('Faces',faces,'Vertices',vertices);

reduceFactor = double(max(size(vertices)));
reduceFactor = double(DesiredNvertices)/reduceFactor;
P = reducepatch(P,reduceFactor);
faces = P.faces;
vertices = P.vertices;
close(BinFigure);

%% Control Mesh

elementSizeWarning = 0;
Nfaces = max(size(faces));
Nvertices = max(size(vertices));
areas = zeros(Nfaces,1);
formFactors = zeros(Nfaces,1);
for i = 1:Nfaces
    a = vertices(faces(i,1),:);
    b = vertices(faces(i,2),:);
    c = vertices(faces(i,3),:);
    ab = norm(a-b);
    ac = norm(a-c);
    bc = norm(b-c);
    p = ab + ac+ bc;
    s = p/2;
    areas(i) = sqrt(s*(s-ab)*(s-ac)*(s-bc)); %area of triangle given by Heron's formula
    f1 = ab/ac;f2 = ac/ab;f3 = ab/bc;f4 = bc/ab;f5 = ac/bc;f6 = bc/ac;
    f = [f1 f2 f3 f4 f5 f6];
    formFactors(i) = max(f);
end

sizeFactor = double(max(areas))/double(min(areas));
if sizeFactor > minSizeFactor
    elementSizeWarning = 1;
end

formFactorWarning = 0;
formFactor = double(max(formFactors))/double(min(formFactors));
if formFactor > minFormFactor
    formFactorWarning = 1;
end

%% Display Mesh Results;

disp(sprintf('\n'));
disp('Mesh succesfully completed.');
disp(strcat(sprintf('Number of elements: %i/%i',Nvertices,DesiredNvertices),'.'));
if elementSizeWarning
    disp('WARNING: Element size condition unmet.');
else
    disp('Element size condition met.');
end
sizeFactor = round(sizeFactor*100);
disp(strcat(sprintf('Size factor: %i',sizeFactor),'%.'));
if formFactorWarning
    disp('WARNING: Form factor condition unmet.');
else
    disp('Form factor condition met.');
end
formFactor = round(formFactor*100);
disp(strcat(sprintf('Form factor: %i',formFactor),'%.'));

