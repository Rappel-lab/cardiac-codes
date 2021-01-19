function [objects] = xmlCoordExtract(filepath,N)
%Extracts mesh data from xml file. Vertices for object i are stored in Objecti, Polygon
%triplets for the mesh are stored in Polygonsi. N indicates which object is
%to be extracted; if N = 0, all objects are extracted. Validity of N is not
%checked!
%% Run through file to count number of objects

Nobj = 0;
fID = fopen(filepath);

line=fgetl(fID);
while isempty(strfind(line,'</DIF>'))
    if not(isempty(strfind(line,'<Vertices')))
        Nobj = Nobj+1;
    end
    line=fgetl(fID);
end
fclose(fID);

%% Extract object vertices

fID = fopen(filepath);
objects = struct('NumberOfObjects',Nobj);

for i = 1:Nobj
    xyz = zeros(1,3);
    
    while isempty(strfind(line,'<Vertices'))
        line=fgetl(fID);
    end
    line=fgetl(fID);
    if N == 0 || i == N
        l = 1;
        while isempty(strfind(line,'Vertices'))
            line = cleanUpLine(line); %delete spaces at beginning of line
            coord = textscan(line,'%f','delimiter','  ');
            coord = coord{1};

            xyz(l,:) = [coord(1) coord(3) coord(5)];

            l = l+1;
            line=fgetl(fID);
        end
        objects = setfield(objects,sprintf('Object%i',i),xyz);
    end
end
fclose(fID);   

%% Extract polygon triplets

fID = fopen(filepath);

for i = 1:Nobj
    xyz = zeros(1,3);
    
    while isempty(strfind(line,'<Polygons'))
        line=fgetl(fID);
    end
    line=fgetl(fID);
    if N == 0 || i == N
        l = 1;
        while isempty(strfind(line,'Polygons'))
            line = cleanUpLine(line); %delete spaces at beginning of line
            coord = textscan(line,'%f','delimiter','  ');
            coord = coord{1};
            xyz(l,:) = [coord(1) coord(2) coord(3)];

            l = l+1;
            line=fgetl(fID);
        end
        objects = setfield(objects,sprintf('Polygons%i',i),xyz);
    end
end
fclose(fID); 