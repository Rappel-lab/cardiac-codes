function electrodes=readshadows(filepath)

fileID=fopen(filepath);
for i=1:26
line=fgetl(fileID);
end

displayed=[];
d=1;

while 1
    
    line=fgetl(fileID);
    
    if strcmp(line(1), ',')
        break
    end
    
    C=textscan(line, '%s%d8', 'delimiter', ',');
    C{2}(1);
    if C{2}(1)
        displayed=[displayed, d];
    end
    d=d+1;
end

line=fgetl(fileID);


for i=1:displayed(1)-1
    line=fgetl(fileID);
end

basket_pos=[];

for i=1:8
line=fgetl(fileID);
C=textscan(line, '%d8%s%f%f%f%f%f%f%f', 8,'delimiter', ',');
basket_pos=[basket_pos; C{6}, C{7}, C{8}];
end
electrodes=basket_pos;



