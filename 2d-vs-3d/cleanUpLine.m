function [line] = cleanUpLine(line)
%Removes spaces and/or commas at the beginning of a line (str)
    if isempty(line(line~=','))
        line = ',';
    else
            lineLength = length(line);
            while line(1) == ' ' || line(1) == ','%length(line)>2
                n = 1;
                s = 0; %stop code for the loop
                while s == 0
                    if n < lineLength
                        if (line(n) == ',' || line(n) == ' ') && (line(n+1) ~= ',' && (line(n+1) ~= ' '))
                            line = line(n+1:end);
                            s = 1;
                         else
                            n = n+1;
                         end
                    else
                        s = 1;
                    end
                end
            end 
    end   
end

