function [proj,projBarycentric,TI] = projectPointOnEShell(electrodes,Eshell,QP)
Nfaces = max(size(Eshell.ConnectivityList));
basketC = mean(electrodes);
vertices = Eshell.Points;
proj = NaN(1,3);
TI = 0;
for i = 1:Nfaces
	TRIplet = Eshell.ConnectivityList(i,:)';
    [~,I] = sort(TRIplet);
	v1 = vertices(TRIplet(1),:);
	v2 = vertices(TRIplet(2),:);
	v3 = vertices(TRIplet(3),:);

	N = cross(v2-v1,v3-v1); % normal
    P = QP + dot(v1-QP,N)/dot(basketC-QP,N)*(basketC-QP); % Potential intersection point
    
    if dot(cross(v2-v1,P-v1),N)>=0 &&  dot(cross(v3-v2,P-v2),N)>=0 &&  dot(cross(v1-v3,P-v3),N)>=0 % check if P lies within tri
        if isnan(proj)
            proj = P;
            TI = i;
        elseif norm(QP-P) <= norm(QP-proj)
                proj = P;
                TI = i;
        end
    end            
end

projBarycentric = cartesianToBarycentric(Eshell,TI,proj);
projBarycentric = projBarycentric(I); % barycentric coordinates are given for points in ascending order of indices in Eshell