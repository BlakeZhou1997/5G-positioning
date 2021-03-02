function plotfig(x,x_est,Na,dist,area,no_link_status,two_hop)
figure;
% x - all nodes, Na - number of anchors, 
% dist - ditance matrix, area - deployment area, 
% no_link_status - equal to 1 if the nodes are 2-hop neighbors
% two_hop - include messages from 2-hop neighbors or not

% Draw true positions of the nodes:
for m=1:length(x)
    if m>Na % targets
        plot(real(x(m)),imag(x(m)),'ok','LineWidth',1,'MarkerFaceColor','k');
        text(real(x(m))+0.02*(max(area)-min(area)),imag(x(m))+0.02*(max(area)-min(area)),int2str(m));         
    else % anchors
        plot(real(x(m)),imag(x(m)),'sr','LineWidth',1,'MarkerFaceColor','r'); if m==1 hold; end 
        text(real(x(m))+0.02*(max(area)-min(area)),imag(x(m))+0.02*(max(area)-min(area)),int2str(m));        
    end
end
% Draw edges:
for m=1:length(x)
    for n=1:length(x)
        if dist(m,n)>0 & m<n
            if (m>Na)||(n>Na)
                plot([real(x(m));real(x(n))],[imag(x(m));imag(x(n))],'-b','LineWidth',1);
            end
        elseif no_link_status(m,n)>0 & m<n & two_hop=='y'
            if (m<=Na)&(n>Na)
                plot([real(x(m));real(x(n))],[imag(x(m));imag(x(n))],'--b','LineWidth',1);
            end            
        end
    end
end
% Plot estimates:
for k=Na+1:length(x)
    plot(x_est(k),'.k','LineWidth',1,'MarkerFaceColor','k');
    plot([real(x_est(k));real(x(k))],[imag(x_est(k));imag(x(k))],'-k','LineWidth',1);        
end
axis([0.97*min(area) 1.02*max(area) 0.97*min(area) 1.02*max(area)])
hold off


