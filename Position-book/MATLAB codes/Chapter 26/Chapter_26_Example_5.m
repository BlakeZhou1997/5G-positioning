%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 3.2: Spanning tree formation using BFS method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
Nn=100; % number of nodes
R=6; % transmission radius
area=[-10 10]; % deployment area
% Warning: if you change parameters (especially, R), make sure that you have 
% connected graph. In case of disconnected graph (e.g., 2 subgraphs), 
% you will get an error.

x_stat=ones(Nn,1); % status of all nodes (1 - not in span. tree, 0 - in span. tree);
root=ceil(Nn*rand(1,1)); % randomly choose a root node
current_root=root; % set current root for the first iteration

% True position of the nodes (random):
for k=1:Nn
    x(k)=((max(area)-min(area))*rand(1,1)+min(area))...
        +j*((max(area)-min(area))*rand(1,1)+min(area));
end
% Define distances (edges in graph):
number_of_neigh_of_x=zeros(Nn,1)'; neigh_of_x(:,1)=zeros(Nn,1)'; 
for m=1:Nn
    i=1;
    for n=1:Nn
        if abs(x(m)-x(n))<=R 
            dist(m,n)=abs(x(m)-x(n));
             if m~=n
                 neigh_of_x(m,i)=n;
                 number_of_neigh_of_x(m)=i; i=i+1;                
             end
        else
            dist(m,n)=-1;  % unobserved distance
        end
        dist_span(m,n)=-1; % initially empty span. tree
    end        
end

%%%%% Breadth First Search (BFS) method %%%%%%%
x_stat(root)=0; % root is already in span. tree
n=1; m=1;
while sum(x_stat)>0
    for n1=1:number_of_neigh_of_x(root) % for all neighbors of root
        m1=neigh_of_x(root,n1); % index of neighbor
        if x_stat(m1)==1 % if node m1 not explored
            x_stat(m1)=0; % mark as explored
            root_queue(n)=m1; n=n+1; % add to root queue
            % Copy distances from original graph
            dist_span(m1,root)=dist(m1,root); dist_span(root,m1)=dist(root,m1);
        end
    end
    root=root_queue(m); % define next root
    m=m+1; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot original graph:
figure,
for k=1:Nn
    plot(real(x(k)),imag(x(k)),'ok','MarkerFaceColor','k'); % plot unknown nodes
    if k==1 hold; end
end
% Draw edges
for m=1:Nn
    for n=1:Nn
        if dist(m,n)>0 & m<n
            plot([real(x(m));real(x(n))],[imag(x(m));imag(x(n))],'-b','LineWidth',1);
        end
    end
end
axis([min(area)-0.2 max(area)+0.2 min(area)-0.2 max(area)+0.2]); hold off

% Plot spanning tree:
figure,
for k=1:Nn
    plot(real(x(k)),imag(x(k)),'ok','MarkerFaceColor','k'); % plot unknown nodes
    if k==1 hold; end
end
% Draw edges
for m=1:Nn
    for n=1:Nn
        if dist_span(m,n)>0 & m<n
            plot([real(x(m));real(x(n))],[imag(x(m));imag(x(n))],'-b','LineWidth',1);
        end
    end
end
axis([min(area)-0.2 max(area)+0.2 min(area)-0.2 max(area)+0.2]); hold off

