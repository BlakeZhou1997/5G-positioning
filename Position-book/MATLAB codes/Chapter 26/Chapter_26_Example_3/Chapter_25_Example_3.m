%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2.3: NBP algorithm for cooperative localization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:       Vladimir Savic (e-mail: vsavic05@gmail.com)
% University:   Universidad Politecnica de Madrid (UPM), Spain
% Title:        NBP algorithm for cooperative localization
% Last change:  November 19, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
tic  % start timer  
fprintf('\n');
init;	% run initialization script

%%% Start NBP algorithm:
for iter=2:Niter+1    
    fprintf('Current iteration:\t%d/%d\n',iter-1,Niter);    
    % Compute messages:
	for m=Na+1:Nn
        for n=Na+1:Nn
            if dist(m,n)>0
                % Update particles:                
                direction=exp(j*2*pi*rand(Np,1)); % random angle (0-360 degrees)                
                msg_part(:,iter,m,n)=M_part(:,iter-1,m)+(dist(m,n)+randn(Np,1)*sigma).*direction;
                % Check if the particles are within the b.box:                
                for p=1:Np
                    check=1;
                    while(check<=20 ...
                    & (real(msg_part(p,iter,m,n))<real(x_min(n)) | imag(msg_part(p,iter,m,n))<imag(x_min(n))...
                    | real(msg_part(p,iter,m,n))>real(x_max(n)) | imag(msg_part(p,iter,m,n))>imag(x_max(n))))
                        msg_part(p,iter,m,n)=M_part(p,iter-1,m)+(dist(m,n)+randn(1,1)*sigma).*exp(j*2*pi*rand(1,1));
                        check=check+1;
                    end
                end
                % Update weights:
                if iter==2 
                    msg_wei(:,iter,m,n)=M_wei(:,iter-1,m); % in first iteration, no-reweighting
                else % reweight by reminder of message update rule (i.e., we need KDE of old message)
                    msg_wei(:,iter,m,n)=M_wei(:,iter-1,m)./ ...
                    part_kde(msg_wei(:,iter-1,n,m),msg_part(:,iter-1,n,m),msg_bw(iter-1,n,m),M_part(:,iter-1,m));
                    msg_wei(:,iter,m,n)=msg_wei(:,iter,m,n)/sum(msg_wei(:,iter,m,n)); % normalize
                end
                % Resample with replacement (due to the sample depletion):
                msg_part(:,iter,m,n)=randsample(msg_part(:,iter,m,n),Np,true,msg_wei(:,iter,m,n));
                % Add small jitter (<sigma) to particles (i.e., avoid same particles);
                msg_part(:,iter,m,n)=msg_part(:,iter,m,n)+sigma*randn(1,1)+j*sigma*randn(1,1);
                msg_wei(:,iter,m,n)=ones(Np,1)/Np; % now all weights are the same
                % Update bandwidth 
                msg_bw(iter,m,n)=Np^(-1/3)*std(real(msg_part(:,iter,m,n)))+j*Np^(-1/3)*std(imag(msg_part(:,iter,m,n)));
            end
        end
    end
   
    % Compute beliefs:
	for k=Na+1:Nn        
        % Compute particles:
        temp_particles=[];
        for m=Na+1:Nn
            % collect all particles from messages
            if m~=k	& dist(m,k)>0 temp_particles=[temp_particles; msg_part(:,iter,m,k)]; end 
        end
        if length(temp_particles)>0 % if there are target-target connections
            indeces=round((length(temp_particles)-1)*rand(Np,1))+1; % choose randomly Np indices
            for n=1:Np
                M_part(n,iter,k)=temp_particles(indeces(n)); % new set of particles
            end
            % Compute weights:
            msg_product=1; msg_sum=0;
            for m=Na+1:Nn
                if m~=k & dist(m,k)>0
                    msg_product=msg_product.* ...
                    part_kde(msg_wei(:,iter,m,k),msg_part(:,iter,m,k),msg_bw(iter,m,k),M_part(:,iter,k));
                    msg_product=msg_product/sum(msg_product);
                    msg_sum=msg_sum+part_kde(msg_wei(:,iter,m,k),msg_part(:,iter,m,k),msg_bw(iter,m,k),M_part(:,iter,k));
                    msg_sum=msg_sum/sum(msg_sum);
                end
            end
            M_wei(:,iter,k)=msg_product./msg_sum; % mixture importance sampling
            M_wei(:,iter,k)=M_wei(:,iter,k)/sum(M_wei(:,iter,k));
        else % only anchor-target connections (so we can use initial particles from b. box)
            M_part(:,iter,k)=(real(x_max(k)-x_min(k))*rand(Np,1) ...
            +real(x_min(k)))+j*(imag(x_max(k)-x_min(k))*rand(Np,1)+imag(x_min(k)));
            M_wei(:,iter,k)=ones(Np,1)/Np;
        end
        % Messages from anchors:
        for n=1:Na
            if dist(n,k)>0 % it's possible that there is no connection with anchor
                M_wei(:,iter,k)=M_wei(:,iter,k).*pairwise_pot(x(n),M_part(:,iter,k),dist(n,k),R,sigma);
            end
        end
        M_wei(:,iter,k)=M_wei(:,iter,k)/sum(M_wei(:,iter,k)); % normalize
      	% Messages from 2-hop anchors (if you don't want it, set two_hop='n'):
        for m=1:Na
            if no_link_status(m,k)==1 & two_hop=='y'
                M_wei(:,iter,k)=M_wei(:,iter,k).*(1-(abs(x(m)-M_part(:,iter,k))<R)+0.00001); % also avoid zero
                M_wei(:,iter,k)=M_wei(:,iter,k)/sum(M_wei(:,iter,k)); % normalize
            end
        end
        M_wei(:,iter,k)=M_wei(:,iter,k)/sum(M_wei(:,iter,k)); % normalize
        % Resample with replacement:
        M_part(:,iter,k)=randsample(M_part(:,iter,k),Np,true,M_wei(:,iter,k));
        M_wei(:,iter,k)=ones(Np,1)/Np;
        % Compute bandwidth:                    
        M_bw(iter,k)=Np^(-1/3)*std(real(M_part(:,iter,k)))+j*Np^(-1/3)*std(imag(M_part(:,iter,k)));        
    end
        
    % Estimates:
    for k=Na+1:Nn
        x_est(iter-1,k)=M_wei(:,iter,k)'*M_part(:,iter,k); % estimated location (MMSE)
        error(iter-1,k)=abs(x_est(iter-1,k)-x(k)); % error (distance between true and estimate)        
        if iter>Niter & k==Nn % compute kde of final belief for just one node
            [X Y]=meshgrid(area,area);            
            z(:,:,iter,k)=part_kde(M_wei(:,iter,k),M_part(:,iter,k),sigma*(1+j),X+j*Y);
        end
    end
    error_rms(iter-1)=sqrt(sum(error(iter-1,Na+1:Nn).^2)/(Nn-Na)); % RMS error    
end
%%% end of NBP

%%% Print & plot:
fprintf('\n');
for k=Na+1:Nn
    fprintf('True position:\t\t\tx(%d)=%.4f%+.4fi\n',k,real(x(k)),imag(x(k))); 
    fprintf('Estimated position:\tx_est(%d)=%.4f%+.4fi\n',k,real(x_est(Niter,k)),imag(x_est(Niter,k)));
end 
fprintf('\nRMS error:\t%.4f\n',error_rms(Niter));
plotfig(x,x_est(Niter,:),Na,dist,area,no_link_status,two_hop); % plot network
if placement=='dete' figure, mesh(X,Y,z(:,:,Niter+1,Nn)); end
toc % measure time

% Clear unnecessary data
clear k kk n m p d_min dist2anchor iter bias indeces msg_product msg_sum draw_edge
clear direction temp_particles kld M_part_temp grid_cell band check part_error
%%% the end

