%2 base nodes localize the three shared reflectors and then use the 3
%shared reflectors to localize target node 

clear all

%the range and angle measurement errors are zero mean Gaussian random
%variables    
sigma_R=5;   %the range measurement error standard deviation  
variance_R=sigma_R^2;

sigma_theta=1/180*pi;     %the angle measurement error standard deviation  
variance_theta=sigma_theta^2;

theta_th=pi/10; %effective shared reflector determin inner angle threshold 

%Repeat times, i.e., the number of geometrical distribution
N=1000;

%the number of base nodes in each set, from which the two base node is
%selected. The base node position israndomly generated, thus two base nodes
%may be close to each other and it is not possiable for real application,
%thus we generate a group of base node position and then select two base
%nodes with reasonable relative position    
M=10;

%number of reflectors that the three shared reflectors would be selected
K=10;

%the number of shared reflectors to be selected
SR=3;

%area size
XY=200*sigma_R;

%base nodes are uniformly distributed in the desired area, i.e., 
%-100sigma_R < x < 100sigma_R 
%-100sigma_R < y < 100sigma_R 
x_base=XY*(rand(N,M)-0.5);
y_base=XY*(rand(N,M)-0.5);

%the reflectors are randomly distributed in the same area
x_reflector=XY*(rand(N,K)-0.5);
y_reflector=XY*(rand(N,K)-0.5);

%the targets are randomly distributed in the same area
x_target=XY*(rand(N,1)-0.5);
y_target=XY*(rand(N,1)-0.5);

%one LOS base node
x_reflector(:,1)=x_target;
y_reflector(:,1)=y_target;

for i=1:N 
    %one set of base node and one set of reflectors and one target node
    X_base=[x_base(i,:);y_base(i,:)];
    X_reflector=[x_reflector(i,:);y_reflector(i,:)];
    X_target=[x_target(i);y_target(i)];
    
    %the true distance between target node and reflectors 
    for r=1:K
        R_RT_t(r)=norm(X_target-X_reflector(:,r));        
    end
       
    %the measured distance and angle of target node with respect to
    %base nodes due to reflectors 
    for b=1:M    %base nodes
        for r=1:K  %reflector
            R_BRT_mea(b,r)=norm(X_base(:,b)-X_reflector(:,r))+...
                R_RT_t(r)+sigma_R*randn(1,1);
            theta_BRT_mea(b,r)=atan2(X_reflector(2,r)-X_base(2,b),...
                X_reflector(1,r)-X_base(1,b))+sigma_theta*randn(1,1);
            if theta_BRT_mea(b,r)<0
                theta_BRT_mea(b,r)=theta_BRT_mea(b,r)+2*pi;
            end
        end
    end

    %localize the shared reflectors, then calculate the range between
    %the shared reflector and the target node
    %get two base nodes to localize shared reflectors
    NN=2;
    
    for refl=1:SR    %reflector index
        x_reflector(1,refl);
        y_reflector(1,refl);
        
        %here is a trick, it is easier to select base node than select
        %shared reflectors. 
        %find the base nodes that share the reflector and have reasonable 
        %geometrical distribution
        bs_selected_number(refl,1)=1; %select the first base node
        
        %selected base node 
        X_BS_selected(:,1)=X_base(:,1);
        
        %the angles of the reflector with respect to the two selected base
        %nodes 
        theta_selected(1)=theta_BRT_mea(1,refl);
                
        %the range between the two selected base nodes and the target node
        %through the reflector 
        range_selected(1)=R_BRT_mea(1,refl);
            
        bs_count=1;    %count the selected base nodes
            
        %from the rest base nodes, select one 
        for bs=1:M-1     %base node
            delat_angle=abs(theta_selected(1)-theta_BRT_mea(bs+1,refl));
            if (theta_th<=delat_angle && delat_angle<pi-theta_th)||...
                    (pi+theta_th<=delat_angle && delat_angle<2*pi-theta_th)
                bs_count=bs_count+1;
                X_BS_selected(:,2)=X_base(:,bs+1);
                
                %the angles of the reflector with respect to the two
                %selected base nodes
                theta_selected(2)=theta_BRT_mea(bs+1,refl);
                
                %the range between the two selected base nodes and the
                %target node through the reflector
                range_selected(2)=R_BRT_mea(bs+1,refl);
            end
            if bs_count>=NN
                break
            end
        end
            
        %calculate the reflector position 
            
        X_B=X_BS_selected;
        Theta=theta_selected';
            
        %generating two aiding points 
        for aid=1:2
            if theta_selected(aid)==pi/2||theta_selected(aid)==3*pi/2
                x_aid(aid)=X_BS_selected(1,aid);
                y_aid(aid)=X_BS_selected(2,aid)+100;
            else
                x_aid(aid)=X_BS_selected(1,aid)+200;
                y_aid(aid)=200*tan(theta_selected(aid))+X_BS_selected(2,aid);
            end
        end

        %AX=C=DE            
        A=[y_aid'-X_BS_selected(2,:)' -(x_aid'-X_BS_selected(1,:)')];
        D=[-x_aid(1) y_aid(1) 0 0; 0 0 -x_aid(2) y_aid(2)];
        E=[X_BS_selected(2,1) X_BS_selected(1,1)...
            X_BS_selected(2,2) X_BS_selected(1,2)]';
        
        %reflector's position
        X_reflector_hat=(A'*A)^(-1)*A'*D*E;           
        x_reflector_position(refl)=X_reflector_hat(1);
        y_reflector_position(refl)=X_reflector_hat(2);
            
        %the distance between reflector and Base node
        R_BR_1r=norm(X_reflector_hat-X_BS_selected(:,1));
        R_BR_2r=norm(X_reflector_hat-X_BS_selected(:,2));
            
        %the distance between reflector and target node
        R_RT_cal(refl)=((range_selected(1)-R_BR_1r)+...
            (range_selected(2)-R_BR_2r))/2;
    end

    %target node localization
    reflector_xy=[x_reflector_position; y_reflector_position]; 
    X_T_hat=[1;10];
    %repeat of the iteration for calculating the target node positioin
    for iii=1:15   
        for jjj=1:SR  %calculate R^
            R_hat(jjj)=sqrt((X_T_hat'-...
                [x_reflector_position(jjj) y_reflector_position(jjj)])*...
                (X_T_hat'-[x_reflector_position(jjj)...
                y_reflector_position(jjj)])');
            for kkk=1:2   %for x and y
                H(jjj,kkk)=(X_T_hat(kkk)-reflector_xy(kkk,jjj))/R_hat(jjj);
            end
        end
        delta_R=R_RT_cal'-R_hat';
        
        if det(H'*H)<1e-6
            delta_X=(H'*H+diag([1e-6 1e-6]))^(-1)*H'*delta_R;
        else
            delta_X=(H'*H)^(-1)*H'*delta_R;
        end
        X_T_hat=X_T_hat+delta_X;
    end
    delta_r_3re_2bs(i)=norm(X_T_hat-X_target);
end
       

cdfplot(delta_r_3re_2bs)
set(findobj(gca,'Type','line','Color',[0 0 1]),...
    'Color','red',...
    'LineWidth',2,...
    'LineStyle',':')
xlim([0 120])
xlabel('NLOS localization error (m)')
ylabel('Localization error CDF')
legend('NLOS localization error CDF with 2 BN and 3 shared reflectors');

