
clear all

A=[10, 31;17, 45;-17, -23;-10, -9;19, 49]/100;
X=[0 0];
sigma=1;
% Lp=[0 1;0 -1]/10;


       
        
        
    

            for i=1:5
           D(i)= norm(X-A(i,:))+ random('normal',0,sigma/100);
       end
   
    
    
     mpol x y
    
     
      for j=3:5
    G=0;
    for i=1:j
%         
          
           
       G=G+((x-A(i,1))^2+(y-A(i,2))^2-D(i)^2)^2;
       
   end
    
 
%     G=G+((delta(i)^2-norm(p(i,:))^2+c2*p(i,1)*x+2*p(i,2)*y+2*delta(i)*D)^2);
   
P=msdp(min(G)); 
[status,obj]=msol(P);
CP=[];

if status==1
    cp{j}=double([x y])*100;
end
% for i=1:size(Ap,1)
      end

    



    
    
    