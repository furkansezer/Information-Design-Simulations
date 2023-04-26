n=5; %number of agents
rng(0,'twister');
J=30; %number of samples


beta_f=0.95; %confidence levels
beta_l=0.95;

nu_f=((1-beta_f)*J)^(-1);
nu_l=((1-beta_l)*J)^(-1);

delta=0.1:0.1:0.9; 

h_variance=0.1:0.1:1;
% h_variance=0.1:0.3:2.8;
% h_variance=0.5;

% app_factor=0.1;

solution=zeros(length(h_variance),length(delta),2*n,2*n);

opt_value=zeros(length(h_variance),length(delta),1);
S=zeros(2*n,2*n);
R=zeros(2*n,2*n);

coefficient_variance=zeros(2*n,2*n);

cov_theta=sqrt(5)*ones(n,n);

for i=1:n
    cov_theta(i,i)=5;
end 

    
no_info_solution=[zeros(n,n),zeros(n,n);zeros(n,n),cov_theta];
% full_info_solution=[Q\cov_theta/inv(Q).',Q\cov_theta;cov_theta/inv(Q).',cov_theta];

% V_ss=zeros(2*n,2*n); %agreement
%     for i=1:n
%         for j=1:n
%             V_ss(i,j)= 1/n;
%         end 
%             V_ss(i,i)=-1 + 1/n;
%     end


    
for g=1:length(h_variance)
    
    for w=1:length(delta)
    
    Q=-delta(w)*ones(n);
        for i=1:n
            Q(i,i)=20;
        end
         
    T=delta(w)*h_variance(g).*randn(n,n,J);
    
    for i=1:n
        T(i,i)=T(i,i)/delta(w);
    end
    
    Q=T+Q;

V_sw=zeros(2*n,2*n,J);
for b=1:J
    
V_sw(:,:,b)=[squeeze(-Q(:,:,b)) delta(w)*eye(n);delta(w)*eye(n) zeros(n)];    %social welfare    

end


cvx_begin quiet
variable X(2*n,2*n) symmetric semidefinite
variable t
variable rho_f
variable rho_l(n)
variable z_f(J) nonnegative
variable z_l(n,J) nonnegative


maximize (t)
subject to

for b=1:J
    
rho_f+ nu_f*sum(z_f)<=0; %social welfare obj

z_f(b)>=t-trace(squeeze(V_sw(:,:,b))*X')-rho_f; %social welfare obj


for k=1:n
    for i=1:n
        R(k,i)=Q(k,i,b)/2;
        R(i,k)=Q(i,k,b)/2;

                
    end
    R(k,k)=Q(k,k,b);
    R(k,n+k)=-0.5;
    R(n+k,k)=-0.5;
    
    
    rho_l(k)+nu_l*sum(z_l(k,b))<=0;            % equilibrium constraint  
    z_l(k,b)>=trace(R*X')-rho_l(k);            % equilibrium constraint   
    

    R=zeros(2*n,2*n);
   
end

end


for i=1:n
    for j=1:i
        if(i==j)
            S(n+i,n+j)=1;
            S(n+i,n+j).*X(n+i,n+j) == cov_theta(i,j); %assignment constraint
            S(n+i,n+j)=0;
            
        end
        if (j<i)
            S(n+i,n+j)=1;
           S(n+i,n+j).*X(n+i,n+j) == cov_theta(i,j); %assignment constraint
       
             S(n+i,n+j)=0;
        end
    end
end

cvx_end 

opt_value(g,w)=cvx_optval;
solution(g,w,:,:)=X;

    end 
    g
end

for i=1:n
    for j=1:n
        if i==j
            X(n+1,1)==X(n+i,j);
            X(n+1,1)==X(i,n+j);
        end
        if i~=j
            X(n+1,2)==X(n+i,j); 
            X(n+1,2)==X(n+j,i); 
        end    
    end
end    

dist_no=zeros(length(h_variance),length(delta));
% dist_full=zeros(10,10);
% dist_full_no=sqrt(sum((full_info_solution-no_info_solution).^2,'all'));

for i=1:length(h_variance)
    for j=1:length(delta)
        dist_no(i,j)=sqrt(sum((squeeze(solution(i,j,:,:))-no_info_solution).^2,'all'));
%         dist_full(i,j)=sqrt(sum((squeeze(solution(i,j,:,:))-full_info_solution).^2,'all'));
    end    
end    
 
normalized_dist_no=dist_no./max(dist_no, [], 'all'); %divide by max element in dist no
contourf(delta,h_variance,normalized_dist_no)

% contourf(delta,h_variance,opt_value)
 

    