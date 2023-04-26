opt_value=zeros(18,52);
dist_no=zeros(18,52);
dist_full=zeros(18,52);

agents=3:20;
sub_parameter=[1:26,-1:-26];

 for y=3:20
% y=6;
    cov_theta=ones(y,y);

    for i=1:y
    cov_theta(i,i)=4;
    end 
    
     for z=1:52
%   z=20;
    if z<=26 %submodular
%     sub_parameter=(-1+sqrt(1-8*(y-1-(y-1).^2)))/2;
%     v1=sub_parameter*ones(y,1);
    v1=z*ones(y,1);
    Q=diag(v1);
    Q(1,:)=1;
    Q(:,1)=1;
    Q(1,1)=z;
    
    else %supermodular
%     sub_parameter=(-1-sqrt(1-8*(y-1-(y-1).^2)))/2; 
%     v2=sub_parameter*ones(y,1);
    v2=(26-z)*ones(y,1);
    Q=diag(v2);
    Q(1,:)=-1;
    Q(:,1)=-1;
    Q(1,1)=abs(26-z);
    end
    
    
S=zeros(2*y,2*y);
R=zeros(2*y,2*y);

V_sw=[-Q,eye(y);eye(y),zeros(y)]; %welfare

no_info_solution=[zeros(y,y),zeros(y,y);zeros(y,y),Q];
full_info_solution=[Q\cov_theta/inv(Q).',Q\cov_theta;cov_theta/inv(Q).',cov_theta];

n=y;
cvx_begin sdp quiet
variable X(2*y,2*y) symmetric
dual variable mu{y,y}

maximize (trace(V_sw*X'))

for i=1:n
    for j=1:i
        if(i==j)
            S(n+i,n+j)=1;
        mu{i,j}: S(n+i,n+j).*X(n+i,n+j) == cov_theta(i,j); %first constraint
            S(n+i,n+j)=0;
            
        end
        if (j<i)
            S(n+i,n+j)=1;
        mu{i,j}: S(n+i,n+j).*X(n+i,n+j) == cov_theta(i,j); %first constraint
       
             S(n+i,n+j)=0;
        end
    end
end

for k=1:n
    for i=1:n
        R(k,i)=Q(k,i)/2;
        R(i,k)=Q(k,i)/2;
    end
    R(k,k)=Q(k,k);
    R(k,n+k)=-0.5;
    R(n+k,k)=-0.5;
    trace(R*X')==0; % second constraint
    R=zeros(2*n,2*n);    
    
end

X>=0; %third constraint

cvx_end 

dist_full_no=sqrt(sum((full_info_solution-no_info_solution).^2,'all'));

dist_no(y-2,z)=sqrt(sum((X-no_info_solution).^2,'all'));
dist_full(y-2,z)=sqrt(sum((X-full_info_solution).^2,'all'));
   
     end
end

normalized_dist_full=dist_full./max(max(dist_full));
 
contourf(agents,sub_parameter(5:25),normalized_dist_full(:,5:25)',10)
% yi=(-1+sqrt(1-8*(agents-1-(agents-1).^2)))/2;
% xi=interp1(sub_parameter(5:22),agents,yi);
% hold on
% plot(xi,yi,'r')


