act_cent=zeros(17,1);
act_perip=zeros(17,1);
n=4;
% u_central=zeros(30,1000,17);
% u_perip=zeros(30,1000,17);
% u_central_no=zeros(30,1000,17);
% u_perip_no=zeros(30,1000,17);
 
% mu=normrnd(1,0.3,30,1);
% mu=sort(mu,'descend');
% gamma=zeros(30,1000);

% for i=1:30
%     gamma(i,:)=normrnd(mu(i),0.1,1000,1);
%  end

for k=1:30

    for j=1:1000

         for i=1:17
         beta=i+3;
%         for i=1:17
%         beta=(21-i);

v=beta*ones(n,1);
h=diag(v);
h(1,:)=1;
h(:,1)=1;
h(1,1)=beta;
inv(h);
x=sum(inv(h),2);


act_cent(i)=x(1);
act_perip(i)=x(2);

u_central(k,j,i)=-beta*(act_cent(i)*gamma(j))^2-2*((n-1)*act_cent(i)*act_perip(i))*gamma(j)^2+2*act_cent(i)*gamma(j)^2;
u_perip(k,j,i)=-beta*(act_perip(i)*gamma(j))^2-2*(act_cent(i)*act_perip(i))*gamma(j)^2+2*act_perip(i)*gamma(j)^2;

u_central_no(k,j,i)=-beta*(act_cent(i)*mu(k))^2-2*((n-1)*act_cent(i)*act_perip(i))*mu(k)^2+2*act_cent(i)*gamma(j)*mu(k);
u_perip_no(k,j,i)=-beta*(act_perip(i)*mu(k))^2-2*(act_cent(i)*act_perip(i))*mu(k)^2+2*act_perip(i)*gamma(j)*mu(k);

        end
    end
end

delta_central=u_central-u_central_no;
delta_perip=u_perip-u_perip_no;
% plot(3:20,act_cent,'r',3:20,act_perip,'b')
% plot(3:20,u_central,'r',3:20,u_perip,'b')
% hold on
% plot(3:20,u_central_no,'g',3:20,u_perip_no,'m')
% hold off

% contourf(4:20,gamma,delta_central')
% contourf(4:20,gamma,delta_perip')

% histogram(mean(delta_central(:,:,1),2)./max(mean(delta_central(:,:,1),2)),10);
% hold on
% for i=2:17
% histogram(mean(delta_central(:,:,i)./max(mean(delta_central(:,:,1),2)),2),10);
% end
% hold off
% mu=round(mu,2);
%  heatmap(4:20,mu,squeeze(mean(delta_central,2))./max(max(squeeze(mean(delta_central,2))))) %sub
%  heatmap(4:20,mu,squeeze(mean(delta_perip,2))./max(max(squeeze(mean(delta_perip,2))))) %sub
%  heatmap(-20:-4,mu,squeeze(mean(delta_central,2))./max(max(squeeze(mean(delta_central,2))))) %super
%  heatmap(-20:-4,mu,squeeze(mean(delta_perip,2))./max(max(squeeze(mean(delta_perip,2))))) %super
 y1=squeeze(mean(delta_central,2));
 y2=squeeze(mean(delta_perip,2));
t=[1 6 11 16 21 26 30];
for i=1:7
plot(4:20,y1(t(i),:),'LineWidth',2) %central
hold on
end

% for i=1:7
% plot(4:20,y2(i,:),'LineWidth',2)%peripheral
% hold on
% end

legend('$\mu=2.0735$','$\mu=1.4469$','$\mu=1.2181$','$\mu=1.1613$','$\mu=0.9811$','$\mu=0.7638$','$\mu=0.3223$')

xlabel('Submodularity parameter, $\rho$')
ylabel('Mean of $\Delta u_{i}(a, \gamma)$')