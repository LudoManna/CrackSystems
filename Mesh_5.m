clc
clear

tic
RAND_DEG=[0,0.25,0.5,1];
ALPHA_0=[0,pi/4,pi/2,3*pi/4];

%% input of geometric parameters
for irand=1:1
for ialpha=1:1
    
    clearvars -except ialpha irand RAND_DEG ALPHA_0
    
n_ell=100;
alpha_0=ALPHA_0(ialpha);
rand_deg=RAND_DEG(irand);

rm=100;
rM=10000;

am=0.1;
aM=0.2;

n_crack=200; %500

rat=rm+rand(n_crack,1)*(rM-rm);
a=am+rand(n_crack,1)*(aM-am);
b=a./rat;
alpha=alpha_0-pi/2*rand_deg+rand_deg*pi*rand(n_crack,1);

% fr=0.5;
fr=0.5;
% box_a=fr*sqrt(n_crack*2*max(a));
box_a=7;
Dis=box_a/1000;
xc=zeros(n_crack,1);
yc=xc;

bound_lim=1;

xc(1)=-box_a+bound_lim*max(a)+2*(box_a-bound_lim*max(a))*rand;
    yc(1)=-box_a+bound_lim*max(a)+2*(box_a-bound_lim*max(a))*rand;

for i=2:n_crack
    xc(i)=-box_a+bound_lim*max(a)+2*(box_a-bound_lim*max(a))*rand;
    yc(i)=-box_a+bound_lim*max(a)+2*(box_a-bound_lim*max(a))*rand;
    while min(sqrt((xc(i)-xc(1:i-1)).^2+(yc(i)-yc(1:i-1)).^2))<=(a(i)+a(i-1))
        xc(i)=-box_a+bound_lim*max(a)+2*(box_a-bound_lim*max(a))*rand;
    yc(i)=-box_a+bound_lim*max(a)+2*(box_a-bound_lim*max(a))*rand;
    end
end

rc=[xc,yc];

fr1=1;
box_a=fr1*box_a;

%% mesh

[nodes,elems]=MeshM5(rc,a,b,n_ell,alpha,n_crack,box_a); %mesh
[nnod,ndim]=size(nodes);
[nel,nnodel]=size(elems);

figure(1)
hold on
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*nodes(:,2))
view(2)
axis equal
xlim([-box_a box_a]/fr1)
ylim([-box_a box_a]/fr1)
shading interp
plot([nodes(1:4,1);nodes(1,1)],[nodes(1:4,2);nodes(1,2)],'k')
for ip=1:n_crack
    plot([nodes(5+(ip-1)*n_ell:4+ip*n_ell,1);nodes(5+(ip-1)*n_ell,1)], ...
       [nodes(5+(ip-1)*n_ell:4+ip*n_ell,2);nodes(5+(ip-1)*n_ell,2)],'k')
end
xlabel('x(m)')
ylabel('y(m)')

ndofel=ndim*nnodel;
index_bottom=find(nodes(:,2)==-box_a);
index_left=find(nodes(:,1)==-box_a);
nbot=length(index_bottom);
nleft=length(index_left);  
index_top=find(nodes(:,2)==box_a);
index_right=find(nodes(:,1)==box_a); 
ntop=length(index_top);
nright=length(index_right);
%READ ABOUT STRUCTURES


%% material parameters
E=6*10^10; %Tal
nu=0.25;
mu=E/(2*(1+nu));
lambda=(nu*E)/((1+nu)*(1-2*nu));

% mu=1;                       %lamè constants
% nu=0.35;
% lambda=2*mu*nu/(1-2*nu);
D=[2*mu+lambda, lambda, 0; lambda, 2*mu+lambda, 0; 0, 0, mu]; 

%% body forces
bx=0; 
by=0;

toc
tic

%% FEM stress matrix and external force computation

[K,f] = k_el_hoop(nodes,elems,D,bx,by);

toc
%% assembly
tic
index_gl = 2*double(elems(:,[1,1,2,2,3,3,4,4,5,5,6,6]))...
          -[1,0,1,0,1,0,1,0,1,0,1,0];
index_gl=index_gl';
I=repmat(index_gl,ndofel,1);
J=reshape(repmat(index_gl(:)',ndofel,1),ndofel^2,[]);
k_gl=sparse(I,J,K);
ELMS=zeros(1,ndofel*nel);
elems=elems';
ELMS(2:2:end)=2*elems(:)';
ELMS(1:2:end-1)=2*elems(:)'-1;
f_gl=accumarray(ELMS',f(:));
elems=elems';

%% neumann bc


%% dirichlet boundary conditions

k_gl(1,:)=zeros(1,ndim*nnod);
k_gl(1,1)=eye(1);
f_gl(1)=0;

k_gl(2*index_top,:)=zeros(ntop,ndim*nnod);
k_gl(2*index_top,2*index_top)=eye(ntop,ntop);
f_gl(2*index_top)=0;

k_gl(2*index_bottom,:)=zeros(nbot,ndim*nnod);
k_gl(2*index_bottom,2*index_bottom)=eye(nbot,nbot);
f_gl(2*index_bottom)=0;

k_gl(2*index_left,:)=zeros(nleft,ndim*nnod);
k_gl(2*index_left,2*index_left)=eye(nleft,nleft);
f_gl(2*index_left)=0;

k_gl(2*index_right,:)=zeros(nright,ndim*nnod);
k_gl(2*index_right,2*index_right)=eye(nright,nright);
f_gl(2*index_right)=0;

k_gl(2*index_top-1,:)=zeros(ntop,ndim*nnod);
k_gl(2*index_top-1,2*index_top-1)=eye(ntop,ntop);
f_gl(2*index_top-1)=Dis;

k_gl(2*index_bottom-1,:)=zeros(nbot,ndim*nnod);
k_gl(2*index_bottom-1,2*index_bottom-1)=eye(nbot,nbot);
f_gl(2*index_bottom-1)=-Dis;


toc
%% displacement computation
tic
uh=k_gl\f_gl; %displacements (velocities in this case)

toc
tic
uh=full(uh); %without this i had some problems with plotting the displacements
uhx=uh(1:ndim:(2*nnod-1));
uhy=uh(2:ndim:2*nnod);
nodes1=nodes+[uhx,uhy];
X=nodes(:,1);
Y=nodes(:,2);
u=sqrt(uhx.^2+uhy.^2); %displacement modulus at each nodal point

%% stress tensor element by element

uv=[0 0; 1 0; 0 1]';
Sigma = Stress_computation(nodes,elems,D,uh,uv);
Sigma=permute(Sigma,[1 3 2]);
Sigma_xx=Sigma(:,:,1);
Sigma_yy=Sigma(:,:,2);
Sigma_xy=Sigma(:,:,3);
Sigma_diff=sqrt((Sigma_xx+Sigma_yy).^2-4*(Sigma_xx.*Sigma_yy-Sigma_xy.^2));
Sigma_m=(Sigma_xx+Sigma_yy)/2;



%% principal stresses
xc=zeros(nel,1);
yc=zeros(nel,1);
sig_xx=zeros(nel,1);
sig_xy=zeros(nel,1);
sig_yy=zeros(nel,1);
sigma1=zeros(nel,1);
sigma2=zeros(nel,1);

for iel=1:nel
    
    xc(iel)=mean(nodes(elems(iel,:),1));
    yc(iel)=mean(nodes(elems(iel,:),2));
    sig_xx(iel)=mean(Sigma_xx(iel,:));
    sig_xy(iel)=mean(Sigma_xy(iel,:));
    sig_yy(iel)=mean(Sigma_yy(iel,:));
    sig=[sig_xx(iel) sig_xy(iel); sig_xy(iel) sig_yy(iel)];
    [eigenvec,eigenval]= eig(sig);
    str_axe1(iel,:)=eigenvec(:,1)';
    str_axe2(iel,:)=eigenvec(:,2)';
    sigma1(iel)=-eigenval(1,1);
    sigma2(iel)=-eigenval(2,2);
end

sigma_diff=sigma1-sigma2;
sigma_m=(sigma1+sigma2)/2;

%%

% M=max(max(Sigma_diff))
% m=min(min(Sigma_diff))
% % SD1(icl)=M;
% [T,WS,stats]=tsearch2(nodes',elems(:,1:3)',[nodes(Loc(n_crack),1)+10^-9; ...
%     nodes(Loc(n_crack),2)]);
% T_node=find(elems(T,1:3)'==Loc(n_crack));
% SD1(icl)=Sigma_diff(T,T_node);
%% points near the tip
toc
tic

[sigma_xx,sigma_yy,sigma_xy,sigma,sigma_1,sigma_2,sigmavec_1,sigmavec_2] ...
    =Nodebynodetri6(nnod,elems,Sigma_xx,Sigma_yy,Sigma_xy);

%% stress intensity factor
% 
% KI(:,iA,iB)=sqrt(2*pi*(nodes(P_label,1)-a)).*sigma_yy(P_label);
% KII(:,iA,iB)=sqrt(2*pi*(nodes(P_label,1)-a)).*sigma_xy(P_label);
% 

toc

IAR=ialpha+(irand-1)*length(ALPHA_0);

tic
% save(string(IAR))
toc

end

end
