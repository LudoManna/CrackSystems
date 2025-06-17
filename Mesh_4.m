clc
clear
DIST=[1/100 1/10 1 10];
% DIST=[1/100];
% ALPHA=[-pi/3 -pi/4 -pi/6 0 pi/6 pi/4 pi/3 pi/2];
% ALPHA=[-pi/12];
% ALPHA=[-146*pi/180];
ALPHA=[-pi/3 -pi/6 0];


for idist=1:1
for ialp=1:1
 
    clearvars -except ialp idist DIST ALPHA
    
    
tic

% fi=0.00003; %porosity, ratio between void area and total area
box_a=60; %horizontal half-dimension of the domain box
box_b=60; %horizontal half-dimension of the domain box
rat=1000;
Dis=box_a/1000;
nh=8;
no=4;

% d=0.1+0.5*rand(nh-1,1);
% d=0.2*ones(nh-1,1);
a=0.005;
d=a*DIST(idist)*ones(nh-1,1);
b=a/rat;
% alpha=pi/12;
n_ell=1000;
L=2*nh*a+sum(d);

rc=zeros(nh,2);

for i=1:nh
    rc(i,:)=[-L/2+a+(i-1)*2*a+sum(d(1:i-1,1)),0];    
end

alpha=ALPHA(ialp);


%% mesh

[nodes,elems]=MeshM4(rc,a,b,n_ell,nh,box_a,box_b,alpha); %mesh

L=0.1150;

figure(2)
hold on
trisurf(elems(:,1:3),nodes(:,1),nodes(:,2),0*nodes(:,2))
plot([-0.6*L 0.6*L 0.6*L -0.6*L -0.6*L],[-L/4 -L/4 L/4 L/4 -L/4],'k','LineWidth',2)
for ic=1:nh
    plot([nodes(5+(ic-1)*n_ell:4+ic*n_ell,1);nodes(5+(ic-1)*n_ell,1)], ...
        [nodes(5+(ic-1)*n_ell:4+ic*n_ell,2);nodes(5+(ic-1)*n_ell,2)],'k','LineWidth',2)
end
view(2)
axis equal
xlim([-0.6*L 0.6*L])
ylim([-L/4 L/4])
shading interp
set(gca,'xtick',[])
set(gca,'ytick',[])

[nnod,ndim]=size(nodes);
[nel,nnodel]=size(elems);

ndofel=ndim*nnodel;
index_bottom=find(nodes(:,2)==-box_b);
index_left=find(nodes(:,1)==-box_a);
nbot=length(index_bottom);
nleft=length(index_left);  
index_top=find(nodes(:,2)==box_b);
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

IAD=ialp+(idist-1)*length(ALPHA);

tic
save(string(IAD))
% filename='56dist001';
% save(filename)
toc

end
end

