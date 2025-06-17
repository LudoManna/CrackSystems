function [K,f] = k_el_hoop(nodes,elems,D,bx,by)
%computes the matrices element by element, they are then assembled in the
%global matrix in the main code
%% inputs
[~,ndim]=size(nodes);
[nel,nnodel]=size(elems);
ndofel=ndim*nnodel;

%% initializations
K=zeros(ndofel^2, nel);
f=zeros(ndofel,nel);

%% integration points
nip=3;
N      = zeros(nnodel,nip);
dNdu = zeros(ndim,nnodel,nip);
dNdx  = zeros(ndim,nnodel);

%% shape functions
U=[0.5;0.5;0];
V=[0.5;0;0.5];
W=1-U-V;
gw=1/2*[1/3 1/3 1/3];
for i=1:nip %ip

        N(1,i) = W(i)*(2*W(i)-1);
        N(2,i) = U(i)*(2*U(i)-1);
        N(3,i) = V(i)*(2*V(i)-1);
        N(4,i) = 4*U(i)*V(i); %nonzero shape functions. N(1) N(2) and N(3)
        N(5,i) = 4*W(i)*V(i); %N(3) are 0 on the integration points, so
        N(6,i) = 4*U(i)*W(i); %i don't even compute them
        dNdu(1,1,i) = 1-4*W(i);  %dN_dxi(1,3) is identically 0
        dNdu(1,2,i) =  4*U(i)-1;
        dNdu(1,3,i) =  0;
        dNdu(1,4,i) =  4*V(i);
        dNdu(1,5,i) = -4*V(i);
        dNdu(1,6,i)= 4*(W(i)-U(i));
        dNdu(2,1,i) = 1-4*W(i);  %dN_dxi(2,2) is identically 0
        dNdu(2,2,i) =  0;
        dNdu(2,3,i) =  4*V(i)-1;
        dNdu(2,4,i) =  4*U(i);
        dNdu(2,5,i) = 4*(W(i)-V(i));
        dNdu(2,6,i)= -4*U(i);
    
end

%% element loop
for iel=1:nel
    nodes_el = nodes(elems(iel,:),:);
    B_el = zeros(3,ndofel);
    K_el = zeros(ndofel);  %ndofel
    f_el = zeros(ndofel,1);
    
    for i=1:nip
        dxdu=dNdu(:,:,i)*nodes_el;
        dudx = inv(dxdu);
        Jac = det(dxdu); %jacobian determinant
        dNdx = dudx*dNdu(:,:,i);
        gwt = Jac*gw(i);
        B_el([1,3],1:ndim:end-1) = dNdx; %use ndim, introduce ndofel
        B_el([3,2],2:ndim:end) = dNdx;
        K_el = K_el + gwt*B_el'*D*B_el;
        f_el(1:ndim:end-1) = f_el(1:ndim:end-1) + gwt*bx*N(:,i);
        f_el(2:ndim:end) = f_el(2:ndim:end) + gwt*by*N(:,i);
         
    end
    K(:,iel)=K_el(:);
    f(:,iel)=f_el(:);
    
end