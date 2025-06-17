function [Sigma,Epsilon,nodecount] = Stress_computation(nodes,elems,D,uh,uv)
%computes the matrices element by element, they are then assembled in the
%global matrix in the main code

%% inputs
[nnod,ndim]=size(nodes);
[nel,nnodel]=size(elems);
ndofel=ndim*nnodel;
nodecount=zeros(nnod,1);


%% displacements
uhx=uh(1:2:(2*nnod-1));
uhy=uh(2:2:2*nnod);
%% stress points
npts=size(uv,2);
dNdu = zeros(ndim,nnodel,npts);
Sigma=zeros(nel,3,npts);
Epsilon=Sigma;

%% shape function derivatives
U=uv(1,:);
V=uv(2,:);
W=1-U-V;
for ipts=1:npts %ip
        dNdu(1,1,ipts) = 1-4*W(ipts);  %dN_dxi(1,3) is identically 0
        dNdu(1,2,ipts) =  4*U(ipts)-1;
        dNdu(1,3,ipts) =  0;
        dNdu(1,4,ipts) =  4*V(ipts);
        dNdu(1,5,ipts) = -4*V(ipts);
        dNdu(1,6,ipts)= 4*(W(ipts)-U(ipts));
        dNdu(2,1,ipts) = 1-4*W(ipts);  %dN_dxi(2,2) is identically 0
        dNdu(2,2,ipts) =  0;
        dNdu(2,3,ipts) =  4*V(ipts)-1;
        dNdu(2,4,ipts) =  4*U(ipts);
        dNdu(2,5,ipts) = 4*(W(ipts)-V(ipts));
        dNdu(2,6,ipts)= -4*U(ipts);
    
end

%% element loop
for iel=1:nel
    index_el=elems(iel,:);
    nodecount(index_el,1)=nodecount(index_el,1)+ones(6,1);
    nodes_el = nodes(elems(iel,:),:);
    B_el = zeros(3,ndofel);  
    u_el=[uhx(index_el(1));uhy(index_el(1));uhx(index_el(2));uhy(index_el(2));uhx(index_el(3));uhy(index_el(3));...
        uhx(index_el(4));uhy(index_el(4));uhx(index_el(5));uhy(index_el(5));uhx(index_el(6));uhy(index_el(6))]; %change it, less explicit, not uhx and uhy but indexes
    for ipts=1:npts
        dxdu=dNdu(:,:,ipts)*nodes_el;
        dudx = inv(dxdu);
        Jac = det(dxdu); %jacobian determinant
        dNdx = dudx*dNdu(:,:,ipts);
        B_el([1,3],1:ndim:end-1) = dNdx; %use ndim, introduce ndofel
        B_el([3,2],2:ndim:end) = dNdx;
    Sigma(iel,:,ipts)=D*B_el*u_el;
    Epsilon(iel,:,ipts)=B_el*u_el;
    end
    
end