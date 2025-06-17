function [nodes,elems]=MeshM4_EE2(rc,a,b,n_ell,nh,box_a,box_b,alpha,xr,yr)


theta=linspace(0,2*pi,n_ell+1);
theta(end)=[];

cont=[-box_a, -box_b; box_a, -box_b; box_a, box_b; -box_a, box_b]';

points1=zeros(nh*n_ell,2);
for i=1:nh
    points1((i-1)*n_ell+1:i*n_ell,:)=[a*cos(theta'),b*sin(theta')];
end

points1=points1*[cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
points1=points1+repelem(rc,n_ell,1);

elim=[];
   for j=1:nh
       elim=[elim;find(((xr-rc(j,1))*cos(-alpha)+(yr-rc(j,2))*sin(-alpha)).^2/a^2 ...
           +((xr-rc(j,1))*sin(-alpha)-(yr-rc(j,2))*cos(-alpha)).^2/b^2 <1)];
   end
elim=unique(elim);
xr(elim)=[];
yr(elim)=[];

points=[cont, points1',[xr';yr']];
segments_box=[1,2,3,4; 2,3,4,1];


segments_ellipse=4+[1:n_ell;2:n_ell,1];
seg_ELL=repmat(segments_ellipse,1,nh-1)+n_ell*repmat(repelem(1:nh-1,n_ell),2,1);
segments= [segments_box,segments_ellipse,seg_ELL];

% figure(1)
% hold on
% scatter(points1(:,1),points1(:,2),1)
% scatter(points2(:,1),points2(:,2),1)
% scatter(points3(:,1),points3(:,2),1)
% scatter(points4(:,1),points4(:,2),1)
% axis equal
% xlim([-2.5 2.5])
% ylim([-0.5 1.5])

opts.element_type     = 'tri6';
opts.min_angle        = 30;
opts.max_tri_area     = 0.001; %5
opts.gen_edges        =1;
tristr.points         = points;
tristr.segments       = uint32(segments);
% tristr.regions        = [96,88,1,-1;85,15,1,-1;80,-50,1,-1;0,-98,1,-1;-40,-25,1,-1;-50,20,1,-1;-25,50,1,-1;-20,80,1,-1;-90,-90,2,-1]';
tristr.holes=[rc'];
MESH = mtriangle(opts, tristr);
nodes=MESH.NODES';
elems=MESH.ELEMS';

end  