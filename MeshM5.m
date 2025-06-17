function [nodes,elems]=MeshM5(rc,a,b,n_ell,alpha,n_crack,box_a)

theta=linspace(0,2*pi,n_ell+1);
theta(end)=[];

cont=[-box_a, -box_a; box_a, -box_a; box_a, box_a; -box_a, box_a]';

points1=zeros(n_crack*n_ell,2);
for i=1:n_crack
    points1((i-1)*n_ell+1:i*n_ell,:)=[a(i)*cos(theta'),b(i)*sin(theta')];
    points1((i-1)*n_ell+1:i*n_ell,:)=points1((i-1)*n_ell+1:i*n_ell,:) ...
        *[cos(alpha(i)) -sin(alpha(i)); sin(alpha(i)) cos(alpha(i))];

end

points1=points1+repelem(rc,n_ell,1);

points=[cont, points1'];
segments_box=[1,2,3,4; 2,3,4,1];


segments_ellipse=4+[1:n_ell;2:n_ell,1];
seg_ELL=repmat(segments_ellipse,1,n_crack-1)+n_ell*repmat(repelem(1:n_crack-1,n_ell),2,1);
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
opts.max_tri_area     = 1; %5
opts.gen_edges        =1;
tristr.points         = points;
tristr.segments       = uint32(segments);
% tristr.regions        = [96,88,1,-1;85,15,1,-1;80,-50,1,-1;0,-98,1,-1;-40,-25,1,-1;-50,20,1,-1;-25,50,1,-1;-20,80,1,-1;-90,-90,2,-1]';
tristr.holes=[rc'];
MESH = mtriangle(opts, tristr);
nodes=MESH.NODES';
elems=MESH.ELEMS';

end  