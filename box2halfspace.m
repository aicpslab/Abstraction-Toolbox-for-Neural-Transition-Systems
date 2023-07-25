function halfspace=box2halfspace(box)
%% transfer a box into halfspace
dim= size(box,1);
% Define normal vectors
ub = box(:,2);
lb =-box(:,1);
G = zeros(2*dim,2);
g = zeros(2*dim,1);

% Coefficient matrix
G(1:dim,1:dim)=eye(dim);
G(dim+1:end,1:dim)=-1*eye(dim);

% Define scalar values
g(1:dim,1) = ub;
g(dim+1:end,1)= lb;

% Form halfspace 
halfspace=HalfSpace(G, g);

end