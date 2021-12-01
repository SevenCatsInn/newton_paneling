clearvars
close all
clc

% International standard atmosphere (from atomsisa matlab function)
rho = 0.4127; % Density
a = 299.4633; % Sound Speed


Vmag = a*5; %Velocity magnitude [m/s]

alpha_deg = 10 %incidence angle [deg]
alpha = alpha_deg*pi/180; %incidence angle [rad]

V = [Vmag*sin(alpha) 0 -Vmag*cos(alpha)]' % Velocity vector

q = 0.5 * rho * norm(V).^2; %Dynamic pressure


N_circ = 20; % Elements along the circumference
N_lenT = 50; % Elements along the axis (Target, not exact, see geometry_func for explanation)
L = [2 6 4 5]; % Vector of sections lengths (see geometry_func)
Diam = [0 1 1 2 2]; % Vector of section diameters (see geometry_func)

% Diam always has one more element than L

% x y z are the coordinates of the mesh midpoints
% N is the vector containing the number of elements
% for each section

[x,y,z,N] = geometry_func(L,Diam,N_circ,N_lenT);


N_len=sum(N); % Real number of elements along the axis

dim = size(x);
NORM=[]; %normals
CENT=[]; %centers
AREA=[]; %areas

% Construction of the normals and centers of the panels
for j=1:dim(2)-1

 % A to D points of every panel labeled counter clockwise
 % a and b are the midpoints,needed to contruct the centers
 %
 % D--b--C
 % |     |
 % A--a--B
 %
 %
 % TMP__ refers to a single vertical row of panels, which is then
 % assembled outside of the inner 'for' cycle to form the complete
 % matrix (temporary, it gets rewritten every j iteration)

  for i = 1:dim(1)-1
    C = [x(i+1,j+1) y(i+1,j+1) z(i+1,j+1)];
    A = [x(i,j) y(i,j) z(i,j)];
    DIAG1 = C-A;

    D = [x(i+1,j) y(i+1,j) z(i+1,j)];
    B = [x(i,j+1) y(i,j+1) z(i,j+1)];
    DIAG2 = D-B;

  % Matrix of centers of that column of panels
        a = (A+B)'/2; % midpoint of lower panel side
        b = (C+D)'/2; % midpoint of upper panel side
        TMPCENT(:,i) = (a+b)/2;

% Matrix of normals of that column of panels
        TMPNORM(:,i) = cross(DIAG1,DIAG2);
        TMPNORM(:,i) = TMPNORM(:,i)/norm(TMPNORM(:,i));
% Vector of panel areas for that column of panels
        TMPAREA(i) = 0.5*norm(cross(DIAG1,DIAG2));
  end

% Assemble matrices containing the normals and the centers
NORM = [NORM TMPNORM];
CENT = [CENT TMPCENT];
AREA = [AREA TMPAREA];

end

% Panel numbering (counter-clockwise)
%   ^     ^
%  |4|   |8|
%  ---   ---
% | 3 | | 7 |
%  ---   ---
% | 2 | | 6 |
%  ---   ---
% | 1 | | 5 |
%  ---   ---



for i=1:length(NORM(1,:))
  theta(i) = acos(dot(V,NORM(:,i)) / (norm(V)*norm(NORM(:,i))));

  % cp_prov(i) = 2 * (sin(theta(i)))^2; % Without shadow effect

  % Implement the shadow effect
  if theta (i) >= pi/2
    cp(i) = 1.8 * (cos(theta(i)))^2;
  else
    cp(i) = 0;
  end

end

%%%%%% POST-PROCESSING %%%%%%

% Reshaping the outputs of cp, area, normals
% Shaping the cp for all the panels in the form of a matrix to
% plot the colormap
CP = reshape(cp,[N_len,N_circ]);
AREA_RES = reshape(AREA,[N_len,N_circ]); % Reshaped matrix of areas
NORM_RES = reshape(NORM,[3,N_len,N_circ]);

% All the reshaped matrices have as columns the vertical columns
% of panel along the axis, the first row corresponds to the bottom
% panels, while the last row are the tip panels

% Resultant of forces calculation
for i = 1:length(cp)
  dF_tot(:,i) = -cp(i) * q * AREA(i) * NORM(:,i);
end

% Sum of the forces
F_tot = sum(dF_tot,2)


% Plot of the panels
Q = surf(x,y,z,CP); axis equal; hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
colormap jet

% Label of the colorbar
n = colorbar;
n.Label.String = 'C_p';

%set(gca,'ColorScale','log'); %Logarithmic colorbar

% set(gca,'ColorScale','log') % Logarithmic color scale
                              % use for improved resolution

% Plot velocity vector (see help for quiver3)
resc = 5e2; %rescaling factor
V=V/resc; %rescale velocity magnitude just for the plot
quiver3([-2],[0],[12],[V(1)],[V(2)],[V(3)],'linewidth',2);
V=V*resc; %scale back to original

% Plot force vector
resc = alpha_deg*1e4; % rescaling factor,
                      % we use alpha_deg because it works well

F_tot=F_tot/resc; %rescale force magnitude just for the plot
quiver3([0],[0],[5],[F_tot(1)],[F_tot(2)],[F_tot(3)],'linewidth',2);
F_tot=F_tot*resc; %scale back to original

% Plot normals to faces
%quiver3(CENT(1,:),CENT(2,:),CENT(3,:),NORM(1,:),NORM(2,:),NORM(3,:));
hold on;

legend('','Velocity','Force','Normals')

versV = V / norm(V);
DRAG = dot(versV,F_tot);


% Forces on sections calculation

prevN = 0; % Number of panel layers underneath the one we're considering

for p=1:length(N) % p = section number
  for j=1:N_circ % j = along the circumference
    for i = 1:N(p) % i = along the length
      dF(:,i+prevN,j) = -CP(i+prevN,j) * q * AREA_RES(i+prevN,j) * NORM_RES(:,i+prevN,j);
    end
  end
  prevN = prevN + N(p);

% Successive buildup
F(:,p) = sum(dF,[2 3]);

end

% Extract the forces on the components
for b = 1:3
  F(:,b+1) = F(:,b+1) - sum(F(:,1:b),2);
end

% Flip the vector for consistency with the numbering used so far

% F = Vector of forces on the components
F = fliplr(F)
