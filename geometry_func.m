function [x,y,z,N] = geom(L,Diam,N_circ,N_lenT)

% Sections are numbered starting from the tip of the rocket onward
% For example a 4-section rocket would be
%
%
%             ^  <-- Diam(1) (= 0)
%   L(1)    /   \
%          /  1  \
%         --------- Diam(2)
%         |       |
%         |       |
%   L(2)  |   2   |
%         |       |
%         |       |
%         |       |
%         ---------    Diam(3)
%        /         \
% L(3)  /     3     \
%      /             \
%      ---------------   Diam(4)
%     |               |
%     |               |
%     |               |
% L(4)|       4       |
%     |               |
%     |               |
%       -------------     Diam(5)
%
%
%
% This script is modular, just add more numbers to the array L(), Diam()
% to add sections
%
%
%
% Note: N_lenT is the target and N_lenR is the actual number
% of elements along the axis, they might slightly differ as
% the number of elements for each section must be a positive integer
% and therefore some rounding is involved



L_tot = sum(L);


% Allocate a number of element along the length proportional to
% the length of that section
q=length(L);


for i = 1:q;
  N(i) = fix(N_lenT*L(i)/sum(L));
end

N_lenR = sum(N);


% Write the distribution of radii R_dis
R_dis=[];
for j = 1:q-1 % j = section
  for i = 1:N(j)
    R(j,i) = Diam(j)/2 + (Diam(j+1) - Diam(j))/(2 * N(j))  * (i-1);
  end
 R_dis = [fliplr(R(j,1:N(j))) R_dis];
end

% Deal with last section separately to add an additional element
for i=1:N(q)+1
R(q,i) = Diam(q)/2 + (Diam(q+1) - Diam(q))/(2 * N(q))  * (i-1);
end

R_dis = [fliplr(R(q,1:N(q)+1)) R_dis];

[x,y,z] = cylinder(R_dis,N_circ);

z = sum(L)*z; %rescale the length
