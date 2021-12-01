clear all
close all
clc

% Sections are numbered starting from the tip of the rocket onward
% For example a 4-section rocket would be
%
%
%             ^  <-- D(1) (= 0)
%   L(1)    /   \
%          /  1  \
%         --------- D(2)
%         |       |
%         |       |
%   L(2)  |   2   |
%         |       |
%         |       |
%         |       |
%         ---------    D(3)
%        /         \
% L(3)  /     3     \
%      /             \
%      ---------------   D(4)
%     |               |
%     |               |
%     |               |
% L(4)|       4       |
%     |               |
%     |               |
%       -------------     D(5)
%
%
%
% This script is modular, just add more numbers to the array L(), D()
% to add sections
%



% Sections length [m]
L(1)= 2.2;
L(2)= 4.123;
L(3)= 3.412;
L(4)= 1.12341;

% Array of diameters [m]
D(1) = 0;
D(2) = 1.451;
D(3) = 1.251;
D(4) = 1.5125;
D(5) = 3.55125;

L_tot = sum(L);

% Number of elements along the length and circumference
N_len = 43;
N_circ = 50;

% Allocate a number of element along the length proportional to
% the length of that section
q=length(L);

for i = 1:q;
  N(i) = fix(N_len*L(i)/sum(L));
end




% Write the distribution of radii R_dis
R_dis=[];
for j = 1:q-1 % j = section
  for i = 1:N(j)
    R(j,i) = D(j)/2 + (D(j+1) - D(j))/(2 * N(j))  * (i-1);
  end
 R_dis = [fliplr(R(j,1:N(j))) R_dis];
end

% Deal with last section separately to add an additional element
for i=1:N(q)+1
R(q,i) = D(q)/2 + (D(q+1) - D(q))/(2 * N(q))  * (i-1);
end

R_dis = [fliplr(R(q,1:N(q)+1)) R_dis];





%R = 0.5*[1.5*ones(1,10) 1*ones(1,30) 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0 ];

[x,y,z] = cylinder(R_dis,N_len);

z = sum(L)*z; %rescale the length

surf(x,y,z); axis equal
