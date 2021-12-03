clear all
close all
clc

load("vars_for_struct_analysis.mat")
load("aerodynamic_forces_entire.mat");


%to be derived from trajectory results
g0=9.81;
R_e = 6378.1e3;
r = R_e + alt;
g = g0*(R_e/r)^2;
gamma=89.5*pi/180;


mass_flare=50.2968;
mass_motor=3528.062;
mass_interstage=18.91;
mass_secondstage_body=472;
mass_nose=1.5; 
m0=4000;
m_propellant= 3153.83;
mass_consumed= m0-mass_whole_rocket(index);

L_flare=0.7197;
L_interstage=0.7007;
L_bodymotor=4.064;
L_upper=1.6;
L_nose=0.9;

L_tot=L_flare+L_interstage+L_bodymotor+L_upper+L_nose;

Dbodymotor=0.8255;
Dupperstage=0.45;
Dexit=1.2604;

%creation of a vector with the masses
massdistr=[mass_flare (mass_motor-mass_consumed)/2 (mass_motor-mass_consumed)/2 mass_interstage mass_secondstage_body/2 mass_secondstage_body/2 mass_nose ];

%position of the center of mass
xcg_flare = L_flare/4*(((Dexit/2)^2+2*Dexit/2*Dbodymotor/2+3*(Dbodymotor/2)^2)/(4*((Dexit/2)^2+Dexit/2*Dbodymotor/2+(Dbodymotor/2)^2)));
xcg_motor1 = L_flare+ L_bodymotor/4;
xcg_motor2 = L_flare + 3*L_bodymotor/4;
xcg_interstage = L_flare + L_bodymotor + L_interstage/4*(((Dbodymotor/2)^2+2*Dbodymotor/2*Dupperstage/2+3*(Dupperstage/2)^2)/(4*((Dbodymotor/2)^2+Dbodymotor/2*Dupperstage/2+(Dupperstage/2)^2)));
xcg_upperstage1= L_flare + L_bodymotor + L_interstage + L_upper/4;
xcg_upperstage2= L_flare + L_bodymotor + L_interstage + L_upper*3/4;
xcg_nose= L_flare + L_bodymotor + L_interstage + L_upper + L_nose/4;



%create a vector with the center of mass
x_discretized_cg= [xcg_flare xcg_motor1 xcg_motor2 xcg_interstage xcg_upperstage1 xcg_upperstage2 xcg_nose];
%compute moment of inertia
%global center of mass
tool=x_discretized_cg.*massdistr;
x_cg=sum(tool)/sum(massdistr);

%x_cg=(massdistr(1)*x_discretized_cg(1)+massdistr(2)*x_discretized_cg(2)+massdistr(3)*x_discretized_cg(3)+massdistr(4)*x_discretized_cg(4)+massdistr(5)*x_discretized_cg(5)+massdistr(6)*x_discretized_cg(6)+massdistr(7)*x_discretized_cg(7))/sum(massdistr);
%define a vector with the distances from the center of mass  x-x_cg
distances=x_cg-x_discretized_cg; %forces on the base will produce a "positive pitching moment"
%moment of inertia around the center of mass
tool2=distances.^2.*massdistr;
I=sum(tool2);
%I=(massdistr(1)*distances(1)^2+massdistr(2)*distances(2)^2+massdistr(3)*distances(3)^2+massdistr(4)*distances(4)^2+massdistr(5)*distances(5)^2+massdistr(6)*distances(6)^2+massdistr(7)*distances(7)^2)/sum(massdistr);

%Normal_force=[F(1,1) F(1,2)/2 F(1,2)/2 F(1,3) F(1,4)/2 F(1,4)/2 F(1,5)]*correctionfactor_normal;
%Drag=-[F(3,1) F(3,2)/2 F(3,2)/2 F(3,3) F(3,4)/2 F(3,4)/2 F(3,5)]*correctionfactor_drag;

%try
Normal_force=ones(1,7)*Lift_whole_rocket(index)/7;
Drag=ones(1,7)*Drag_whole_rocket(index)/7;
T=T_whole_rocket(index);

axial_acceleration=(T-sum(Drag))/sum(massdistr)-g*sin(gamma);
normal_acceleration=sum(Normal_force)/sum(massdistr);
ext_moment=sum(distances.*Normal_force);
angular_acceleration=ext_moment/I;


axial_refined(:,1)= - sum(Drag(1:1)) - sum(massdistr(1:1)).*axial_acceleration-g*sum(massdistr(1:1))*sin(gamma);
axial_refined_matrix=[ axial_refined(:,1)];
for i=2:7
axial_refined(:,i)= +T - sum(Drag(1:i)) - sum(massdistr(1:i)).*(axial_acceleration)-g*sum(massdistr(1:i))*sin(gamma);
axial_refined_matrix=[axial_refined_matrix axial_refined(:,i)];
end

axial_load_abs=abs(axial_refined_matrix);
max_ax_load_refined=max(axial_load_abs,[],'all');

%assume the yield stress of aluminum alloy 7075 T6 of 450 Mpa
%assume a derating factor of 1.25

DF=0.75;
sf=1.5;
sigma_yield=450e6;
sigma_material= DF*sigma_yield;
%thickness_ax=max(axial_load_abs,[],'all')/(2*pi*D/2*sigma_material);

%if a safety factor of 1.5 is included

% thickness_ax=sf*max_ax_load/(2*pi*D/2*sigma_material);
thickness_ax=sf*max_ax_load_refined/(2*pi*Dbodymotor/2*sigma_material);

% shear and bending
% shear_refined_matrix=[];
%  for k=1:7 
%  shear_refined(:,k)=sum(Normal_force(1:k))-sum(massdistr(1:k))*normal_acceleration;
%  shear_refined_matrix=[shear_refined_matrix shear_refined(:,k)];
%  end

shear_refined(:,1)=Normal_force(1)-massdistr(1)*normal_acceleration- massdistr(1)*distances(1)*angular_acceleration;
shear_refined(:,2)=shear_refined(:,1)+Normal_force(2) -massdistr(2)*normal_acceleration -massdistr(2)*distances(2)*angular_acceleration;
shear_refined(:,3)=shear_refined(:,2)+Normal_force(3) -massdistr(3)*normal_acceleration -massdistr(3)*distances(3)*angular_acceleration;
shear_refined(:,4)=shear_refined(:,3)+Normal_force(4) -massdistr(4)*normal_acceleration -massdistr(4)*distances(4)*angular_acceleration;
shear_refined(:,5)=shear_refined(:,4)+Normal_force(5) -massdistr(5)*normal_acceleration -massdistr(5)*distances(5)*angular_acceleration;
shear_refined(:,6)=shear_refined(:,5)+Normal_force(6) -massdistr(6)*normal_acceleration -massdistr(6)*distances(6)*angular_acceleration;
shear_refined(:,7)=shear_refined(:,6)+Normal_force(7) -massdistr(7)*normal_acceleration -massdistr(7)*distances(7)*angular_acceleration;


bending_refined(:,1)= Normal_force(:,1)-Normal_force(:,1);
bending_refined(:,2)=shear_refined(:,1)*(x_discretized_cg(2)-x_discretized_cg(1));
bending_refined(:,3)=bending_refined(:,2)+shear_refined(:,2)*(x_discretized_cg(3)-x_discretized_cg(2));
bending_refined(:,4)=bending_refined(:,3)+shear_refined(:,3)*(x_discretized_cg(4)-x_discretized_cg(3));
bending_refined(:,5)=bending_refined(:,4)+shear_refined(:,4)*(x_discretized_cg(5)-x_discretized_cg(4));
bending_refined(:,6)=bending_refined(:,5)+shear_refined(:,5)*(x_discretized_cg(6)-x_discretized_cg(5));
bending_refined(:,7)=bending_refined(:,6)+shear_refined(:,6)*(x_discretized_cg(7)-x_discretized_cg(6));

 
bending_load_abs=abs(bending_refined);
max_bend_load=max(bending_load_abs,[],'all');
thickness_bending=sf*max_bend_load/(pi*Dbodymotor^2/4*sigma_material);

total_thickness =thickness_bending+thickness_ax;



%compute the axial load on the base when on ground
axial_static=sum(massdistr)*g0;
thickness_ax_static=sf*axial_static/(2*pi*Dbodymotor/2*sigma_material);

