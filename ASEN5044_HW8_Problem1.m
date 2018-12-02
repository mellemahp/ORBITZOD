clc
clear all
close all

Omega=0.045;
delta_t=0.5;

time=linspace(0,300,300);
A=[0 1 0 0;0 0 0 -Omega;0 0 0 1;0 Omega 0 0];

Omega_a=0.045;
Omega_b=-0.045;


Fa=[1 sin(Omega_a*delta_t)/Omega_a 0 -(1-cos(Omega_a*delta_t))/Omega_a;
    0 cos(Omega_a*delta_t) 0 -sin(Omega_a*delta_t);
    0  (1-cos(Omega_a*delta_t))/Omega_a 1 sin(Omega_a*delta_t)/Omega_a;
    0 sin(Omega_a*delta_t) 0 cos(Omega_a*delta_t)]

Fb=[1 sin(Omega_b*delta_t)/Omega_b 0 -(1-cos(Omega_b*delta_t))/Omega_b;
    0 cos(Omega_b*delta_t) 0 -sin(Omega_b*delta_t);
    0  (1-cos(Omega_b*delta_t))/Omega_b 1 sin(Omega_b*delta_t)/Omega_b;
    0 sin(Omega_b*delta_t) 0 cos(Omega_b*delta_t)]

GamA=[0 0;1 0;0 0;0 1];
GamB=[0 0;1 0;0 0;0 1];

qw=10;
W=qw*[2 0.05;0.05 0.5];

Za=delta_t*[-A GamA*W*GamA';zeros(4,4) A'];
Zb=delta_t*[-A GamB*W*GamB';zeros(4,4) A'];

e_za=expm(Za);
e_zb=expm(Zb);

FinvQA=e_za(1:4,5:8);
Qa=(Fa')'*FinvQA

FinvQB=e_zb(1:4,5:8);
Qb=(Fb')'*FinvQB


%% PART B - i

rng(100)


Ha=[1 0 0 0;0 0 1 0];
Ra=[20 0.05;0.05 20];

load('hw8problem1_data')
xa=xasingle_truth;



R_block=[];
H_block=[];



for i=1:length(xasingle_truth)
    noise=mvnrnd([0 0],Ra);
    y_k(:,i)=Ha*xasingle_truth(:,i)+noise';

end

figure(1)
time=linspace(0,20,20);
plot(time,y_k(:,1:20))

xlabel('Time (s)')
ylabel('y_a(k)')
grid on


mu_a=[0 85*cos(pi/4) 0 -85*sin(pi/4)]';
P0_a = 900*diag([10 2 10 2]);

xm=mu_a;
Pm=P0_a;

for i=1:length(xasingle_truth)
    xm(:,i+1)=Fa*xm(:,i)
    Pm(i+1)=Fa*Pm(i)*Fa'+Qa;
    K(i+1)=Pm(i+1)*Ha'*inv(Ha*Pm(i+1)*Ha'+Ra);
  
    xp(i+1)=xm(i+1)+K(i+1)*(y_k(:,i+1)-Ha*xm(i+1));
    Pp(i+1)=(eye(2)-K(i+1)*Ha)*Pm(i+1);
end


