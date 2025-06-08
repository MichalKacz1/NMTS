a = 0;
b = 3;
A_ = 0.02;
% Przypadek a
z0_a = -abs(b - a);
G_a = tf([1, -z0_a], conv([1 a], [1 b]));

z0_b = 0.5 * (a + b);
G_b = tf([1, -z0_b], conv([1 a], [1 b]));

chosen_G = 2; %Choose transmitance a - 1 b - 2
if chosen_G==1
    G=G_a;
else
    G=G_b;
end
des_GM=4;
des_PM=deg2rad(40);
omega3=1000;
Ms=1.01:0.01:4;
GM=Ms./(Ms-1);
PM=2*asin(0.5./Ms);
subplot(2,1,1)
plot(Ms,mag2db(GM),Ms,zeros(1,length(Ms))+des_GM)
subplot(2,1,2)
plot(Ms,PM,Ms,des_PM.*ones(1,length(Ms)))

%select omega
if chosen_G==1
    omega3_a=1;
else
    omega3_b=0.5*z0_b;
end

Ms=1.45;
Mt=Ms+1;
%weight functions for cost calculation
wt=tf((Mt*s+omega3_a)/(Mt*omega3_a));
ws=tf((s/Ms+omega3_a)/(s+A_*omega3_a));

% P1 = augw(G_a,ws,wt,[]);


% GCP
P=[ws -ws*G;0 wt*G;1 -G];
% state space
k=1/Ms;
A=[0 1 0;0 0 0; 0 0 0];
A(2,1)=-prod(pole(G));
A(2,2)=sum(pole(G));
A(3,1)=zero(G);
A(3,2)=-1;
A(3,3)=pole(ws);
B1=[0;0;1];
B2=[0;1;0];
C1=[k*zero(G) -k k*pole(ws)-zero(ws); -zero(G)-prod(pole(G)) sum(pole(G)-zero(G)+zero(wt)) 0];
C2=[-zero(G) 1 0];
D11=[0; k];
D12=[0; 1];
D21=[1];
D22=[0];

%Checking conditions
r=rank(ctrb(A, B2));
if r<min(size(A))
    O = ctrb(A, B2);
    N = null(O);
    A_uncontr = N' * A * N;
    disp(A_uncontr)
end


r=rank(obsv(A, C2));
if r<min(size(A))
    O = obsv(A, C2);
    N = null(O);
    A_unobs = N' * A * N;
    disp(A_unobs)
end

syms w real

n = size(A, 1);
jwI = 1i*w*eye(n);   % jωI
M = [jwI-A -B2;
     C1 D12];

r = rank(M);         % Symbolic rank
disp('Symbolic rank of the first matrix is:');
disp(r);
disp('And size of it is:');
disp(size(M));

M = [jwI-A -B1;
     C2 D21];

r = rank(M);         % Symbolic rank
disp('Symbolic rank of the second matrix is:');
disp(r);
disp('And size of it is:');
disp(size(M));

% find gamma
gamma=0;
while true
    gamma=gamma+0.01;
    A_wave=A-B2*transpose(D12)*C1;
    C_wave=sqrtm(transpose(transpose(C1)*(eye(size(D12,1))-D12*transpose(D12))*C1));
    X=icare(A_wave, [B2 B1],C_wave*transpose(C_wave),[eye(size(B2,2)) zeros(size(B2,2),size(B1,2)); zeros(size(B1,2),size(B2,2)),-gamma^(-2)*eye(size(B1,2))]);
    if isempty(X)
        continue
    end
    if max(eig(A).*(real(A_wave)-(B2*transpose(B2)-gamma^(-2)*(B1*transpose(B1))*X)))>=0
        continue
    end
    A_roof=A-B1*transpose(D21)*C2;
    B_roof=sqrtm(transpose(B1)*(eye(size(D21,2))-transpose(D21)*D21)*B1);
    Y=icare(A_roof, transpose([C2; C1]),B_roof*transpose(B_roof),[eye(size(C2,1)) zeros(size(C2,1),size(C1,1)); zeros(size(C1,1),size(C2,1)),-gamma^(-2)*eye(size(C1,1))]);
    if isempty(Y)
        continue
    end
    if max(eig(A).*(real(A_roof)-(transpose(C2)*C2-gamma^(-2)*(transpose(C1)*C1)*Y)))>=0
        continue
    end
    if max(abs(eig(X*Y)))>=gamma^2
        continue
    end
    break
end

%finally do the H inf
