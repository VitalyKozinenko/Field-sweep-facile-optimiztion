clear;
%-------------------------
%Setting the total number N of spins in the molecule
%-------------------------
N=9;
%-------------------------
%Operators for N spins 1/2
%-------------------------
[Ix,Iy,Iz,Ip,Im] = Operators(N,1/2);
%-------------------------
%Density matrix with spins 1 and 2 being pH2 protons in singlet state and other nuclei non polarized
%-------------------------
rhoP=eye(2^N)/2^N-(Iz{1}*Iz{2}+Ix{1}*Ix{2}+Iy{1}*Iy{2})/(2^(N-2));
%-------------------------
%Constructing Zeeman Hamiltonian
%-------------------------
ppm = {5.25,'1H'; 5.36,'1H';5.96,'1H'; 4.53,'1H'; 4.53,'1H'; 2.94,'1H'; 2.94,'1H';2.94,'1H'; 160, '13C'};
[Hz1] = ZeemanHam(N,ppm);
%-------------------------
%Constructing J-coupling Hamiltonian
%-------------------------
J = [0 10.5 1.4 1.7 1.7 0 0 0 0.2; 0 0 17.2 5.6 5.6 0 0 0 -0.18; 0 0 0 1.7 1.7 0 0 0 0.08; 0 0 0 0 0 0 0 0 3.1; 0 0 0 0 0 0 0 0 3.1; 0 0 0 0 0 0 0 0 1.6; 0 0 0 0 0 0 0 0 1.6; 0 0 0 0 0 0 0 0 1.6; 0 0 0 0 0 0 0 0 0];

[Hj1] = JHam(N,J,'strong');
%-------------------------
%Initiating arrays for field and magnetization values 
%-------------------------
SweepStart=10e-9; %First sweep point in T
SweepEnd=1e-6; %Last sweep point in T
B=linspace(SweepStart,SweepEnd,1000);
Net13C=zeros(length(B),1);
%-------------------------
%Calculation of the 13C PHIP field dependence
%-------------------------
for b=1:length(B)
   rho=rhoP;
%-------------------------
[vv,ee]=eig((Hz1*B(b)+Hj1)/2/pi);
%-------------------------
vv1=inv(vv);
rho0=vv*diag(diag(vv\rho*vv))/vv;
Net13C(b)=real(trace(rho0*Iz{N}/2^(N-1))/trace(Iz{N}/2^(N-1)*Iz{N}/2^(N-1)));
%-------------------------
end
%-------------------------
%Smoothing field dependence
%-------------------------
Net13C_s = smooth(Net13C,51,'sgolay',2);
%-------------------------
%Calculating time profile of the optimal sweep
%-------------------------
tau_sw=cumsum(real(Net13C_s));
tau_sw=tau_sw/tau_sw(length(B));
%-------------------------
%Calculating field profile of the optimal sweep
%-------------------------
x=tau_sw;
v=B;
xq=linspace(0,1,1000);
xq=transpose(xq);
B_opt=interp1(x,v,xq,'spline');
%-------------------------
%Calculating 13C polarization dependence on sweep duration for linear and
%optimal profiles
%-------------------------
%-------------------------
%This block is the most time consuming, if only interested in optimal
%profile - better to comment everything from here to Protting results
%-------------------------
delta_t=0.001;
k=10;
T=linspace(0,60,k); %Sweep durations in s
Netlin=zeros(length(T),1);
Netopt=zeros(length(T),1);
parfor t=1:k
    rho1lin=rhoP;
    rho1opt=rhoP;
    T_sw=T(t);
    n=T_sw/delta_t;
    Blin=linspace(SweepStart,SweepEnd,n);
    xq1=linspace(0,1,n);
    xq1=transpose(xq1);
    x=tau_sw;
    v=B;
    Bopt = interp1(x,v,xq1,'spline');
    Bopt=transpose(Bopt);

    for b=1:n
       rho1lin=expm(-(Hz1*Blin(b)+Hj1)*1i*delta_t)*rho1lin*expm((Hz1*Blin(b)+Hj1)*1i*delta_t); 
        rho1opt=expm(-(Hz1*Bopt(b)+Hj1)*1i*delta_t)*rho1opt*expm((Hz1*Bopt(b)+Hj1)*1i*delta_t); 
          %-------------------------
    end 
  Netlin(t)=real(trace(rho1lin*Iz{N}/2^(N-1))/trace(Iz{N}/2^(N-1)*Iz{N}/2^(N-1)));
  Netopt(t)=real(trace(rho1opt*Iz{N}/2^(N-1))/trace(Iz{N}/2^(N-1)*Iz{N}/2^(N-1)));
end
%-------------------------
%Plotting results
%-------------------------
subplot(4,1,1);
plot(B,Net13C,'g',B,Net13C_s);
legend('13C','13C smoothed')
subplot(4,1,2);
plot(B,tau_sw);
subplot(4,1,3);
plot(xq,B_opt,'r',xq,B);
legend('Optimal','Linear')
subplot(4,1,4);
plot(T,Netopt,'r',T,Netlin);
legend('Optimal','Linear')

