function [Hz] = ZeemanHam(N,ppm)
%This function is made for defining Zeeman Hamiltonian
%-------------------------
%Operators for N spin 1/2
[Ix,Iy,Iz,Ip,Im] = Operators(N,1/2);

Hz=zeros(2^N);

for i=1:N
    
   switch ppm{i,2}
   case '1H'
      Hz=Hz+(-267)*(1+ppm{i,1}*(1e-6))*Iz{i};
   case '13C'
      Hz=Hz+(-67)*(1+ppm{i,1}*(1e-6))*Iz{i};
   case '15N'
      Hz=Hz+(40)*(1+ppm{i,1}*(1e-6))*Iz{i};
    ...
   otherwise
      Hz=0;
   end

end

Hz=Hz*(1e+6);

%Gyromagnetic ratio is in 10^6 rad /s /T



end
