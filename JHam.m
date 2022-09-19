function [Hj] = ZeemanHam(N,J,regime)
%This function is made for defining J-coupling Hamiltonian
%
%J must have the following form: J=[J12 J13 J14 .. J1N; 0 0 .. J23 .. J2N; ..; 0 0 .. 0 0]
%
%-------------------------
%Operators for N spin 1/2
[Ix,Iy,Iz,Ip,Im] = Operators(N,1/2);

Hj=zeros(2^N);


for l=1:N-1
    for m=l+1:N
        if regime == "weak" %Weak coupling for strong fields
             Hj=Hj+J(l,m)*Iz{l}*Iz{m};
        elseif regime == "strong" %Weak coupling for arbitrary fields
             Hj=Hj+J(l,m)*(Iz{l}*Iz{m}+Ix{l}*Ix{m}+Iy{l}*Iy{m});
        end
    end     
end

Hj=Hj*2*pi;

end
