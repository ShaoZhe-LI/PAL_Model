% =========================================================================
% INTRODUCTION
% - Implement the miu-order Hankel transformation on the columns of the 
% input matrix h
% -------------------------------------------------------------------------
% INPUT
% h - target matrix (matrix before transformation)
% N - number of transformation points (the rows of the target matrix)
% Nz - number of transformations (the columns of the target matrix)
% Nh - truncation length of h
% NH - truncation length of H
% alpha - exponential coefficients of geometric series
% x0 - first element of x1
% x1 - normalized sampling points' coordinates in \rho direction
% k0 - interpolation coefficient
% miu - order of transformation
% OUTPUT
% H - matrix after transformation
% =========================================================================

% h=vs1;N=N_FHT;Nz=1;Nh=Nh_v;NH=NH_v;alpha=a_solve;miu=m;
function H=m_FHT(h,N,Nz,Nh,NH,alpha,x0,x1,k0,miu)
    x1 = x1(:).';   % x1 as row vector
    assert(numel(x1) == N, ...
    'm_FHT: x1 must have length N (%d), but got %d.', N, numel(x1));
    assert(size(h,1) == N, ...
    'm_FHT: h must have N rows (%d).', N);
    Nf=Nh*NH;
    phi=zeros(2*N,Nz);
    dh1=(exp(alpha*([0:N-2]+1-N).')./((x1(1:end-1)).')).^miu;
    dh2=(exp(alpha*([0:N-2]+1-N).')./((x1(2:end)).')).^miu;
    h1=h(1:end-1,:).*exp(alpha*([0:N-2]+1-N).').*dh1;
    h2=h(2:end,:).*exp(alpha*([0:N-2]+1-N).').*dh2;
    phi(1:N-1,:)=(h1-h2);
    phi(1,:)=k0*phi(1,:);
    phi(N,:)=h(end,:)./((x1(end))^miu);Phi=fft(phi);
    jnn=Nf*x0*exp(alpha*([0:2*N-1]+1-N));jn0=(besselj(miu+1,jnn)).';
    J=ifft(jn0*ones(1,Nz));
    G=fft(Phi.*J);
    H=G(1:N,:)./(Nf*x1.')*Nh^2;
end