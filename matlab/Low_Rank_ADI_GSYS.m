clear all;
clc;

%load('beam.mat')
%load('beamm.mat')
%load('build.mat')
load('CDPlayer.mat')
%load('EandBn.mat')
% load('iss.mat')
%load('issm.mat')
%load('mod.mat')
%load('eady.mat')

E=eye(size(A,1));
maxiter1 = 150;
maxiter2 = 100;
restol = 10^-8; 
l = 10;
% 
% [Zc,res1]=adaptiveshift_lradi(E,A,B,l,maxiter1,restol);
% z = C * Zc;
% zz = z'*z
[Zc,res1]=adaptiveshift_lradi(E,A,B,l,maxiter1,restol);
[Zo,res1]=adaptiveshift_lradi(E',A',C',l,maxiter1,restol);
% sys = ss(A, B, C, 0);
% norm(sys)
% tic;norm(sys);toc
% 
z = C * Zc;
z1 = B' * Zo;
H2c = sqrt(trace(z'*z))
H2o = sqrt(trace(z1*z1'))
% [Zo,res2]=adaptiveshift_lradi(E',A',C',l,maxiter2,restol);

% sys = ss(A, B, C, 0);
% norm(sys)
% tic;norm(sys);toc
function [Zc,res1]=adaptiveshift_lradi(E,A,B,l,maxiter1,restol)
i=1;
Zc=[];
ip=0;
res1=zeros(1,maxiter1+2);
 m3=size(B,2);
%l=length(p);
W=B;
bnorm=norm(W'*W,'fro');         % new method to compute ADI res
%% Initial shift computation
  Bn=sprand(size(A,1),100,.8);
 [Q,~]=qr(full(Bn),0);
   p = adaptive_shift(E,A,Q,l); % compute l ADI shift parameters
   l=length(p);
while i<=maxiter1
 if ip<l
        ip=ip+1;
  else
   m2=size(Zc,2);
    Zcn=Zc(:,m2+1-l*m3:end);
   % size(Zcn,2)
    [Q,~]=qr(full(Zcn),0);
    p = adaptive_shift(E,A,Q,l);
    l=length(p);
    ip =1;
 end %
    pc=p(ip);
    V=(A+pc*E)\W;
  if isreal(pc)
          Zc=[Zc sqrt(-2*real(pc))*real(V)];
      W = W-2*pc*E*V;
  else
         beta=real(pc)/imag(pc);
         gam=2*sqrt(-real(pc));
        Zc=[Zc,...
             gam*(real(V)+beta*imag(V)),...
             gam*sqrt(beta^2+1)*imag(V)];
         V=real(V)+beta*imag(V);
         pc=conj(pc);
         W = W-4*real(pc)*E*V;
         i=i+1;
         res1(i)=norm(W'*W,'fro')/bnorm;
   fprintf(1,'step: %4d  normalized residual: %d\n',i,res1(i));
   if res1(i)<restol % ||(norm(T{i},2)/normest(Zb)<n*eps);
   break;
   end
   end
   i=i+1;
end
  
return
end

function p = adaptive_shift(E,A,V,l)
   An=V'*A*V;
   En=V'*E*V;
   rw=eig(full(An),full(En));
   rw0 = rw;
   rw = [];
   for j = 1:length(rw0)
    if real(rw0(j))<0
      rw = [rw;rw0(j)];
    end
   end
    p=lp_mnmx(rw,l);
    
end

function p = lp_mnmx(rw,l0)
if length(rw)<l0
  error('length(rw) must be at least l0.');
end

max_rr = +Inf;                      % Choose initial parameter (pair)

for i = 1:length(rw)
  max_r = lp_s(rw(i),rw); 
  if max_r < max_rr
    p0 = rw(i);
    max_rr = max_r;
  end
end  

if imag(p0)
  p = [ p0; conj(p0) ];
else
  p = p0;                            
end  

[max_r,i] = lp_s(p,rw);         % Choose further parameters.

while size(p,1) < l0 
   
  p0 = rw(i);
  if imag(p0)
    p = [ p; p0; conj(p0) ];
  else
    p = [ p; p0]; 
  end
  
  [max_r,i] = lp_s(p,rw);
    
end
end

function [max_r,ind] = lp_s(p,set)
max_r = -1;
ind = 0;
  
for i = 1:length(set)
  
  x = set(i);
  
  rr = 1;
  for j = 1:length(p)

    rr = rr*abs(p(j)-x)/abs(p(j)+x);
    
  end  
    
  if rr > max_r
    
    max_r = rr;
    ind = i;
   
  end
  
end  
end