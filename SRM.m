function [s0,s1,s2]=SRM(img,m0,m1,m2,lmd0,lmd1,lmd2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% img: DoFP image.                                                        %
% m0, m1, and m2: Modulation parameters                                   %
% lmd0, lmd1, and lmd2: Regularization parameters.                        %
% s0, s1, and s2: Stokes parameters.                                      %
% Using pcg to solve A*s = b.                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y] = size(img);

[y,x] = meshgrid(0:Y-1,0:X-1);
k1=mean(mean((m1+m2).*cos(pi*x)))^2;
k2=mean(mean((m1-m2).*cos(pi*y)))^2;
lmd1 = lmd1*k1;
lmd2 = lmd2*k2;

D = [sparse(1:X*Y,1:X*Y,m0(:)),sparse(1:X*Y,1:X*Y,m1(:)),sparse(1:X*Y,1:X*Y,m2(:))];
Hyy = cv2mtx([1 -2 1],X,Y);
Hxx = cv2mtx([1 -2 1]',X,Y);
Hxy = cv2mtx([1 -1;-1 1],X,Y);
H = [Hxx*sqrt(X/(X-2));Hyy*sqrt(Y/(Y-2));Hxy*sqrt(X/(X-1)*Y/(Y-1)*2)];
R = kron([lmd0 0 0;0 (lmd1+lmd2)/4 (lmd1-lmd2)/4;0 (lmd1-lmd2)/4 (lmd1+lmd2)/4], H'*H);

A = D'*D+R;
b = D'*img(:);
tol = 1e-14;
maxit = 500;

try
    L = ichol(A);
    [s,flag1] = pcg(A,b,tol,maxit,L,L');
    if flag1==1
        disp('Flag1: PCG iterated maxit times but did not converge.');
    end
catch
    try
        L = ichol(A,struct('michol','on'));
        [s,flag2] = pcg(A,b,tol,maxit,L,L');
        if flag2==1
            disp('Flag2: PCG iterated maxit times but did not converge.');
        end
    catch
        [s,flag3] = pcg(A,b,tol,maxit);
        if flag3==1
            disp('Flag3: PCG iterated maxit times but did not converge.');
        end
    end
end

s0 = reshape(s(1+X*Y*0:X*Y*1),X,Y);
s1 = reshape(s(1+X*Y*1:X*Y*2),X,Y);
s2 = reshape(s(1+X*Y*2:X*Y*3),X,Y);

end

function T = cv2mtx(H,M,N)

[P, Q] = size(H);

if (nargin == 2)
    if (numel(M) ~= 2) || (min(size(M)) ~= 1)
        error(message('images:convmtx2:invalidInput'))
    else
        N = M(2);
        M = M(1);
    end
end

inputNames = {'H' 'M' 'N'};
for p=1:1:3
    validateattributes(H,{'double'},{},mfilename,inputNames{p},p);
end

Tf = convmtx2(H,M,N);
indM = zeros(P+M-1,Q+N-1);
indM(:) = 1:(P+M-1)*(Q+N-1);
ind = indM(size(H,1):end-size(H,1)+1,size(H,2):end-size(H,2)+1);
T = Tf(ind,:);

end
