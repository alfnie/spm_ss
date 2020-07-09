function [z,iC,a,b,w,k]=spm_ss_fls(f,x,Niter)
% SPM_SS_FLS maximum likelihood estimation of between-subjects covariance
%

% [nill,iC]=spm_ss_fls({x,n});
%  Maximum Likelihood covariance estimation. Returns whitening vector iC.
%  the data x (mx1) is assumed to come from a distribution with covariance x_i = a*(b+1/n_i) 
%  Maximum likelihood estimation of these parameters is obtained by
%  reducing the ML search to a line-search over the parameter b. The
%  default range of b where the optimal is search for is 10^-4 to 10^4
%  (range = [-4,4] in log10 units). This behavior can be changed using
%  spm_ss_fls({x,n},range). The default number of iterations is 20. This
%  behavior can be changed using spm_ss_fls({x,n},range,Niter).
%

% Fibonacci Line-Search algorithm
% x_final = spm_ss_fls(f,x_initial,Niter)
%  where f(x) is a function to be maximized over the initial range
%  x_initial (assuming uni-modality within this range). x_final is the
%  range where the optimal is found to be after Niter iterations
%

if nargin<3||isempty(Niter), Niter=20; end  %number of iterations (Niter=20 => final parameter range < initial parameter range / 1e4)
if nargin<2||isempty(x), x=[-4,4]; end      %parameter search initial range

if iscell(f),
    N=f{2};
    spm_ss_LL(f);
    f=@spm_ss_LL;
end

n=[1,1,zeros(1,Niter-2)];for niter=3:Niter,n(niter)=n(niter-2)+n(niter-1);end; k=n(1:end-1)./n(2:end);
%K=(1+sqrt(5)/2);K=K/(K+1);k=K+zeros(size(1,Niter)); golden-section-search

z=x;
dx=k(end)*(x(2)-x(1));
x=[x(2)-dx,x(1)+dx];
y=[f(x(1)),f(x(2))];
for niter=1:Niter-2,
    %plot(10.^x,y,'ko'); hold on; 
    dx=k(end-niter)*dx;
    if y(1)>y(2),   %left
        z(2)=x(2);              %parameter search current range
        x=[x(2)-dx,x(1)];       %new points to evaluate
        y=[f(x(1)),y(1)];
    else            %right
        z(1)=x(1);
        x=[x(2),x(1)+dx];
        y=[y(2),f(x(2))];
    end
end

if nargout>1,
    [nill,a,b,w]=spm_ss_LL(mean(z));
    iC=zeros(size(N));iC(N>1e-10)=sqrt(w/max(eps,a));
end

end

% covariance model C_ii = a*(b+1/N_i)
function [L,a,b,w,k]=spm_ss_LL(z)
persistent iN x;
if isnumeric(z),
    b=10.^z;
    w=1./(b+iN);
    k=(w'*x)/max(eps,sum(w));
    a=0;for n1=1:size(x,2),a=a+mean(w.*(x(:,n1)-k(n1)).^2);end
    L=mean(log(max(eps,w)))-log(max(eps,a))-(z/2).^2;
else
    idx=find(z{2}>1e-10);
    x=z{1}(idx,:);
    iN=1./max(eps,z{2}(idx));
end
end
