function q=spm_ss_fdr(p,dim,nc)
% FDR False discovery rate
% Q=FDR(P); returns vector Q of estimated false discovery rates (set-level q-values) from 
% a vector P of multiple-test false positive levels (uncorrected p-values)
% Q=FDR(P,dim); where P is a matrix computes the fdr along the dimension dim of P
%

if nargin<2, 
    if sum(size(p)>1)==1,dim=find(size(p)>1);
    else dim=1; end
end
if nargin<3||isempty(nc),nc=0;end

nd=length(size(p)); 
if dim~=1, p=permute(p,[dim,1:dim-1,dim+1:nd]); end

sp=size(p);
q=nan+ones(sp);
N1=sp(1);
N2=prod(sp(2:end));
for n2=1:N2,
    [sp,idx]=sort(p(:,n2));
    N1=sum(~isnan(p(:,n2)));
    if N1>0,
        if nc,  qt=min(1,N1*sum(1./(1:N1))*sp(1:N1)./(1:N1)'); 
        else    qt=min(1,N1*sp(1:N1)./(1:N1)'); end
        min1=nan;
        for n=N1:-1:1,
            min1=min(min1,qt(n));
            q(idx(n),n2)=min1;
        end
    end
end
if dim~=1, q=ipermute(q,[dim,1:dim-1,dim+1:nd]); end



