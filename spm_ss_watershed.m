function [D,P] = spm_ss_watershed(A,IDX)
% SPM_SS_WATERSHED watershed segmentation
%
% C=spm_ss_watershed(A);
% C=spm_ss_watershed(A,idx);
%

% note: assumes continuous volume data (this implementation does not work well with discrete data). In practice this means having sufficiently-smoothed volume data
%


sA=size(A);

%zero-pad&sort
if nargin<2, IDX=find(~isnan(A)); IDX=IDX(:); else IDX=IDX(:); end
[a,idx]=sort(A(IDX)); idx=IDX(idx); 
[pidx{1:numel(sA)}]=ind2sub(sA,idx(:));
pidx=mat2cell(1+cat(2,pidx{:}),numel(pidx{1}),ones(1,numel(sA)));
eidx=sub2ind(sA+2,pidx{:});
sA=sA+2;
N=numel(eidx);

%neighbours (max-connected; i.e. 26-connected for 3d)
[dd{1:numel(sA)}]=ndgrid(1:3);
d=sub2ind(sA,dd{:});
d=d-d((numel(d)+1)/2);d(~d)=[];

%assigns labels
C=zeros(sA);P=zeros(sA-2);
m=1;
for n1=1:N,
    c=C(eidx(n1)+d);
    c=c(c>0);
    if isempty(c),
        C(eidx(n1))=m;P(idx(n1))=m;m=m+1;
    elseif ~any(diff(c))
        C(eidx(n1))=c(1);
    end
end
D=zeros(size(A));D(idx)=C(eidx);

% %find lowest-valued neighbours
% n=(1:N)';
% B=nan(sA);C=B;
% for n1=1:numel(d), B(idx+d(n1))=min(B(idx+d(n1)),n); end
% %B(~isnan(B))=idx(B(~isnan(B)));
% %assigns labels
% m=1;
% for n1=1:N,
%     if isnan(C(idx(B(idx(n1))))), C(idx(n1))=m;m=m+1; 
%     else C(idx(n1))=C(idx(B(idx(n1)))); end;
% end
