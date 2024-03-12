function varargout=spm_ss_glm(option,varargin)
% SPM_SS_GLM general linear model

% [B,rss]=spm_ss_glm('estimate',X,Y) 
% Fits the multi-variate model Y = X*B + E   (with noise E following a multivariate normal distribution)
%
% [H,f,p,dof]=spm_ss_glm('evaluate',X,B,rss,C,M) 
% Test the hypothesis: H=C*B*M'=0;
%
%  X: (n x r) matrix
%  Y: (n x m) matrix
%  C: (k x r) matrix
%  M: (j x m) matrix
%
%  B: (r x m) estimated regressor matrix
%  rss: (m x m) residual sum squares matrix
%
%  H: (k x j) estimated contrast matrix
%  f: F/T/X2 statistic
%  p: false-positive value
%  dof: degrees of freedom
%

switch(lower(option)),
    case 'estimate',
        x=varargin{1};y=varargin{2};        
        ix=pinv(x'*x);
        b=ix*(x'*y);
        e=y-x*b;
        ee=e'*e; 
        varargout={b,ee};
    case 'evaluate',
        x=varargin{1};b=varargin{2};ee=varargin{3};
        if nargin<5||isempty(varargin{4}), C=eye(size(x,2));else C=varargin{4}; end
        if nargin<6||isempty(varargin{5}), M=eye(size(b,2));else M=varargin{5}; end
        
        if nargin<7||isempty(varargin{6}), Ny=size(x,1);    else Ny=varargin{6}; end
        if nargin<8||isempty(varargin{7}), Nx=rank(x);      else Nx=varargin{7}; end
        if nargin<9||isempty(varargin{8}), Nc0=rank(x*C');  else Nc0=varargin{8}; end
        if nargin<10||isempty(varargin{9}), Nh=[size(C,1),size(M,1)];    else Nh=varargin{9}; end
        if nargin<11||isempty(varargin{10}), if numel(M)==1,Ns=1;else Ns=rank(M);end; else Ns=varargin{10}; end
    
        ix=pinv(x'*x);
        dof=Ny-Nx;
        h=C*b*M';
        r=C*ix*C';
        ee2=M*ee*M';
        
        if Nh(1)==1&&Nh(2)==1,% T-stats
            k=sqrt(r*ee2);
            F=real(h./max(eps,k))*sqrt(dof);
            if isnan(F)||dof<=0,F=nan;p=nan;else p=1-spm_Tcdf(F,dof);end
            edof=dof;
        elseif Nh(1)>1&&Nh(2)==1 % F-stats
            bb=h'*pinv(r)*h;
            F=real(bb./ee2)*dof/Nc0;
            if isnan(F)||dof<=0||Nc0<=0,F=nan;p=nan;else p=1-spm_Fcdf(F,Nc0,dof); end
            edof=[Nc0,dof];
        elseif Nh(1)==1&&Nh(2)>1 % F-stats
            bb=h*pinv(ee2)*h';
            F=real(bb./r)*(dof-Ns+1)/Ns;
            if isnan(F)||Ns<=0||dof-Ns+1<=0,F=nan;p=nan;else p=1-spm_Fcdf(F,Ns,dof-Ns+1); end
            edof=[Ns,dof-Ns+1];
        elseif Nh(1)>1&&Nh(2)>1 % Wilk's Lambda stat,Bartlett's Chi2
            bb=h'*pinv(r)*h;
            F=-(dof-1/2*(Ns-Nc0+1))*real(log(real(det(ee2)./det(ee2+bb))));
            if isnan(F)||Nc0<=0||Ns<=0,F=nan;p=nan;else p=1-spm_Xcdf(F,Ns*Nc0); end
            edof=[Ns*Nc0];
        end
        if nargout>4,
            stderr=sqrt(diag(r)*diag(ee2).')/sqrt(dof);
            varargout={h,F,p,edof,stderr};
        else
            varargout={h,F,p,edof};
        end
end

end
