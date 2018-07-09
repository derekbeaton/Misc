function [l,t,di,fi,ri,ci,dj,fj,rj,cj,m,w,Y] = corresp(X,nfk) ;
%usage:  [l,t,di,fi,ri,ci,dj,fj,rj,cj] = Mcorresp(X,nfact2keep) ;
% Multiple Correspondence Analysis:
%           X: matrix to analyze
%  nfact2keep: # of factors (def=all)
% l is the eigenvalue vector,
% t is the percentage of inertia vector
%  ***  The following matrices represent the rows -> ending i ***
% fi is the matrix of the row-coordinates
% di is the vector of the (Chi-squared) distance to the centroid
% ri is the matrix of the Correlation between the i set and the axis
% ci is the matrix of the Contributions
%  ***  The following matrices represent the columns -> ending j ***
% fj is the matrix of the column-coordinates
% dj is the vector of the (Chi-squared) distance to the centroid
% rj is the matrix of the Correlation between the j set and the axis
% cj is the matrix of the Contributions
%%%    Hervé Abdi April 2004
%  See also mca for multiple correspondence analysis

%% Compute CA as a bilinear model

le_flip=0;
[I,J]=size(X);
if J<I; X=X';le_flip=1;[I,J]=size(X);end
if exist('nfk') ~=1;nfk=I;end
xtot=sum(sum(X)');
xpj=sum(X);
xip=sum(X,2);
c=sum(X)/xtot;
m=sum(X')/xtot;
w=ones(1,J) ./ c ;
Y= (X./(xip*ones(1,J)))-ones(I,1)*c;
%% use an eigenvalue decomposition to save memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[P,l]=eigen(((Y.*repmat(w,I,1))*Y').*( (m.^(1/2)')*(m.^(1/2)) ) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nf=length(l);all_l=l;
if nf > nfk;nf=nfk;
     P=P(:,1:nf);l=l(1:nf);
end
P= repmat((m'.^(-1/2)),1,nf).*P;
d=l.^(1/2);
fi=P.*repmat(d',I,1);
t=(l/sum(all_l))*100;
di=(Y.^2)*w';
ri=repmat((1./di),1,nf ).*( fi.^2);
ci=repmat(m',1,nf).*(fi.^2)./repmat(l',I,1);
% Compute the solution for the J set using the transition formula
  Z=(X./repmat(xpj,I,1)  )';
  fj=Z*P;
  dj=( (Z-repmat(m,J,1)).^2)*(ones(1,I)./m)';
  rj=repmat((1./dj),1,nf ).*( fj.^2);
  cj=repmat(c',1,nf).*(fj.^2)./repmat(l',J,1);
%
%unflip
if le_flip==1;
    tm=fi;fi=fj;fj=tm;
    tm=di;di=dj;dj=tm;
    tm=ri;ri=rj;rj=tm;
    tm=ci;ci=cj;cj=tm;
end

%%%%%%%%%%%%%%% private functions here
%%%%%%%%%%%%%%% Eigen
