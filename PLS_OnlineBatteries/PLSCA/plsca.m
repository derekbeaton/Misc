function [l,t,di,fi,ri,ci,dj,fj,rj,cj,fii,fjj,m,w,Smat,R] = plsca(X,Y)
    
    R = X'*Y;
    [l,t,di,fi,ri,ci,dj,fj,rj,cj,m,w,Smat] = corresp(R);
    fii = afc_sup(X,fi,l);
    fjj = afc_sup(Y,fj,l);

end