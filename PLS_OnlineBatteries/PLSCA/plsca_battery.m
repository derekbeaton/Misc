function [X_bsrs,Y_bsrs,Y_boots,Y_bsrs_from_boot,comp_perms,omni_perm,l,t,di,fi,ri,ci,dj,fj,rj,cj,fii,fjj,w,m,Smat,R] = plsca_battery(X,Y,iters,DESIGN)
%function [X_bsrs,Y_bsrs,comp_perms,omni_perm,l,t,di,fi,ri,ci,dj,fj,rj,cj,fii,fjj,w,m,Smat,R] = plsca_battery(X,Y,iters,DESIGN)

    if exist('iters') ~=1;iters=100;end
    if exist('DESIGN') ~=1;DESIGN=0;end

    [l,t,di,fi,ri,ci,dj,fj,rj,cj,fii,fjj,m,w,Smat,R] = plsca(X,Y);

    X_bsrs = zeros(size(X,2),length(l));
    X_bsrs_mean = X_bsrs;
    X_bsrs_M2 = X_bsrs;

    Y_bsrs = zeros(size(Y,2),length(l));
    Y_bsrs_mean = Y_bsrs;
    Y_bsrs_M2 = Y_bsrs;    

    Y_boots = zeros(size(Y,2),length(l),iters); %use this to compute more accurate BSRs for Y and to get confidence intervals.


    comp_perms = zeros(iters,length(l));
    omni_perm = zeros(iters,1);
    for i=1:iters
        %%bootstrap as supplemental.
        %boot_indices = 1:size(Y,1); %use this to test that project is
        %correct.
        
        if DESIGN==0
            boot_indices = randsample(size(X,1),size(X,1),'true');
        else 
            boot_indices = createNewY(DESIGN);
        end
        
        BootY = Y(boot_indices,:);
        BootX = X(boot_indices,:);
        BootR = BootX'*BootY;
 
        X_sup = afc_sup(BootR,fj,l);
        Y_sup = afc_sup(BootR',fi,l);
        
        if isempty(find(isnan(X_sup)))~=1 ;  %%replace NaNs
            X_sup(find(isnan(X_sup))) = 0;
        end       
        if isempty(find(isnan(Y_sup)))~=1 ; %%replace NaNs
            Y_sup(find(isnan(Y_sup))) = 0;
        end       

        X_delta = X_sup - X_bsrs_mean;
        X_bsrs_mean = X_bsrs_mean + (X_delta/i);
        X_bsrs_M2 = X_bsrs_M2 + X_delta .* (X_sup-X_bsrs_mean);
        %%the problem: when X_bsrs_M2 is still 0 after the bootstraps are
        %%done, it will go to Inf or -Inf, depending on sign of Q vectors.

        Y_delta = Y_sup - Y_bsrs_mean;
        Y_bsrs_mean = Y_bsrs_mean + (Y_delta/i);
        Y_bsrs_M2 = Y_bsrs_M2 + Y_delta .* (Y_sup-Y_bsrs_mean);

        Y_boots(:,:,i) = Y_sup;

        comp_perms(i,:) = plsca(X,Y(randsample(size(Y,1),size(Y,1)),:));
        omni_perm(i,:) = sum(comp_perms(i,:));
    end
    X_variance = X_bsrs_M2/iters;
    X_bsrs = X_bsrs_mean./sqrt(X_variance);  %if the top is not zero, and the bottom is, it's Inf
    %this is only a problem if the data are very sparse. That is, when
    %columns are relatively empty, such as seen in the example data.

    Y_variance = Y_bsrs_M2/iters;
    Y_bsrs = Y_bsrs_mean./sqrt(Y_variance);

        
    Mean_BootY = mean(Y_boots,3);
    S_BootY = mean( (Y_boots - repmat(Mean_BootY,[1,1,iters])).^2,3).^(1/2);
    Y_bsrs_from_boot = Mean_BootY ./ S_BootY ;
     
end