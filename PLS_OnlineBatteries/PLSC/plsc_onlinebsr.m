clear
load BADIA_CHAPTER_DATA

Xnorm = zscore(X)/sqrt(size(X,1)-1); %normed
Ynorm = zscore(Y)/sqrt(size(Y,1)-1); %normed
R = Xnorm'*Ynorm; %this is the correlation matrix between X & Y.
[U,s,V] = svd(R,0); %and PLSC

%%and now, the bootstrap, and ratios.
%%this is the online version (Weldford/Knuth algorithm. see Wikipedia).
%%%i.e., a running tally. Confidence intervals cannot be computed anymore.

%step 1: initialize some structures for computing ratios
%X_bsrs = zeros(size(X,2),length(s)); %%this will hold the ratios
X_bsrs_mean = zeros(size(X,2),length(s)); %this holds an online mean
X_bsrs_M2 = zeros(size(X,2),length(s)); %this helps us get the variance.

%Y_bsrs = zeros(size(Y,2),length(s));
Y_bsrs_mean = zeros(size(Y,2),length(s));
Y_bsrs_M2 = zeros(size(Y,2),length(s));
iters=100;
for i=1:iters
    %if DESIGN==0
        indices = randsample(size(X,1),size(X,1),'true');
    %else 
    %    indices = createNewY(DESIGN); %%if you have a design matrix, use
    %    this one so you can resample within groups.
    %end
    
    %BootY = Ynorm(indices,:);
    %BootX = Xnorm(indices,:);
    %BootR = BootX'*BootY;
    BootR = Xnorm(indices,:)'*Ynorm(indices,:);

    %compute BootR as supplementary variables -- skipping some center/scale
    %assumptions for sake of example.
    X_sup = BootR * V;
    Y_sup = BootR' * U;
    
    Xnans = find(isnan(X_sup));
    Ynans = find(isnan(Y_sup));
    if isempty(Xnans)~=1 ;  %%replace NaNs
        X_sup(Xnans) = 0;
    end       
    if isempty(Ynans)~=1 ; %%replace NaNs
        Y_sup(Ynans) = 0;
    end    

    %where the magic happens. 
    X_delta = X_sup - X_bsrs_mean;
    X_bsrs_mean = X_bsrs_mean + (X_delta/i);
    X_bsrs_M2 = X_bsrs_M2 + X_delta .* (X_sup-X_bsrs_mean);

    Y_delta = Y_sup - Y_bsrs_mean;
    Y_bsrs_mean = Y_bsrs_mean + (Y_delta/i);
    Y_bsrs_M2 = Y_bsrs_M2 + Y_delta .* (Y_sup-Y_bsrs_mean);
    
end
%and now compute the BSRatios.
X_variance = X_bsrs_M2/iters;
X_bsrs = X_bsrs_mean./sqrt(X_variance);

Y_variance = Y_bsrs_M2/iters;
Y_bsrs = Y_bsrs_mean./sqrt(Y_variance);

