function [ dsfmobj ] = dsfm_iv(IV,Km,Kt,Zeta_start,dim_mon,dim_ttm,ikmon,ikttm,tol,maxiter,h,beta1)
%THE CODE IS "CRUDE": USE AS A SCRIPT, NOT AS FUNCTION OF CREATE YOUR
%VERSION

% num_digw = 3;
% IV = round(IV.*(10^num_digw))./(10^num_digw);

IV = IVDATA;

ttm_data    = IV(:,1);
mon_data    = IV(:,2);
date_data   = IV(:,4);

Dates       = unique(date_data);
DatNum      = datenum(date_data); 
matnum      = [0; find(diff(DatNum) ~= 0)];
T           = length(matnum);

%delete the duplicates
IVCLEAN = [];

for i = 1:T
    %obtain data for each J_t
    if i == length(matnum)
        ind    = matnum(i)+1:length(mon_data);
        TTM_Jt = ttm_data(ind);
    else
        ind    = matnum(i)+1:matnum(i+1);
        TTM_Jt = ttm_data(ind);
    end
    
    DataMtx  = IV(ind,:); 
    intmtn   = [0; find(diff(TTM_Jt) ~= 0)];
    DUMTX    = [];
    
    for j = 1:length(intmtn)
    
        if j == length(intmtn)
            PaDaMtx = DataMtx(intmtn(j)+1:end,:);
        else
            PaDaMtx = DataMtx(intmtn(j)+1:intmtn(j+1),:);
        end
        
        %extract the unique matrix for the step i
        ColStrike      = 5;  
        DataCol        = PaDaMtx(:,ColStrike);
        [~,indun,~]    = unique(DataCol);
        DataUMtx       = PaDaMtx(indun,:);
        DUMTX          = [DUMTX;DataUMtx];
    end
    
    IVCLEAN = [IVCLEAN;DUMTX];
end

IV          = IVCLEAN;

ttm_data    = IV(:,1);
mon_data    = IV(:,2);
iv_data     = IV(:,3);
date_data   = IV(:,4);
strike_data = IV(:,5);
pricu_data  = IV(:,6);
rate_data   = IV(:,7)./100;
prico_data  = IV(:,8);

Dates       = unique(date_data);
DatNum      = datenum(date_data); 
matnum      = [0; find(diff(DatNum) ~= 0)];
T           = length(matnum);

%marginally transformed data
ttm_data    = ksdensity(ttm_data,ttm_data,'function','cdf');
mon_data    = ksdensity(mon_data,mon_data,'function','cdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mon_grid = unique(mon_data);
%ttm_grid = unique(ttm_data);
%iv_grid  = unique(iv_data);

%OBTAIN THE SCALED MONEYNESS
%(1) Calculate average volatilities
IVT = [];
for i = 1:T
    if i == T
        TTM_Jt = ttm_data(matnum(i)+1:end);
        IV_Jt  = iv_data(matnum(i)+1:end);
        ttmnum = [0; find(diff(TTM_Jt) ~= 0)];
        Tt     = length(ttmnum);
        ivt = [];
        for j = 1:Tt            
            if j == Tt
                iv_Tt  = IV_Jt(ttmnum(j)+1:end);
                mivt   = mean(iv_Tt);
                repivt = ones(length(iv_Tt),1).*mivt;
            else
                iv_Tt  = IV_Jt(ttmnum(j)+1:ttmnum(j+1));
                mivt   = mean(iv_Tt);
                repivt = ones(length(iv_Tt),1).*mivt;
            end
            ivt = [ivt;repivt];
        end        
    else       
        TTM_Jt = ttm_data(matnum(i)+1:matnum(i+1));
        IV_Jt  = iv_data(matnum(i)+1:matnum(i+1));
        ttmnum = [0; find(diff(TTM_Jt) ~= 0)];
        Tt     = length(ttmnum);
        ivt = [];
        for j = 1:Tt            
            if j == Tt
                iv_Tt  = IV_Jt(ttmnum(j)+1:end);
                mivt   = mean(iv_Tt);
                repivt = ones(length(iv_Tt),1).*mivt;
            else
                iv_Tt  = IV_Jt(ttmnum(j)+1:ttmnum(j+1));
                mivt   = mean(iv_Tt);
                repivt = ones(length(iv_Tt),1).*mivt;
            end
            ivt = [ivt;repivt];
        end

    end
    IVT = [IVT;ivt];
end

%(2) Do moneyness scaling
sc_mon_data = exp( (-beta1*(beta1 - 1).*IVT.*ttm_data) + beta1.*log(mon_data) );
sc_mon_data = ksdensity(sc_mon_data,sc_mon_data,'function','cdf');

mon_data    = sc_mon_data;

%Produce knots' sequences for B-spline estimation
%knots_mon      = aptknt(mon_grid,K);
%knots_mon      = [min(mon_data),min(mon_data),min(mon_data):10:max(mon_data),max(mon_data),max(mon_data)];
knots_mon = [min(mon_data),min(mon_data),linspace(min(mon_data),max(mon_data),ikmon),max(mon_data),max(mon_data)];
knots_ttm = [min(ttm_data),min(ttm_data),linspace(min(ttm_data),max(ttm_data),ikttm),max(ttm_data),max(ttm_data)];
%knots_ttm      = aptknt(ttm_grid,K);


%KK = (length(knots_mon)-K-1)*(length(knots_ttm)-K-1); 
KK = (length(knots_mon)-Km)*(length(knots_ttm)-Kt);

PHI     = cell(1,T);
YY      = cell(1,T);
JT      = zeros(1,T); 

for i = 1:T
    %obtain data for each J_t
    if i == length(matnum)
        ind    = matnum(i)+1:length(mon_data);
        MON_Jt = mon_data(ind);
        Jt     = length(MON_Jt);
        TTM_Jt = ttm_data(ind);
        IV_Jt  = iv_data(ind);
    else
        ind    = matnum(i)+1:matnum(i+1);
        Jt     = length(ind);
        MON_Jt = mon_data(ind);
        TTM_Jt = ttm_data(ind);
        IV_Jt  = iv_data(ind);
    end
    JT(i) = Jt;
    
    %estimate splines at x points: MATLAB BUILT-IN SPCOL AND THE HUNYADI
    %BSPLINE_BASIS FUNCTION GIVE THE SAME ANSWERS, BUT SPCOL IS MUCH FASTER!!!
%     Phi     = zeros(KK, Jt);
%     for j = 1:Jt
%         MonBas = zeros(1,length(knots_mon)-Km);
%         for k = 1:length(knots_mon)-Km
%             MonBas(k) = bspline_basis(k-1,Km,knots_mon,MON_Jt(j)); 
%         end
%         TtmBas = zeros(1,length(knots_ttm)-Kt);
%         for k = 1:length(knots_ttm)-Kt
%             TtmBas(k) = bspline_basis(k-1,Kt,knots_ttm,TTM_Jt(j)); 
%         end
%         
%         tenscomb = combvec(MonBas,TtmBas);
%         tenscol  = prod(tenscomb,1);
%         Phi(:,j) = tenscol';
%     end
    
    %estimate splines
    UMON_Jt = unique(MON_Jt);
    UTTM_Jt = unique(TTM_Jt);
    MonMat  = spcol(knots_mon,Km,UMON_Jt);
    TtmMat  = spcol(knots_ttm,Kt,UTTM_Jt);
    Phi     = zeros(KK, Jt);
    for j = 1:Jt
        ind_mon  = find(UMON_Jt == MON_Jt(j));
        ind_ttm  = find(UTTM_Jt == TTM_Jt(j));
        combin   = combvec(MonMat(ind_mon,:), TtmMat(ind_ttm,:)); %#ok<FNDSB>
        Phi(:,j) = prod(combin,1);       
    end
    PHI{i} = Phi;
    YY{i}  = IV_Jt;
end



%Zeta       = [ones(1,T);Zeta_start'];
L          = size(Zeta_start,1);
ZETA_start = reshape(Zeta_start,L*T,1).*100000;

%set up the initial guess solution
alpha_start = rand(KK*(L+1),1);
SOL_OLD = [alpha_start;ZETA_start];

diffV = 10; % make sure the loop runs at least once
iter  = 0;

while (diffV > tol)
    
    if iter > maxiter
        error('dsfm_iv:maxiter','The algorithm does not converge with the given maximum number of iterations')
        %break
    end
    
    alpha = SOL_OLD(1:KK*(L+1),1);
    ZETA  = SOL_OLD(KK*(L+1)+1:end,1);
    Zeta  = [ones(1,T);reshape(ZETA,L,T)];
    
    Alpha = reshape(alpha,L+1,KK);
    AA    = Alpha(2:end,:); 
    first_el     = zeros(KK*(L+1),T);
    secd_el      = zeros(KK*(L+1),T);
    first_el_noa = zeros(KK*(L+1),KK*(L+1),T);
    %F01          = zeros(T*(L+1),1);
    F01          = cell(T,1);
    %jump_length  = 0;
    F02          = zeros(T*L);
    II = [zeros(1,L);eye(L)];
    F11          = cell(1,T);
    for i = 1:T
        first_el(:,i)       = kron( PHI{i}*PHI{i}', Zeta(:,i)*Zeta(:,i)' )*alpha;
        secd_el(:,i)        = kron( PHI{i}*YY{i}, Zeta(:,i)); 
        first_el_noa(:,:,i) = kron( PHI{i}*PHI{i}', Zeta(:,i)*Zeta(:,i)' );

    %     F01(jump_length+1:jump_length+L,1) = ( Zeta(:,i)'*Alpha*PHI{i}*PHI{i}'*AA' - YY{i}'*PHI{i}'*AA' )';
    %     jump_length                        = jump_length + L; 
        F01{i}              = ( Zeta(:,i)'*Alpha*PHI{i}*PHI{i}'*AA' - YY{i}'*PHI{i}'*AA' )';

        f02                 = zeros(T);
        f02(i,i)            = 1;
        f02_block           = AA*PHI{i}*PHI{i}'*AA';
        f02_kron            = kron(f02,f02_block);
        F02                 = F02 + f02_kron;

        F11{i}              = kron( PHI{i}*PHI{i}'*AA',Zeta(:,i) ) + kron( PHI{i}*PHI{i}'*Alpha'*Zeta(:,i), II ) ...
                              - kron( PHI{i}*YY{i}, II);
    end

    F10 = 2.*(sum(first_el,2) - sum(secd_el,2));
    F20 = 2.*(sum(first_el_noa,3));
    F01 = 2.*cell2mat(F01);
    %F02 = F02;
    F11 = 2.*cell2mat(F11);

    FBIG  = [F10;F01];
    DFBIG = [F20, F11;
             F11',F02];
    
    %compute the solution iteration

    SOL_NEW = SOL_OLD - (pinv(DFBIG))*FBIG;
    
    diffV = max(abs(SOL_NEW - SOL_OLD))
    
    SOL_OLD = SOL_NEW;
    
    iter = iter + 1
end



SOL_FIN = SOL_NEW;
ALPHA   = SOL_FIN(1:KK*(L+1),1);
ZZETA   = SOL_FIN(KK*(L+1)+1:end,1);
ALPHA   = reshape(ALPHA,L+1,KK);


ttm_data  = IV(:,1);
mon_data  = IV(:,2);

mongrid   = linspace(min(mon_data),max(mon_data),dim_mon);
ttmgrid   = linspace(min(ttm_data),max(ttm_data),dim_ttm);


MONMAT         = spcol(knots_mon,Km,mongrid);
TTMMAT         = spcol(knots_ttm,Kt,ttmgrid);

COEF             = zeros(length(knots_mon)-Km,length(knots_ttm)-Kt,L+1);
MHAT             = zeros(length(mongrid),length(ttmgrid),L+1);

for i = 1:L+1
    COEF(:,:,i)   = reshape(ALPHA(i,:)',length(knots_mon)-Km,length(knots_ttm)-Kt);
    %obtain the estimated factor functions
    MHAT(:,:,i)   = MONMAT*COEF(:,:,i)*TTMMAT';
end

%Norming and orthogonalization of factor functions mhat and coefficients
%zeta
du             = (mongrid(2)-mongrid(1))*(ttmgrid(2)-ttmgrid(1));

tempmat              = 0*MHAT(:,:,1)+1;
tempmat(2:(end-1),:) = 2*tempmat(2:(end-1),:);
tempmat(:,2:(end-1)) = 2*tempmat(:,2:(end-1));

%Norming matrices 
GAMMA = zeros(L);
gamma = zeros(L,1);


%Numeric integration
for i = 1:L
    gamma(i)             = sum(sum(tempmat.*MHAT(:,:,1).*MHAT(:,:,i+1)))*du/4;
    for j = 1:L
        GAMMA(i,j)             = sum(sum(tempmat.*MHAT(:,:,j+1).*MHAT(:,:,i+1)))*du/4;
    end
end

%Vectorize factor functions
MHATMat             = zeros(size(MHAT(:,:,1),1)*size(MHAT(:,:,1),2),L+1);
for i = 1:(L+1)
      MHATMat(:,i)             = reshape(MHAT(:,:,i),size(MHAT(:,:,1),1)*size(MHAT(:,:,1),2),1); 

end

%Obtain normed coefficients Zeta (as in Fengler, p. 166, eq. 5.74)
Zeta_new             = zeros(L,T);
Zeta_est             = reshape(ZZETA,L,T);
for i = 1:T
    Zeta_new(:,i)             = (GAMMA^0.5)*( Zeta_est(:,i) + (GAMMA^(-1))*gamma );
end

%Normalizing the factor functions
% MHATMatNew = zeros(size(MHAT(:,:,1),1)*size(MHAT(:,:,1),2),L+1);
% for i = 1:L+1
%     if i == 1
%         MHATMatNew(:,i) = MHATMat(:,1)' -gamma'*GAMMA^(-1)*MHATMat(:,2:end)';
%     end
%     MHATMatNew(:,i) = GAMMA^(-0.5)*MHATMat(:,2:end)';   
% end

MHATMatZero  = MHATMat(:,1)' -gamma'*GAMMA^(-1)*MHATMat(:,2:end)';
MHATMatShort = GAMMA^(-0.5)*MHATMat(:,2:end)';

%Create the B matrix for PCA transformation
B                 = Zeta_new*Zeta_new';

[Z,~]             = eigs(B,L);

ZETA_FIN          = zeros(L,T);

for i = 1:T
    ZETA_FIN(:,i) = Z'*Zeta_new(:,i);
end

MHATMatFin = Z'*MHATMatShort;
MHATMatFin = MHATMatFin';


%Obtain final factor functions 
MHAT_FIN             = zeros(size(MHAT(:,:,1),1), size(MHAT(:,:,1),2), L+1);
for i = 1:L+1
    if i == 1
        MHAT_FIN(:,:,i)             = reshape(MHATMatZero, size(MHAT(:,:,1),1), size(MHAT(:,:,1),2));
        continue
    end
    MHAT_FIN(:,:,i)             = reshape(MHATMatFin(:,i-1), size(MHAT(:,:,1),1), size(MHAT(:,:,1),2));  
end



%PRODUCE THE DYNAMIC IV SURFACES
COEFF        = zeros(length(knots_mon)-Km,length(knots_ttm)-Kt,L+1);
COEFF(:,:,1) = ((MONMAT'*MONMAT)^(-1))*MONMAT'*MHAT_FIN(:,:,1)*TTMMAT*(TTMMAT'*TTMMAT)^(-1);
COEFF(:,:,2) = ((MONMAT'*MONMAT)^(-1))*MONMAT'*MHAT_FIN(:,:,2)*TTMMAT*(TTMMAT'*TTMMAT)^(-1);
COEFF(:,:,3) = ((MONMAT'*MONMAT)^(-1))*MONMAT'*MHAT_FIN(:,:,3)*TTMMAT*(TTMMAT'*TTMMAT)^(-1);
COEFF(:,:,4) = ((MONMAT'*MONMAT)^(-1))*MONMAT'*MHAT_FIN(:,:,4)*TTMMAT*(TTMMAT'*TTMMAT)^(-1);

cc = zeros(length(knots_mon)-Km,length(knots_ttm)-Kt,L+1);
spfm = cell(1,L+1);
for l = 1:L+1
    cc(:,:,l) = COEFF(:,:,l);
    spfm{l} = spmak({knots_mon,knots_ttm},cc(:,:,l));
end

%NOW PUT IN "TRUE" GRIDS

monsurf  = zeros(length(mongrid),length(ttmgrid),L+1);
mon2surf = zeros(length(mongrid),length(ttmgrid),L+1);
ttmsurf  = zeros(length(mongrid),length(ttmgrid),L+1);
for l = 1:L+1
    for i = 1:length(mongrid)
        for j = 1:length(ttmgrid)
            monsurf(i,j,l)  = fnval(fnder(spfm{l},[1,0]),[mongrid(i);ttmgrid(j)] );
            mon2surf(i,j,l) = fnval(fnder(spfm{l},[2,0]),[mongrid(i);ttmgrid(j)] );
            ttmsurf(i,j,l)  = fnval(fnder(spfm{l},[0,1]),[mongrid(i);ttmgrid(j)] );
        end
    end
end


IV_DYN             = zeros(size(MHAT(:,:,1),1), size(MHAT(:,:,1),2), T);

for t = 1:T
    IV_DYN(:,:,t)    = MHAT_FIN(:,:,1);
    for l = 1:L
        IV_DYN(:,:,t) = IV_DYN(:,:,t) + ZETA_FIN(l,t).*MHAT_FIN(:,:,l+1);
    end      
end


D_IV_mon   = zeros(size(MHAT(:,:,1),1), size(MHAT(:,:,1),2), T);
D_IV_ttm   = zeros(size(MHAT(:,:,1),1), size(MHAT(:,:,1),2), T);
D2_IV_mon  = zeros(size(MHAT(:,:,1),1), size(MHAT(:,:,1),2), T);
for t = 1:T
    D_IV_mon(:,:,t)  = monsurf(:,:,1);
    D_IV_ttm(:,:,t)  = ttmsurf(:,:,1);
    D2_IV_mon(:,:,t) = mon2surf(:,:,1);

    for l = 1:L
        
        D_IV_mon(:,:,t)  = D_IV_mon(:,:,t) + ZETA_FIN(l,t).*monsurf(:,:,l+1); 
        D_IV_ttm(:,:,t)  = D_IV_ttm(:,:,t) + ZETA_FIN(l,t).*ttmsurf(:,:,l+1);
        D2_IV_mon(:,:,t) = D2_IV_mon(:,:,t) + ZETA_FIN(l,t).*mon2surf(:,:,l+1);
    end      
end

%local volatility for ETF options
if nargin < 11
    beta = 1;
else
    beta = beta1;
end

LV_DYN = zeros(size(MHAT(:,:,1),1), size(MHAT(:,:,1),2), T);
d1     = zeros(size(MHAT(:,:,1),1), size(MHAT(:,:,1),2), T);
d2     = zeros(size(MHAT(:,:,1),1), size(MHAT(:,:,1),2), T);
for t = 1:T 
    
    for i = 1:length(mongrid)
        for j = 1:length(ttmgrid)
            d1(i,j,t)     = ( -log(mongrid(i)) + 0.5*(abs(beta)^2)*(IV_DYN(i,j,t)^2)*ttmgrid(j) )/(beta*IV_DYN(i,j,t)*sqrt(ttmgrid(j))); 
            d2(i,j,t)     = d1(i,j,t) - abs(beta)*IV_DYN(i,j,t)*sqrt(ttmgrid(j));
            
            LV_DYN(i,j,t) = ( IV_DYN(i,j,t)^2 + 2*ttmgrid(j)*IV_DYN(i,j,t)*D_IV_ttm(i,j,t) )/( 1 + 2*...
                abs(beta)*mongrid(i)*sqrt(ttmgrid(j))*d1(i,j,t)*D_IV_mon(i,j,t) + (beta^2)*...
                (mongrid(i)^2)*ttmgrid(j)*( d1(i,j,t)*d2(i,j,t)*D_IV_mon(i,j,t)^2 + IV_DYN(i,j,t)*D2_IV_mon(i,j,t) ) );
        end
    end
    
end



dsfmobj.mhat = MHAT_FIN;
dsfmobj.zhat = ZETA_FIN;
dsfmobj.iv   = IV_DYN;
dsfmobj.date = Dates;
dsfmobj.TTM  = ttmgrid;
dsfmobj.MON  = mongrid;

dsfmobj.dmon  = D_IV_mon;
dsfmobj.dttm  = D_IV_ttm;
dsfmobj.d2mon = D2_IV_mon;

dsfmobj.lv    = LV_DYN;

% if MON_Jt(j) < int_knot{k}{1}(2) && MON_Jt(j) > int_knot{k}{1}(1) && ...
%                     TTM_Jt(j) < int_knot{k}{2}(2) && TTM_Jt(j) > int_knot{k}{2}(1)
%                 %obtain the value of the Phi matrix
%                 phi_j(k) = spval(sp,{MON_Jt(j),TTM_Jt(j)}); 
%                 %spcol(knots_mon,K,mon_grid)*coefs*spcol(knots_ttm,K,ttm_grid).'
%                 else
%                 phi_j(k) = 0;
% end

end


