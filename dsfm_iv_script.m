%%
%INPUT DATA
clc
clear

load IVDATAFULLLAST

%change these parameters to your purpose
IV         = IVDATA;
Km         = 3;
Kt         = 3;
dim_mon    = 60;
dim_ttm    = 60;
ikmon      = 5;
ikttm      = 3;
tol        = 1e-06;
maxiter    = 100000;
h          = 1e-06;

L          = 4;
T          = 195;

%create starting values for Z
CoefMat2    = [0.95  0;
               0.2  -0.3];

CoefMat3    = [0.95 -0.2 0;
               0     0.8 0.1;
               0.1   0   0.6];
          
CoefMat4    = [0.8 -0.2 0 0.3;
              0 0.8 0.1 -0.3 ;
              0.1 0 -0.3 0.2;
              -0.7 -0.5 0.6 0];
          
CoefMat5    = [0.95 -0.2 0 0.3 -0.25;
              0 0.8 0.1 0 -0.3 ;
              0.1 0 -0.3 0.2 0;
              -0.1 0.5 0 0.4 0
              0.1 0 0.6 0.7 0];
          
CoefMat6    = [0.85 -0.2  0    0.3 -0.25 0;
               0     0.8  0.1 -0.5 -0.3  0.7;
              -0.5   0   -0.3  0.2  0.1 -0.1;
              -0.4   0.5 -0.3  0.4  0    0.3;
               0.1  -0.1  0.6  0.7 -0.6  0;
              -0.3   0  0.7  0    0.8 -0.7];
          
CoefMat7    = [0.85 -0.2 0   0.3 -0.25 0.1 0.3;
               0     0.8 0.1 -0.5   -0.3  0.7   0.1;
              0.1 0 -0.3 0.2 0.1 -0.1 0.2;
              -0.1 0.5 -0.3 0.4 0 0.3 -0.1;
              0.1 -0.9 0.6 0.7 -0.6 0 -0.3;
              -0.7 0.1 0.2 0 -0.2 0 0.6;
              0.2 0.8 -0.3 0.4 0.8 -0.7 0];
          
          

%CoefMat = randn(L);
% up_lim  = 1;
% low_lim = -1;
% rng(0,'twister');
% CoefMat3 = (up_lim-low_lim).*rand(L,L) + low_lim

% CoefMat3(1,3) = 0
% CoefMat3(1,5) = 0
% CoefMat3(1,7) = 0
% CoefMat3(2,1) = 0
% CoefMat3(2,2) = 0
% CoefMat3(2,6) = 0
% CoefMat3(3,3) = 0
% CoefMat3(3,6) = 0
% CoefMat3(4,1) = 0
% CoefMat3(4,5) = 0
% CoefMat3(4,7) = 0
% CoefMat3(5,3) = 0
% CoefMat3(5,6) = 0
% CoefMat3(6,2) = 0
% CoefMat3(6,4) = 0
% CoefMat3(7,1) = 0
% CoefMat3(7,3) = 0
% CoefMat3(7,7) = 0

% num_digw = 2
% CoefM3 = round(CoefMat3.*(10^num_digw))./(10^num_digw);

varspec    = vgxset('n',L,'nAR',1,'AR',CoefMat3,'Q',10e-5.*eye(L));
Zeta_start = vgxsim(varspec,T);

Cmat       = eye(L) - (1/L)*ones(L,1)*ones(1,L);
Zeta_start = Zeta_start*Cmat;

Zeta_start = Zeta_start';

% r=rank(Zeta_start);
% [Q,R,E]=qr(Zeta_start);
% newZeta_start=[Zeta_start;transpose(E(:,end-r+1:end))];

Zeta_start(2,50) = -0.0037;

rank(Zeta_start)

%%
%COEFFICIENTS' PLOTTING
linewidth = 2;
plot(Dates,ZETA_FIN(1,:), 'LineWidth',linewidth)
hold on
datetick('x','mmmyy','keepticks')
plot(Dates,ZETA_FIN(2,:),'Color','r', 'LineWidth',linewidth)
plot(Dates,ZETA_FIN(3,:),'Color',[0,.5,0], 'LineWidth',linewidth)
%%
%MONEYNESS' STRINGS' PLOTTING
start_point = 100;
date_data = IV(:,4);
DatNum    = datenum(date_data); 
MATNUM    = [0; find(diff(DatNum) ~= 0)];

for i = start_point:T
    %obtain data for each J_t
    if i == length(MATNUM)
        MON_Jt = mon_data(MATNUM(i)+1:end);
        Jt     = length(MON_Jt);
        TTM_Jt = ttm_data(MATNUM(i)+1:end);
        IV_Jt  = iv_data(MATNUM(i)+1:end);
    else
        Jt     = length(MATNUM(i)+1:MATNUM(i+1));
        MON_Jt = mon_data(MATNUM(i)+1:MATNUM(i+1));
        TTM_Jt = ttm_data(MATNUM(i)+1:MATNUM(i+1));
        IV_Jt  = iv_data(MATNUM(i)+1:MATNUM(i+1));
    end
    M      = MON_Jt;       % moneyness (the larger, the less in the money)
    IVV    = IV_Jt;        % Implied Volatility
    TT     = TTM_Jt;
    matnum = [0;find(diff(TT) ~= 0)];
    IVMmat = {};
    for j = 1:length(matnum)
        if j == length(matnum)
            mmat   = M(matnum(j)+1:end);
            ivmat  = IVV(matnum(j)+1:end);
            ivmmat = [mmat,ivmat];
            ivmmat = ivmmat(~any(isnan(ivmmat),2),:);
            IVMmat = [IVMmat,ivmmat];
            continue
        end
        mmat   = M(matnum(j)+1:matnum(j+1));
        ivmat  = IVV(matnum(j)+1:matnum(j+1));
        ivmmat = [mmat,ivmat];
        ivmmat = ivmmat(~any(isnan(ivmmat),2),:);
        IVMmat = [IVMmat,ivmmat];
    end

    Tpoints = unique(TT);
    fig = figure;
    
    for j = 1:length(Tpoints)
        plot3(Tpoints(j).*ones(length(IVMmat{1,j}),1), IVMmat{1,j}(:,1),...
            IVMmat{1,j}(:,2), '.', 'Color', 'r', 'MarkerSize',15)
        grid on
        hold on
        set(gca,'YDir','reverse');
    end
    Date = datevec(Dates(i));
    if numel(num2str(Date(2))) < 2
        if numel(num2str(Date(3))) < 2         
            title(strcat('IV ticks on', {' '}, num2str(Date(1)), '0', num2str(Date(2)), '0', num2str(Date(3)) ))
        else
            title(strcat('IV ticks on', {' '}, num2str(Date(1)), '0', num2str(Date(2)), num2str(Date(3)) ))
        end
    else
        if numel(num2str(Date(3))) < 2         
            title(strcat('IV ticks on', {' '}, num2str(Date(1)), num2str(Date(2)), '0', num2str(Date(3)) ))
        else
            title(strcat('IV ticks on', {' '}, num2str(Date(1)), num2str(Date(2)), num2str(Date(3)) ))
        end
    end
    
    xlabel('TTM');
    ylabel('Moneyness');
    zlabel('IV');
    print (fig, '-r600', '-depsc', strcat('ivtick',num2str(i-start_point+1)));
    eps2pdf(strcat('ivtick',num2str(i-start_point+1),'.eps'), strcat('ivtick',num2str(i-start_point+1),'.pdf'), 1, 0)
    hold off
end




%%
%SURFACES' PLOTTING
tpoint = 100;
Date = Dates(tpoint);
values = IV_DYN(:,:,tpoint);

figure
surf(mongrid,ttmgrid,values')
hold on
xlabel('Moneyness')
ylabel('TTM')
zlabel('IV')
title( strcat('IV surface on', {' '}, datestr(Date)) )
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FOR MONEYNESS-SCALED IV ONLY!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ttmgrid  = linspace(min(IV(:,1)),max(IV(:,1)),dim_ttm);
mongrid  = linspace(min(IV(:,2)),max(IV(:,2)),dim_mon);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%correct to "scivsurf" for moneyness-scaled IV

start_point = 129;

for i = start_point:T
    fig = figure;
    mesh(mongrid,ttmgrid,IV_DYN(:,:,i)')
    hold on
    xlabel('Moneyness')
    ylabel('TTM')
    zlabel('IV')
    Date = datevec(Dates(i));
    if numel(num2str(Date(2))) < 2
        if numel(num2str(Date(3))) < 2         
            title(strcat('IV surface on', {' '}, num2str(Date(1)), '0', num2str(Date(2)), '0', num2str(Date(3)) ))
        else
            title(strcat('IV surface on', {' '}, num2str(Date(1)), '0', num2str(Date(2)), num2str(Date(3)) ))
        end
    else
        if numel(num2str(Date(3))) < 2         
            title(strcat('IV surface on', {' '}, num2str(Date(1)), num2str(Date(2)), '0', num2str(Date(3)) ))
        else
            title(strcat('IV surface on', {' '}, num2str(Date(1)), num2str(Date(2)), num2str(Date(3)) ))
        end
    end
    print (fig, '-r600', '-depsc', strcat('scivsurf',num2str(i-start_point+1)));
    eps2pdf(strcat('scivsurf',num2str(i-start_point+1),'.eps'), strcat('scivsurf',num2str(i-start_point+1),'.pdf'), 1, 0)
    hold off
end

values = D_IV_mon(:,:,15);
surf(mongrid(1:30),ttmgrid(1:30),values(1:30,1:30)')
surf(mongrid,ttmgrid,MHAT_FIN(:,:,6)')
%%
%LVS plotting
values2 = (LV_DYN(:,:,148));
% Rcut = 38;
% Ccut = 54;
surf(mongrid,ttmgrid,((values2)') )

lcut = 5;
rcut = 55;
start_point = 100;

for i = start_point:T
    fig = figure;
    values2 = LV_DYN(:,:,i);
    surf(mongrid(lcut:rcut),ttmgrid(lcut:rcut),(sqrt(values2(lcut:rcut,lcut:rcut))'))
    hold on
    xlabel('Moneyness')
    ylabel('TTM')
    zlabel('LV')
    Date = datevec(Dates(i));
    if numel(num2str(Date(2))) < 2
        if numel(num2str(Date(3))) < 2         
            title(strcat('LV surface on', {' '}, num2str(Date(1)), '0', num2str(Date(2)), '0', num2str(Date(3)) ))
        else
            title(strcat('LV surface on', {' '}, num2str(Date(1)), '0', num2str(Date(2)), num2str(Date(3)) ))
        end
    else
        if numel(num2str(Date(3))) < 2         
            title(strcat('LV surface on', {' '}, num2str(Date(1)), num2str(Date(2)), '0', num2str(Date(3)) ))
        else
            title(strcat('LV surface on', {' '}, num2str(Date(1)), num2str(Date(2)), num2str(Date(3)) ))
        end
    end
    print (fig, '-r600', '-depsc', strcat('sclvsurf',num2str(i-start_point+1-2)));
    eps2pdf(strcat('sclvsurf',num2str(i-start_point+1-2),'.eps'), strcat('sclvsurf',num2str(i-start_point+1-2),'.pdf'), 1, 0)
    hold off
end
%%
%RV CRITERION COMPUTATION
% ZZETA = ZETA_FIN;
% ZZETA = [ones(1,T);ZZETA];
ZZETA = SOL_FIN(KK*(L+1)+1:end,1);
ZZETA = [ones(1,T);reshape(ZZETA,L,T)];
YHAT  = zeros(T,max(JT));
for i = 1:T
    for j = 1:JT(i)
        YHAT(i,j) = ZZETA(:,i)'*(ALPHA*PHI{i}(:,j));
    end
end

YHAT = reshape(YHAT',T*max(JT),1);
ind  = find(YHAT ~= 0);
YHAT = YHAT(ind);

RV   = sum( (iv_data - YHAT).^2 )/sum( (iv_data - mean(iv_data)).^2 );

%EXPLAINED VARIATION
ExplVar = 1 - RV

%%
%FACTOR FUNCTIONS' PLOTTING
figure
surf(mongrid,ttmgrid,MHAT_FIN(:,:,4)')
hold on
xlabel('Moneyness')
ylabel('TTM')
%zlabel('IV')
title( '$$\hat{m}_3$$', 'Interpreter', 'Latex' )
hold off
%%
%STRINGS' FITTING
upper = 8425;
lower = 8474;
plot( sort( mon_data(upper:lower),'descend' ),sort( iv_data(upper:lower), 'ascend'),'.','Color','r','MarkerSize',25)
hold on
plot(sort( mon_data(upper:lower),'descend' ),sort( YHAT(upper:lower) ,'ascend'),'LineWidth',3)

figure
plot( ( mon_data(upper:lower) ),( iv_data(upper:lower)),'.','Color','r','MarkerSize',25)
hold on
plot(( mon_data(upper:lower)),( YHAT(upper:lower) ),'LineWidth',3)
xlabel('Moneyness')
ylabel('IV')
title( 'DSFM individual string fit 20150326, 108 days to exp.' )
hold off

%%
%NEGATIVE LV INTERPOLATION ATTEMPT
lcut = 6;
rcut = 54;
val = LV_DYN(:,:,120);
V = val;
V(V<0) = V(V<0)*0;
[row,col] = find(V == 0);
X   = mongrid;
Y   = ttmgrid;
Vq  = interp2(X,Y,V,X(row),Y(col),'cubic');

for i = 1:length(row)
    for j = 1:length(col)
        V(row(i),col(j)) = Vq(i);
    end
end

figure
surf(X,Y,sqrt(V)')

Xq = linspace(min(X),max(X),100);
Yq = linspace(min(Y),max(Y),100);

[Xq,Yq] = meshgrid(Xq,Yq);

VVq = interp2(X,Y,V,Xq,Yq,'cubic');

figure
mesh(Xq,Yq,VVq);
hold on
colormap summer

%%
%DYNAMIC LV PLOTTING (CORRECT DERIVATIVES)

cc      = COEF(:,:,1);
spfm1   = spmak({knots_mon,knots_ttm},cc);
dermat1 = fnder(spfm{1},[1,0])
ncc     = [dermat1.coefs(:,:,1);dermat1.coefs(:,:,2);dermat1.coefs(:,:,3);dermat1.coefs(:,:,4)]'
ders = fnval(fndir(spfm1,[2 0;0 0]),[mongrid(10);ttmgrid(10)]);
spfm2   = spmak({dermat1.knots{1},dermat1.knots{2}},ncc)

timep = 120;
monsurf  = zeros(length(mongrid),length(ttmgrid),L+1);
mon2surf = zeros(length(mongrid),length(ttmgrid),L+1);
ttmsurf  = zeros(length(mongrid),length(ttmgrid),L+1);
for l = 1:L+1
    for i = 1:length(mongrid)
        for j = 1:length(ttmgrid)
            monderiv(i,j)  = fnval(fnder(spfm2,[1,0]),[mongrid(i);ttmgrid(j)] );
            %mon2surf(i,j,l) = fnval(fnder(spfm1,[2,0]),[mongrid(i);ttmgrid(j)] );
            %ttmsurf(i,j,l)  = fnval(fnder(spfm1,[0,1]),[mongrid(i);ttmgrid(j)] );
        end
    end
end

%dersurf = dersurf1 + dersurf2 + dersurf3 + dersurf4;
monsurf  = monsurf1 + ZETA_FIN(1,timep).*monsurf2 + ZETA_FIN(2,timep).*monsurf3 + ZETA_FIN(3,timep).*monsurf4;
mon2surf = mon2surf1 + ZETA_FIN(1,timep).*mon2surf2 + ZETA_FIN(2,timep).*mon2surf3 + ZETA_FIN(3,timep).*mon2surf4;
ttmsurf  = ttmsurf1 + ZETA_FIN(1,timep).*ttmsurf2 + ZETA_FIN(2,timep).*ttmsurf3 + ZETA_FIN(3,timep).*ttmsurf4;

D_IV_mon(:,:,timep)  = monsurf;
D_IV_ttm(:,:,timep)  = ttmsurf; 
D2_IV_mon(:,:,timep) = mon2surf;

vval     = dersurf;
surf(mongrid,ttmgrid,ans)

dermat1 = fnder(spfm1,[1,0]);
nkn_mon = knots_mon(2:end);
nkn_ttm = knots_ttm(2:end);
MM      = spcol(nkn_mon,Km,mongrid);
TT      = spcol(nkn_ttm,Kt,ttmgrid);
ncc     = [dermat1.coefs(:,:,1);dermat1.coefs(:,:,2);dermat1.coefs(:,:,3);dermat1.coefs(:,:,3)]';

MH1     = MM*ncc*TT';

Rcut    = 22;
Ccut    = 10;
val     = MH1(1:end-Rcut,1:end-Ccut);

surf(mongrid(1:end-Rcut),ttmgrid(1:end-Ccut),val')

km = 3;
kt = 3;

LVD_FIN = cell(1,T);
LVD_NN  = cell(1,T);
TTIND   = cell(1,T);
MV      = cell(1,T);
TV      = cell(1,T);

for i = 100:T
    lvd        = LV_DYN(:,:,i);
    lvd        = lvd(1:64,1:90);
    ttind      = find(any(lvd < 0) == 0);
    TTIND{i}   = ttind;
    lvd        = lvd(:,ttind);
    LVD_NN{i}  = lvd;
    mm         = mongrid(1:size(lvd,1));
    tt         = ttmgrid(1:size(lvd,2));
    knotst     = augknt(linspace(min(tt),max(tt),5),kt);
    sp         = spap2(knotst,kt,tt,sqrt(lvd));
    coefst     = fnbrk(sp,'coefs');
    %vals   = fnval(sp,tt);
    %mesh(mm,tt,vals.') 
    knotsm     = augknt(linspace(min(mm),max(mm),5),km);
    sp2        = spap2(knotsm,km,mm,coefst.');
    coefs      = fnbrk(sp2,'coefs').';
    mv         = linspace(min(mm),max(mm),50); 
    tv         = linspace(min(tt),max(tt),50);
    MV{i}      = mv;
    TV{i}      = tv;
    vals       = spcol(knotsm,km,mv)*coefs*spcol(knotst,kt,tv).';
    LVD_FIN{i} = vals;
end

%the same for the IMPLIED VOLATILITY

IVD_FIN = cell(1,T);
IVD_NN  = cell(1,T);
TTIND   = cell(1,T);
MV      = cell(1,T);
TV      = cell(1,T);

for i = 100:T
    ivd        = IV_DYN(:,:,i);
    ivd        = ivd(1:68,1:90);
    ttind      = find(any(ivd < 0) == 0);
    TTIND{i}   = ttind;
    ivd        = ivd(:,ttind);
    IVD_NN{i}  = ivd;
    mm         = mongrid(1:size(ivd,1));
    tt         = ttmgrid(1:size(ivd,2));
    knotst     = augknt(linspace(min(tt),max(tt),5),kt);
    sp         = spap2(knotst,kt,tt,sqrt(ivd));
    coefst     = fnbrk(sp,'coefs');
    knotsm     = augknt(linspace(min(mm),max(mm),5),km);
    sp2        = spap2(knotsm,km,mm,coefst.');
    coefs      = fnbrk(sp2,'coefs').';
    mv         = linspace(min(mm),max(mm),50); 
    tv         = linspace(min(tt),max(tt),50);
    MV{i}      = mv;
    TV{i}      = tv;
    vals       = spcol(knotsm,km,mv)*coefs*spcol(knotst,kt,tv).';
    IVD_FIN{i} = vals;
end

start_point = 129;

for i = start_point:T
    fig = figure;
    %values2 = LV_DYN_NN_restr{i};
    mesh(MV{i},TV{i},LVD_FIN{i}')
    hold on
    grid off
    xlabel('Moneyness')
    ylabel('TTM')
    zlabel('LV')
    Date = datevec(Dates(i));
    if numel(num2str(Date(2))) < 2
        if numel(num2str(Date(3))) < 2
            title(strcat('LV surface on', {' '}, num2str(Date(1)), '0', num2str(Date(2)), '0', num2str(Date(3)) ))
        else
            title(strcat('LV surface on', {' '}, num2str(Date(1)), '0', num2str(Date(2)), num2str(Date(3)) ))
        end
    else
        if numel(num2str(Date(3))) < 2
            title(strcat('LV surface on', {' '}, num2str(Date(1)), num2str(Date(2)), '0', num2str(Date(3)) ))
        else
            title(strcat('LV surface on', {' '}, num2str(Date(1)), num2str(Date(2)), num2str(Date(3)) ))
        end
    end
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');

    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    print ( fig, strcat('smsclvsurf',num2str(i-start_point+1)), '-dpdf' );
    %eps2pdf(strcat('smlvsurf',num2str(i-start_point+1-2),'.eps'), strcat('smlvsurf',num2str(i-start_point+1-2),'.pdf'), 1, 0)
    hold off
end

%plotting the SMOOTHED IMPLIED VOLATILITY now

for i = start_point:T
    fig = figure;
    %values2 = LV_DYN_NN_restr{i};
    mesh(MV{i},TV{i},IVD_FIN{i}')
    hold on
    xlabel('Moneyness')
    ylabel('TTM')
    zlabel('LV')
    Date = datevec(Dates(i));
    if numel(num2str(Date(2))) < 2
        if numel(num2str(Date(3))) < 2
            title(strcat('LV surface on', {' '}, num2str(Date(1)), '0', num2str(Date(2)), '0', num2str(Date(3)) ))
        else
            title(strcat('LV surface on', {' '}, num2str(Date(1)), '0', num2str(Date(2)), num2str(Date(3)) ))
        end
    else
        if numel(num2str(Date(3))) < 2
            title(strcat('LV surface on', {' '}, num2str(Date(1)), num2str(Date(2)), '0', num2str(Date(3)) ))
        else
            title(strcat('LV surface on', {' '}, num2str(Date(1)), num2str(Date(2)), num2str(Date(3)) ))
        end
    end
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');

    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    print ( fig, strcat('smivsurf',num2str(i-start_point+1)), '-dpdf' );
    %eps2pdf(strcat('smlvsurf',num2str(i-start_point+1-2),'.eps'), strcat('smlvsurf',num2str(i-start_point+1-2),'.pdf'), 1, 0)
    hold off
end


mons  = zeros(length(mongrid),length(ttmgrid),L+1);
for l = 1:L+1
    for i = 1:length(mongrid)
        for j = 1:length(ttmgrid)
            mons(i,j,l)  = fnval(spfm{l},[mongrid(i);ttmgrid(j)] );
%             mon2surf(i,j,l) = fnval(spfm{l},[mongrid(i);ttmgrid(j)] );
%             ttmsurf(i,j,l)  = fnval(spfm{l},[mongrid(i);ttmgrid(j)] );
        end
    end
end

mm = mongrid(1:68);
tt = ttmgrid(1:90);
km = 3;
kt = 3;
knotst = augknt(linspace(min(tt),max(tt),5),kt);
sp = spap2(knotst,kt,tt,sqrt(val));
coefst = fnbrk(sp,'coefs');
vals = fnval(sp,tt);
mesh(mm,tt,vals.') 
knotsm = augknt(linspace(min(mm),max(mm),5),km);
sp2 = spap2(knotsm,km,mm,coefst.');
coefs = fnbrk(sp2,'coefs').';
mv = linspace(min(mm),max(mm),100); 
tv = linspace(min(tt),max(tt),100);
values = spcol(knotsm,km,mv)*coefs*spcol(knotst,kt,tv).';
mesh(mv,tv,values')

%%
%ARBITRAGE CHECK WITH TOTAL VARIANCE CROSSING
date_data = IV(:,4);
DatNum    = datenum(date_data); 
MATNUM    = [0; find(diff(DatNum) ~= 0)];

tpoint = 117;

if tpoint == length(MATNUM)
    MON_Jt = mon_data(MATNUM(tpoint)+1:end);
    Jt     = length(MON_Jt);
    TTM_Jt = ttm_data(MATNUM(tpoint)+1:end);
    IV_Jt  = iv_data(MATNUM(tpoint)+1:end);
else
    Jt     = length(MATNUM(tpoint)+1:MATNUM(tpoint+1));
    MON_Jt = mon_data(MATNUM(tpoint)+1:MATNUM(tpoint+1));
    TTM_Jt = ttm_data(MATNUM(tpoint)+1:MATNUM(tpoint+1));
    IV_Jt  = iv_data(MATNUM(tpoint)+1:MATNUM(tpoint+1));
end
M      = MON_Jt;       % moneyness (the larger, the less in the money)
IVV    = IV_Jt;        % Implied Volatility
TT     = TTM_Jt;
matnum = [0;find(diff(TT) ~= 0)];
IVMmat = {};
for j = 1:length(matnum)
    if j == length(matnum)
        mmat   = M(matnum(j)+1:end);
        ivmat  = IVV(matnum(j)+1:end);
        ivmmat = [mmat,ivmat];
        ivmmat = ivmmat(~any(isnan(ivmmat),2),:);
        IVMmat = [IVMmat,ivmmat];
        continue
    end
    mmat   = M(matnum(j)+1:matnum(j+1));
    ivmat  = IVV(matnum(j)+1:matnum(j+1));
    ivmmat = [mmat,ivmat];
    ivmmat = ivmmat(~any(isnan(ivmmat),2),:);
    IVMmat = [IVMmat,ivmmat];
end

Tpoints = unique(TT);
fig     = figure;
cmap    = hsv(10);
%cols    = ['b', 'r', 'g', 'm'];

for j = 1:length(Tpoints)
    
    plot(IVMmat{1,j}(:,1),IVMmat{1,j}(:,2), '.', 'Color', cmap(j,:), 'MarkerSize',15)
    grid on
    hold on
    %set(gca,'YDir','reverse');
    
end

