function [G1,G2] = G_Estimate(lrms,P,opts)
% Inputs:
%           lrms:            MS image;
%           P:               PAN image;
%           opts:            Parameters.
% Outputs:
%           G1_1,G2_1:       Adaptive coefficients.
% 
% Reference:
% [1] J.-L. Xiao, T.-Z. Huang, L.-J. Deng, Z.-C. Wu and G. Vivone, 
% A New Context-Aware Details Injection Fidelity with Adaptive Coefficients
% Estimation for Variational Pansharpening,
% IEEE Trans. Geosci. Remote Sens., doi:10.1109/TGRS.2022.3154480.
%==========================================================================

%%  Initialization
Nways = opts.Nways;
sensor = opts.sensor;
sz = size(P);
YY=lrms;
sf = opts.sf;
G1=zeros(Nways);
G2=zeros(Nways);
nclusters=opts.nclusters;
for i=1:Nways(3)
    G1_1 = zeros(sz);
    G2_1 = zeros(sz);
    YYB = interp23tap(YY,sf);
    YYB1 = cat(3,YYB,P);
    PB = MTF_PAN(P,sensor,sf);
    Q = reshape(YYB1,Nways(1)*Nways(2),Nways(3)+1);
    I = kmeans(Q,nclusters);% cluster context-based regions
    I_mat = reshape(I,sz);
    %% Estimation by regression
    Y=lrms(:,:,i);
    XB = interp23tap(Y,sf);
    GXB1 = gradient(XB,1);
    GXB2 = gradient(XB,3);
    GPB1 = gradient(PB,1);
    GPB2 = gradient(PB,3);
    for ii = 1:nclusters
        P_Detail = GPB1(I_mat == ii);
        X_Detail = GXB1(I_mat == ii);
        h = robustfit(P_Detail,X_Detail,'ols');
        A = repmat(h(2,1),sz);
        G1_1 = (I_mat == ii) .* A + (I_mat ~= ii) .*G1_1; 
    end
    for ii = 1:nclusters
        P_Detail = GPB2(I == ii);
        X_Detail = GXB2(I == ii);
        h = robustfit(P_Detail,X_Detail,'ols');
        A = repmat(h(2,1),sz);
        G2_1 = (I_mat == ii) .* A + (I_mat ~= ii) .*G2_1; 
    end
    G1(:,:,i)=G1_1;
    G2(:,:,i)=G2_1;
end
end

