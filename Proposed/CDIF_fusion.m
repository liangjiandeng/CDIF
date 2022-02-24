function [X_out]     =  CDIF_fusion(Y, P,G1_1,G2_1,opts)
% Inputs:
%           Y:               MS image;
%           P:               PAN image;
%           G1_1,G2_1:       Adaptive coefficients;
%           opts:            Parameters.
% Output:
%           X_out:           Pasharpened image.
% 
% Reference:
% [1] J.-L. Xiao, T.-Z. Huang, L.-J. Deng, Z.-C. Wu and G. Vivone, 
% A New Context-Aware Details Injection Fidelity with Adaptive Coefficients
% Estimation for Variational Pansharpening,
% IEEE Trans. Geosci. Remote Sens., doi:10.1109/TGRS.2022.3154480.
%==========================================================================

%% Initiation
maxit = opts.maxit;
tol   = opts.tol;              
Nways = opts.Nways;
sz    = opts.sz;
lambda = opts.lambda;
lambda2 = opts.lambda2;
lambda3 = opts.lambda3;
eta_1 = opts.eta_1 ;
eta_2 = opts.eta_2 ;
eta_3 = opts.eta_3 ;
eta_4 = opts.eta_4 ;
sf    = opts.sf;
sensor = opts.sensor;
X      = imresize(Y, sf); % Initialize the X
par  = FFT_kernel(sf, sensor ,Nways);
M     = zeros(Nways); 
H1= zeros(Nways);
H2= zeros(Nways);
H3= zeros(Nways);
Thet_1= zeros(Nways);  
Thet_2= zeros(Nways);
Thet_3= zeros(Nways);
Thet_4= zeros(Nways);
%%====Operators========%%
fx=[1,-1];
fy=[1;-1];
otfFx = psf2otf(fx,[Nways(1),Nways(2)]);
otfFy = psf2otf(fy,[Nways(1),Nways(2)]);
eig = abs(otfFx).^2+abs(otfFy).^2;
filter.x(1,:,:) = 1;      filter.x(2,:,:) = -1;
filter.y(:,1,:) = 1;      filter.y(:,2,:) = -1;
filter.z(:,:,1) = 1;      filter.z(:,:,2) = -1;
eigsDxTDx = abs(psf2otf(filter.x,Nways)).^2;
eigsDyTDy = abs(psf2otf(filter.y,Nways)).^2;
eigsDzTDz = abs(psf2otf(filter.z,Nways)).^2;
%%  The ADMM-based solver
X_k = X;
for it = 1:maxit
    % update X
    for sp=1:Nways(3)
        Q1 = G1_1(:,:,sp).*gradient(P,1);
        Q2 = G2_1(:,:,sp).*gradient(P,3);
        K1 = ForwardDxT(H1);
        K2 = ForwardDyT(H2);
        K3 = ForwardDzT(H3);
        P1 = ForwardDxT(Thet_2);
        P2 = ForwardDyT(Thet_3);
        P3 = ForwardDzT(Thet_4);
        FFT_Up = 2*lambda*conj(otfFx).*fft2(Q1)+...
            2*lambda*conj(otfFy).*fft2(Q2)+eta_1*M(:,:,sp).*par.fft_BT(:,:,sp)-fft2(Thet_1(:,:,sp)).*par.fft_BT(:,:,sp)...
            +eta_2*fft2(K1(:,:,sp))+eta_3*fft2(K2(:,:,sp))+eta_4*fft2(K3(:,:,sp))...
            -fft2(P1(:,:,sp))-fft2(P2(:,:,sp))-fft2(P3(:,:,sp));
        FFT_Down = eta_1*par.fft_B(:,:,sp).*par.fft_BT(:,:,sp) + (eta_2)*eigsDxTDx(:,:,sp)+...
            (eta_3)*eigsDyTDy(:,:,sp)+eta_4*eigsDzTDz(:,:,sp)+2*lambda*eig;
        X(:,:,sp) = real(ifft2(FFT_Up./FFT_Down));
    end  
    X(X<0)=0;
    X(X>1)=1;

    % update M
    temp_UP  =  2*par.ST(Y)+eta_1*par.B(X)+Thet_1;

    SST      =  zeros(sz);
    s0       =  3;
    SST(s0:sf:end,s0:sf:end) = ones(sz/sf);

    temp_DOWN=  2*SST+eta_1;
    temp_DOWN=  repmat(temp_DOWN, [1 1 Nways(3)]);
    M = temp_UP./temp_DOWN;
    
    % update H
    H1 = wthresh(ForwardDx(X)+Thet_2/eta_2,'s', lambda2/eta_2);
    H2 = wthresh(ForwardDy(X)+Thet_3/eta_3,'s', lambda2/eta_3);
    H3 = wthresh(ForwardDz(X)+Thet_4/eta_4,'s', lambda3/eta_4);
    
    % update multipliers
    Thet_1 = Thet_1+eta_1*(MTF_im(X, sensor, '', sf)-M); 
    Thet_2 = Thet_2+eta_2*(ForwardDx(X)-H1);
    Thet_3 = Thet_3+eta_3*(ForwardDy(X)-H2);
    Thet_4 = Thet_4+eta_4*(ForwardDz(X)-H3);
   
    % relative errors
    Rel_Err = norm(Unfold(X-X_k,Nways,3) ,'fro')/norm(Unfold(X_k,Nways,3),'fro');
    X_k = X;
    if Rel_Err < tol  
        break;
    end   
end
 X_out=X;
end

