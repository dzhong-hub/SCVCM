function [amgtws, amgstd, amgws, amstd, itera] = temporalAssimilateMasconBySCVCM2(mgtws, mgstd, mttws, mtstd, mgws, mstd,weight)
% function for temporal assimilating GRACE TWSA into LSM TWSA by the iterative SCVCM
%
% Input:
%   mgtws -- the GRACE TWSA value
%   mgstd -- the standard deviation of the GRACE TWSA value
%   mttws -- the EALCO TWSA data
%   mtstd -- the standard deviation of the EALCO TWSA data
%   mgws  -- the GSWSA
%   mstd  -- the standard deviation of the GSWSA
%   weight -- option for weighting method
%             1 - equal weighting
%             2 - weighting with estimated variances
%   
% Output:
%   amgtws -- the assimilated TWSA
%   amgstd -- the standard deviation of the assimilated TWSA
%   amgws  -- the estimated groudwater storage
%   amstd  -- the standard deviation of the estimated groudwater storage
%   itera  -- the total iteration number
%
% Detang Zhong, CCRS, NRCAn 2020-07-20

%% assign the maximal iteration number
max_itera = 5;

%% assign observations
lsm = mttws; clear mttws;
grc = [mgtws; mean(mgtws)]; clear mgtws;
gws = mgws; clear mgws;

%% build up the coefficient matrices         
m = length(lsm);
so = [1, 1, 1, 1];

%% build up the variances for the stochastic model 
if weight == 1 % apply the equal weighting 
    Vlsm = ones(m,1)*(mean(mtstd)^2);
    Vgrc = ones(m+1,1)*mean(mgstd)^2;
    Vgws = ones(m,1)*mean(mstd)^2;
else
    Vlsm = mtstd.^2;
    Vgrc = [mgstd.^2;mean(mgstd)^2];
    Vgws = mstd.^2;
end

%% start iteration loop
loop_continue = 1;
xo = [grc(1:m,1); gws(1:m,1)];

%% coefficient matrices
Alsm = [eye(m), -eye(m)];
Agrc = [eye(m), zeros(m,m); ones(1,m)/m, zeros(1,m)];
Agws = [zeros(m,m)  eye(m)];
itera = 1;

while loop_continue  
    %% observations
    llsm = lsm - Alsm*xo;
    lgrc = grc - Agrc*xo;
    lgws = gws - Agws*xo;

    %% weight matrices
    Plsm = eye(m);
    for i = 1:m
        Plsm(i,i)=1/(Vlsm(i)+0.0001);
    end
    Pgrc = eye(m+1);
    for i = 1:m+1
        Pgrc(i,i)=1/Vgrc(i);
    end
    Pgws = eye(m);
    for i = 1:m
        Pgws(i,i)=1/Vgws(i);
    end

    %% normal equation
    % Nlsm = Alsm'*Plsm*Alsm; blsm = Alsm'*Plsm*llsm;
    % Ngrc = Agrc'*Pgrc*Agrc; bgrc = Agrc'*Pgrc*lgrc;
    % Ngws = Agws'*Pgws*Agws; bgws = Agws'*Pgws*lgws;    
    % N = Nlsm+Ngrc+Ngws; b = blsm+bgrc+bgws;   
    N = Alsm'*Plsm*Alsm + Agrc'*Pgrc*Agrc + Agws'*Pgws*Agws; 
    b = Alsm'*Plsm*llsm + Agrc'*Pgrc*lgrc + Agws'*Pgws*lgws;
    % Q = inv(N); 
    % dx = Q*b;
    dx = N\b;
    x = xo + dx;
    %% residuals
    vlsm = Alsm*dx - llsm;
    vgrc = Agrc*dx - lgrc;
    vgws = Agws*dx - lgws;

    %% pvv
    pvvlsm = vlsm'*Plsm*vlsm;
    pvvgrc = vgrc'*Pgrc*vgrc;
    pvvgws = vgws'*Pgws*vgws;

    %% redudancy
    rlsm = trace(eye(m)-Alsm*(N\(Alsm'*Plsm)));
    rgrc = trace(eye(m+1)-Agrc*(N\(Agrc'*Pgrc)));
    rgws = trace(eye(m)-Agws*(N\(Agws'*Pgws)));
    % rlsm = trace(eye(m)-Alsm*(Q*(Alsm'*Plsm)));
    % rgrc = trace(eye(m+1)-Agrc*(Q*(Agrc'*Pgrc)));
    % rgws = trace(eye(m)-Agws*(Q*(Agws'*Pgws)));

    %% variance component estimation
    so2 = (pvvlsm+pvvgrc+pvvgws)/(rlsm+rgrc+rgws);
    slsm2 = pvvlsm/rlsm;
    sgrc2 = pvvgrc/rgrc;
    sgws2 = pvvgws/rgws;

    %% check convergency
    s = sqrt([so2, slsm2, sgrc2, sgws2]);
    max_ds = max(abs(s - so));
    max_dx = max(abs(dx));
    if (max_ds < 0.05 || max_dx < 1) || itera >= max_itera
        loop_continue = 0;
    else
        itera = itera + 1;
        xo = x;
        % so = s;
        Vlsm = slsm2*Vlsm;
        Vgrc = sgrc2*Vgrc;
        Vgws = sgws2*Vgws;
    end    
end

%% get the solutions for output
Q = inv(N);
sx = sqrt(so2*diag(Q));
amgtws = x(1:m);   
amgws = x(m+1:2*m);
amgstd = sx(1:m);
amstd = sx(m+1:2*m);




