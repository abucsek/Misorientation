% Calculates misorientation axis from spot misorientations
% Needs as many spot info as possible - in this example there are six (all I could get with monoclinic crap data)
% Spot input info from unit pole figures
%
% NOTE: You need to be consistent on which is the "spot end" and "spot
% start" (i.e., order matters). It's not always easy to tell since
% different (hkl)'s will flip them. If certain (hkl)s are totally changing
% your answer when you include them, it's probably b/c you need to switch
% "start" and "end. 
%
% NOTE: You can use radius,eta,omega instead of the unit pole figure info.
% Just make eta,omega your inputs instead of something calculated from the
% unit pole figures. 
%
% NOTE: I calculate the "n" (terminology from Pagan paper) using the g
% vector corresponding to the spot maximum or centroid. This is safer than
% using a nominal (hkl), since your orientation may be slightly off. This
% way you don't have to fully trust your orientation. '

clear; clc


%%
%% Inputs
%%
Orientation = [0.74219 0.65472 0.14252 0.013634];  % orientation of grain
R = [.0781  .0873  .1406  .1406  .104  .104];  % spot radii (m)
HKL = [1 0 0;  1 0 -1;  -1 -2 1;  -1 2 1;  -1 -1 1;  -1 1 1];  % spot (hkl)'s
% position on unit pole figure corresponding to spot maximum or centroid
spotMax_PF = [-0.9720     -0.2029      -0.1107;  -0.7516   -0.1468   -0.6426;  -0.2834   -0.8528   -0.4385;
    0.6644   -0.6598    0.3511;  -0.5196   -0.6370   -0.5694;  0.7720   -0.3968    0.4966];
% position on unit pole figure corresponding to one "end" of spot
spotStart_PF = [-.9773 -.2077 .04279;  -.8307 -.166 -.5314;  -.3271 -.8703 -.3683
    .6957 -.651 .3035;  -.5897 -.6619 -.4628;  .8311 -.3787 .4073];
% position on unit pole figure corresponding to other "end" of spot
spotEnd_PF = [-.9139 -.181 -.3634;  -.7166 -.1367 -.684;  -.2609 -.848 -.4614;
    .643 -.6682 .3742;  -.471 -.6295 -.618;  .7266 -.4094 .5517];

latticeParams = [2.9011 4.1333 4.64771 90 97.281 90];

D = 1012.36 * 1e-3;  % S2D dist (m)
beamEnergy = 55.618;

% direct lattice vectors
m1 = latticeParams(1) * [1;0;0];
m2 = latticeParams(2) * [0;1;0];
m3 = latticeParams(3) * [cosd(latticeParams(5));0;sind(latticeParams(5))];
[r1,r2,r3] = reciprocal(m1,m2,m3);  % reciprocal lattice vectors
r1r=quat2rot(Orientation)*r1;  r2r=quat2rot(Orientation)*r2;  r3r=quat2rot(Orientation)*r3;  % rotated reciprocal lattice vectors


%%
%% Convert pole figure info to g vectors in the lab coordinate system, then use to calculate n and v (terminology from Pagan paper)
%%

% Calculate incoming x-ray wave vector
h = 6.626E-34;  c = 3.000E+08;  e = 1.602E-19;
Wavelength = h * c / ( 1000 * beamEnergy * e ) * 1e10;
ki_L = 2*pi / Wavelength * (-[0;0;1]);

for gg = 1 : size(spotStart_PF,1)
    
    %% Calculate nominal n vector from grain orientation and (hkl) - used to make sure sign of g vectors are correct
    n_Nominal = unit(HKL(gg,1)*quat2rot(Orientation)*r1r + HKL(gg,2)*quat2rot(Orientation)*r2r + HKL(gg,3)*quat2rot(Orientation)*r3r);
    
    %% Convert spot maximum unit pole figure info to g vector in the lab coordinate system
    % Convert pole figure info --> omega, eta
    [ome, eta, ~] = cart2sph(spotMax_PF(gg,1),spotMax_PF(gg,2),spotMax_PF(gg,3));
    ome = ome - pi;   ome = ome + 2*pi;  % cart2sph has been tricky for me so these are "fixes"
    eta = pi - eta;
    
    % Convert omega, eta, radius --> outgoing x-ray wave vector in the lab coord. sys.
    zeta_x = R(gg) * cos(eta);  zeta_y = R(gg) * sin(eta);
    kounit_L = unit([zeta_x; zeta_y; -D]);
    ko_L = norm(ki_L) * kounit_L;
    
    % Convert outgoing x-ray wave vector --> g vector in the lab coord. sys
    ghkl_L = unit(ko_L - ki_L);
    n = Rodrot(ome,[0;1;0])' * ghkl_L;
    if sign(n(1)) ~= sign(n_Nominal(1))  % fix sign if input (hkl) is actually -(hkl)
        n = -n;
        sig = 1;
    else
        sig = 0;
    end
    
    
    %% Convert start/end spot unit pole figure info to g vectors in the lab coordinate system
    % Convert pole figure info --> omega, eta
    xS = spotStart_PF(gg,1);  yS = spotStart_PF(gg,2);  zS = spotStart_PF(gg,3);
    xE = spotEnd_PF(gg,1);  yE = spotEnd_PF(gg,2);  zE = spotEnd_PF(gg,3);
    for ii = 1:2
        if ii == 1
            x=xS;  y=yS;  z=zS;
        else
            x=xE;  y=yE;  z=zE;
        end
        [ome, eta, ~] = cart2sph(x,y,z);
        ome = ome - pi;  ome = ome + 2*pi;
        eta = pi - eta;
        
        % Convert omega, eta, radius --> outgoing x-ray wave vector in the lab coord. sys.
        zeta_x = R(gg) * cos(eta);  zeta_y = R(gg) * sin(eta);
        kounit_L = unit([zeta_x; zeta_y; -D]);
        ko_L = norm(ki_L) * kounit_L;
        
        % Convert outgoing x-ray wave vector --> g vectora in the lab coord. sys
        ghkl_L = unit(ko_L - ki_L);
        
        % Convert g vectora in the lab coord. sys --> n_E and n_S (terminology from Pagan paper)
        if ii == 1
            nS = Rodrot(ome,[0;1;0])' * ghkl_L;
            if sig == 1
               nS = -nS;
            end
        else
            nE = Rodrot(ome,[0;1;0])' * ghkl_L;
            if sig == 1
               nE = -nE;
            end            
        end
    end
    
    %% Store n and v vectors
    v = (nE - nS);
    vStore(:,gg) = v;
    nStore(:,gg) = n;
end


%%
%% My attempt at a bullshit least squares minimization, works if I run a bunch of times and just use the answer with the smallest "quantityToMinimize"
%%

quantityToMinimize = 100;
for ii = 1 : 10000
    
    % Try random w vector
    h = (rand - 0.5)*2;
    k = (rand - 0.5)*2;
    l = (rand - 0.5)*2;
    wTrial = unit([h k l])*rand*25*pi/180;
    
    % Calculate quantityToMinimize
    trialQuantityToMinimize = 0;
    for gg = 1 : size(spotStart_PF,1)
        % From Pagan paper, App. B
        u = norm(vStore(:,gg) - cross(wTrial,nStore(:,gg))');
        zeta = acos( dot(vStore(:,gg), cross(wTrial,nStore(:,gg))') / (norm(vStore(:,gg)) * norm(cross(wTrial,nStore(:,gg))')) );
        trialQuantityToMinimize = trialQuantityToMinimize + (u^2 + zeta^2);
    end
    if trialQuantityToMinimize < quantityToMinimize
        quantityToMinimize = trialQuantityToMinimize;
        w = wTrial;
    end
end


%% Solutions
quantityToMinimize
unit(w)  % misorientation axis
hkl = (eye(3)/[r1r r2r r3r]) * unit(w)';  h=hkl(1);  k=hkl(2);  l=hkl(3);
hkl_w = unit(hkl)  % misorientation axis in (hkl) form w/ respect to input grain orientation
mori = norm(w) * 180/pi  % misorientation amount
