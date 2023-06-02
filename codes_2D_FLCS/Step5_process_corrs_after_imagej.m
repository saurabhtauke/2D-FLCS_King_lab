d11 = c11(:,2);
d12 = c12(:,2);
d21 = c21(:,2);
d22 = c22(:,2);

% d11(26) = NaN;
% d12(26) = NaN;
% d21(26) = NaN;
% d22(26) = NaN;

tot = d11(end)+d12(end)+d21(end)+d22(end);
d11 = d11/tot;
d12 = d12/tot;
d21 = d21/tot;
d22 = d22/tot;

occupancies = [d11(end),d12(end),d21(end),d22(end)]

g11 = d11/d11(end);
g12 = d12/d12(end); 
g21 = d21/d21(end);
g22 = d22/d22(end);
% g22 = zeros(size(d22));

p1 = [1E-6 10E-6];
p2 = [10E-6 20E-6];
p3 = [20E-6 30E-6];
p4 = [30E-6 40E-6];
p5 = [40E-6 50E-6];
p6 = [50E-6 60E-6];
p7 = [60E-6 70E-6];
p8 = [70E-6 80E-6];
p9 = [80E-6 90E-6];
p10 = [90E-6 100E-6];
p11 = [100E-6 200E-6];
p12 = [200E-6 300E-6];
p13 = [300E-6 400E-6];
p14 = [400E-6 500E-6];
p15 = [500E-6 600E-6];
p16 = [600E-6 700E-6];
p17 = [700E-6 800E-6];
p18 = [800E-6 900E-6];
p19 = [900E-6 1E-3];
p20 = [1E-3 2E-3];
p21 = [2E-3 5E-3];
p22 = [5E-3 10E-3];
p23 = [10E-3 20E-3];
p24 = [20E-3 50E-3];
p25 = [50E-3 100E-3];
p26 = [100E-5 200E-3];
p27 = [200E-3 300E-3];
p28 = [300E-3 500E-3];
p29 = [500E-3 700E-3];
p30 = [700E-3 1];
p31 = [1 2];
p32 = [2 3];


frames = [p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11;p12;p13;p14;p15;p16;p17;p18;p19;p20;p21;p22;p23;p24;p25;p26;p27;p28;p29;p30;p31;p32];
    dframe = frames(:,2) - frames(:,1);
dt = mean(frames,2);


figure(619)
semilogx(dt,d11,'LineWidth',2)
hold on
semilogx(dt,d12,'LineWidth',2)
semilogx(dt,d21,'LineWidth',2)
semilogx(dt,d22,'LineWidth',2)
hold off

xlabel('\Delta T')
ylabel('corr')
legend('corr 1-1','corr 1-2', 'corr 2-1,','corr 2-2')
title('egfp-W-16hrs')