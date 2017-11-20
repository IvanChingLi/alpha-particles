format shortg

data1=importdata('0_5fine100xcourse6Torr.dat');
data2=importdata('0_5fine100xcourse8Torr.dat');
data3=importdata('0_5fine100xcourse12Torr.dat');
data4=importdata('0_5fine100xcourse25Torr.dat');
data5=importdata('0_5fine100xcourse46Torr.dat');
data6=importdata('0_5fine100xcourse102Torr.dat');
data7=importdata('0_5fine100xcourse149Torr.dat');
data8=importdata('0_5fine100xcourse199Torr.dat');
data9=importdata('0_5fine100xcourse299Torr.dat');
data10=importdata('0_5fine100xcourse348Torr.dat');
data11=importdata('0_5fine100xcourse405Torr.dat');
data12=importdata('0_5fine100xcourse449Torr.dat');

counts_data=[data1(:,2) data2(:,2) data3(:,2) data4(:,2) data5(:,2) data6(:,2) data7(:,2) data8(:,2) data9(:,2) data10(:,2) data11(:,2) data12(:,2)];

channel=data1(:,1);
p1 =    0.001262  ; %calibration parameters
p2 =      -1.871  ;
energy=p1*channel+p2;

subtract1=data1(:,2);
subtract11=data11(:,2);

for i=1:11
    channel=data1(:,1);
    counts=counts_data(:,i);
    plot(energy,counts,'k')
    hold on
    f = fit(channel,counts,'gauss1');
    lowerbound=max(round(f.b1-5*f.c1),1);
    upperbound=round(f.b1+5*f.c1);
    channel=channel(lowerbound:upperbound);
    counts=counts(lowerbound:upperbound);
    if i<6
        subtract_counts=subtract11(lowerbound:upperbound);
    else
        subtract_counts=subtract1(lowerbound:upperbound);
    end
    int_counts(i)=sum(counts)-sum(subtract_counts);
    [f gof] = fit(channel,counts,'gauss1');
    peakloc(i)=f.b1;
    peakloc_err(i)=f.c1/sqrt(2);
    adjrsq(i)=gof.adjrsquare;
end

pressure=[6 8 12 25 46 102 149 199 299 348 405];

plot(pressure,int_counts,'ko')


output=[peakloc' peakloc_err' pressure' adjrsq']

%%%%%%%%%% plot and fit energy vs. pressure %%%%%%%
p1 =    0.001262  ;
p1_err=max(p1-0.001193, 0.00133-p1);
p2 =      -1.871  ;
p2_err=max(p2-(-2.27), -1.473-p2);

for i=1:11
    channel=output(i,1);
    channel_err=output(i,2);
    energy(i)=p1*channel+p2;
    energy_err(i)=sqrt((p1*channel_err)^2+(channel*p1_err)^2+(p2_err)^2);
    weights(i)=1/energy_err(i)^2;
end
% ploterr(pressure,energy,[],energy_err,'ko')
% hold on


% f=fit(pressure',energy','poly2','Weights',weights)
% 
% %%%%%%%% find dE/dP %%%%%%
% p1=f.p1;
% p1_err=max(p1-(-1.705e-05), -1.111e-05-p1);
% p2=f.p2;
% p2_err=max(p2-(-0.007714), -0.005093-p2);

% for i=1:12
%     dEdP(i)=2*p1*pressure(i)+p2;
%     dEdP_err(i)=sqrt((2*pressure(i)*p1_err)^2+(p2_err)^2);
% end
% 
% ploterr(pressure,-dEdP,[],dEdP_err,'k.')
% %%%%%%%% find -dE/dx %%%%%%%%
%%%%%%%% assume r=5.5cm %%%%%%%
r=5.5;
% 
% for i=1:12
%     neg_dEdx(i)=dEdP(i)*2*pressure(i)/r;
%     neg_dEdx_err(i)=dEdP_err(i)*2*pressure(i)/r;
% end
% plot(energy,neg_dEdx,'ko')

%%%%%% above failed...try range instead, remember to RESET constant c NOT to zero %%%%%%
fun = @(E,P,c) - c*E./(P.*log((10^5)*E));

    %%%% input is MeV and Torr %%
    %% convert to SI units so output is m!! %%
MeV2J=1.6e-13;
Torr2Pa=133;
c=7.01513e46;
for i=1:11
    range(i) = integral(@(x)fun(x,pressure(i)*Torr2Pa,c),5.486*MeV2J,energy(i)*MeV2J);
    upperbound = integral(@(x)fun(x,pressure(i)*Torr2Pa,c),5.486*MeV2J,(energy(i)+energy_err(i))*MeV2J);
    lowerbound = integral(@(x)fun(x,pressure(i)*Torr2Pa,c),5.486*MeV2J,(energy(i)-energy_err(i))*MeV2J);
    range_err(i)=max(abs(upperbound-range(i)),abs(range(i)-lowerbound));
end

% MeV2J=1;
% Torr2Pa=1;
% c=1;
% for i=1:12
%     range(i) = integral(@(x)fun(x,pressure(i)*Torr2Pa,c),5.486*MeV2J,energy(i)*MeV2J);
%     upperbound = integral(@(x)fun(x,pressure(i)*Torr2Pa,c),5.486*MeV2J,(energy(i)+energy_err(i))*MeV2J);
%     lowerbound = integral(@(x)fun(x,pressure(i)*Torr2Pa,c),5.486*MeV2J,(energy(i)-energy_err(i))*MeV2J);
%     range_err(i)=max(upperbound,lowerbound);
% end

ploterr(pressure,range,[],range_err,'ko')

