function [cvar circvar R Rd Y]=showvar_tmp(y,fs,display)
% Estimation of the von Mises parameters
framelength=24/1000*fs;
inc=framelength/8;
win=@blackman;
rNdft=2;
Ndftmax=1024;
Ndft=256;

framelength=round(framelength/2-0.5)*2+1;

if nargin>5
    [f0,t0,p]=fxpefac(s,fs,inc/fs);
else
    [f0,t0,p]=fxpefac(y,fs,inc/fs);
end
f0s=interp1(t0*fs,f0,1:length(y),'linear','extrap');
f0s(f0s<50)=50;
flm=round(framelength/(fs/median(f0s)));    % framelength in periods
t=(framelength+1)/2:inc:(length(y)-(framelength+1)/2);
%
Y=zeros(Ndftmax,length(t));
rY=zeros(2,length(t));

% analysis
% first frame
idx=1:(t(1)*2-1);
rY(:,1)=[idx(1),idx(end)];
i=1;
Y(1:length(idx),i)=fft(y(idx).*win(length(idx)));
for i=2:length(t)-1
    idxspan=round((flm*fs/f0s(t(i))+1)/2)*2-1;
    % limit idxspan
    idxspan=min(idxspan,framelength*2+1);
    idxspan=max(idxspan,(framelength+1)/2);
    rY(:,i)=[max(t(i)-(idxspan-1)/2,1),min(t(i)+(idxspan-1)/2,length(y))];
    yc=y(rY(1,i):rY(2,i));
    Y(1:rNdft*(diff(rY(:,i))+1),i)=fft(yc.*win(diff(rY(:,i))+1),(diff(rY(:,i))+1)*rNdft)...
        .*delay2spec(diff(rY(:,i))/2,(diff(rY(:,i))+1)*rNdft);
end
% last frame
idx=(length(y)-((length(y)-t(end))*2-1)):length(y);
rY(:,end)=[idx(1),idx(end)];
Y(1:length(idx),end)=fft(y(idx).*win(length(idx)));
% end of analysis

% linear phase
xp=(1+flm*rNdft):flm*rNdft:(size(Y,1));
temp=zeros(size(Y,1),1);
temp(round((1+flm*rNdft)/2):size(Y,1))=interp1(xp,1:length(xp),round((1+flm*rNdft)/2):size(Y,1),'nearest','extrap');
z=filter(1,[1 -1],f0s*2*pi/fs);
linphasehm=temp*interp1((0:length(y)-1),z,t,'linear');

%unwrapped phase
R=wrap(angle(Y)-linphasehm);

% detrend phase with a mean of 40ms
Nd=10/1000;
Nd=round((Nd*fs/inc+1)/2)*2-1;
temp2=filter(ones(Nd,1)/Nd,1,exp(1i*[R(1+flm*rNdft,:),zeros(1,(Nd-1)/2)]),[],2);
temp2=angle(temp2(:,(Nd+1)/2:end));
Rd=wrap(angle(Y)-linphasehm-temp*temp2);

% circular mean and variance
rkappa=40/1000; % 40ms range
Nr=round((rkappa*fs/inc+1)/2)*2-1;
circ=filter(ones(Nr,1)/Nr,1,exp(1i*[Rd,zeros(size(Y,1),(Nr-1)/2)]),[],2);
circmean=angle(circ(:,(Nr+1)/2:end));
circvar=1-abs(circ(:,(Nr+1)/2:end));

% interpolation in frequency
cvar=zeros(Ndft/2+1,length(t));
for i=1:length(t)
    Nidx=min((fs/2/f0s(round(t(i)))*(xp(1)-1)),size(circmean,1)-1);
    cvar(:,i)=interp1(((0:round(Nidx))/2)/Nidx,circvar(1:round(Nidx)+1,i),((0:Ndft/2))/Ndft,'linear','extrap');
end

if display
    imagesc(flipud(cvar));
    xticklabels = t/fs;
    % xticklabels
    % xticks = linspace(1, size(cvar, 2), numel(xticklabels));
    xticks_my = [0:length(xticklabels)/3.75:length(xticklabels)];
    xticklabels_my = [0 0.5 1 1.5];
    set(gca, 'XTick', xticks_my, 'XTickLabel', xticklabels_my);
    xlabel('Time (s)');
    yticklabels = fs/2/1000:-1:0;
    yticks = 1:size(cvar,1)/4:size(cvar,1);
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels(:));
    ylabel('Frequency (kHz)');
    % title('Circular Variance');
end
end