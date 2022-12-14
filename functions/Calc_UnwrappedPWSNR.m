function [UnRMSE , UnHPSNR]=Calc_UnwrappedPWSNR(RefSpeech,ModSpeech,param)
% calculate only at voicing regions 
framelength = 24*0.001;   % in ms
stepsize =     3*0.001;   % in ms
window=@blackman;
pvthreshold = 0.99;
[xr,t]=enframe(RefSpeech,window(framelength*param.fs),stepsize*param.fs);
[fx, tx, pv]=fxpefac(RefSpeech,param.fs,stepsize);
pv=interp1(tx,pv,t/param.fs,'linear','extrap');

% calculate unwrapped phase
[var_Ref circvar_Ref R_Ref Rd_Ref Y_Ref] = showvar_tmp(RefSpeech,param.fs,0);
[var_Mod circvar_Mod R_Mod Rd_Mod Y_Mod] = showvar_tmp(ModSpeech,param.fs,0);

if size(circvar_Ref,1) ~= size(circvar_Mod,1)
    if size(circvar_Ref,1)<size(circvar_Mod,1)
        circvar_Mod = circvar_Mod(1:size(circvar_Ref,1),:);
        R_Mod= R_Mod(1:size(circvar_Ref,1),:);
        Rd_Mod = Rd_Mod(1:size(circvar_Ref,1),:);
        Y_Mod = Y_Mod(1:size(circvar_Ref,1),:);
    else
        circvar_Mod = [circvar_Mod;zeros(size(circvar_Ref,1)-size(circvar_Mod,1),size(circvar_Mod,2))];
        R_Mod = [R_Mod;zeros(size(R_Ref,1)-size(R_Mod,1),size(circvar_Mod,2))];
        Rd_Mod = [Rd_Mod;zeros(size(Rd_Ref,1)-size(Rd_Mod,1),size(circvar_Mod,2))];
        Y_Mod = [Y_Mod;zeros(size(Y_Ref,1)-size(Y_Mod,1),size(circvar_Mod,2))];
    end
end

if param.weightmode == 0
    weight = abs(Y_Mod).^2;
elseif param.weightmode == 1
    weight = abs(Y_Ref).^2;
elseif param.weightmode == -1
    weight = 1;
end
% FrameNo = size(R_Mod,2);
% FreqNo = size(R_Mod,1);
if ~param.trend
    % UnRMSE = db(mean(sqrt(sum(weight.*abs(R_Mod-R_Ref).^2)./sum(weight))));
    UnRMSE = ((sqrt(sum(weight.*abs(R_Mod-R_Ref).^2)./sum(weight))));
    % % 	UnPWSNR = mean(sum(weight)./sum((weight).*abs(R_Mod-R_Ref).^2));
elseif param.trend
    % UnRMSE = db(mean(sqrt(sum(weight.*abs(Rd_Mod-Rd_Ref).^2)./sum(weight))));
    UnRMSE = ((sqrt(sum(weight.*abs(Rd_Mod-Rd_Ref).^2)./sum(weight))));
    % % 	UnPWSNR = mean(sum(weight)./sum((weight).*abs(Rd_Mod-Rd_Ref).^2));
end

if UnRMSE==0
    UnRMSE = 0;
else
    UnRMSE_VUV=10*log10(UnRMSE)';
    if length(pv)>length(UnRMSE_VUV)
        pv=pv(1:length(UnRMSE_VUV));
    end
    UnRMSE=mean(UnRMSE_VUV(~isinf(UnRMSE_VUV) & pv>pvthreshold));
end

% HPSNR in unwrapped domain
if param.weightmode == 0
    num=sum(abs(Y_Mod).*abs(Y_Mod));
    if ~param.trend
        den=sum(abs(Y_Mod).*abs(Y_Mod).*(1-cos(R_Mod-R_Ref)));
    elseif param.trend
        den=sum(abs(Y_Mod).*abs(Y_Mod).*(1-cos(Rd_Mod-Rd_Ref)));
    end
elseif param.weightmode == 1
    num=sum(abs(Y_Ref).*abs(Y_Ref));
    if ~param.trend
        den=sum(abs(Y_Ref).*abs(Y_Ref).*(1-cos(R_Mod-R_Ref)));
    elseif param.trend
        den=sum(abs(Y_Ref).*abs(Y_Ref).*(1-cos(Rd_Mod-Rd_Ref)));
    end
end
if den ==0
    UnHPSNR = Inf;
else
    HPSNR_VUV=10*log10(num./den)';
    if length(pv)>length(HPSNR_VUV)
        pv=pv(1:length(HPSNR_VUV));
    end
    UnHPSNR=mean(HPSNR_VUV(~isinf(HPSNR_VUV) & pv>pvthreshold));
end
end