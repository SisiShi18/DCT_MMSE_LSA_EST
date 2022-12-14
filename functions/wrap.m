function phase = wrap(phase)

phase = phase - round(phase/2/pi)*2*pi;

if phase>pi;        phase=phase-2*pi;
elseif phase<-pi;   phase=phase+2*pi;
end
return
end