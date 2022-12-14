clear all




% tabulate for all a posteriori and a priori SNRs of interest:

nu=0.2;
Rprior= -40:1:40; % in dB
Rpost=-40:1:40; % in dB
[table_logmmse]=TabulateGainGamma2logmmse(Rprior,Rpost,nu);



%%%% example usage
a_post_snr=10.^((-40:1:40)/10) ; % not in dB
a_priori_snr=10.^((-40:1:40)/10);% not in dB
%%extract the right gain-values from the gain_table.
[G_ampl_log]=lookup_gain_in_table_2(table_logmmse,a_post_snr,a_priori_snr',Rpost,Rprior);

%%%% example usage2
a_post_snr=10.^((-40:1:40)/10) ; % not in dB
a_priori_snr=10.^((5)/10);% not in dB
%%extract the right gain-values from the gain_table.
[G_ampl_log]=lookup_gain_in_table_2(table_logmmse,a_post_snr,a_priori_snr,Rpost,Rprior);


figure
plot(10*log10(a_post_snr),G_ampl_log)
%%%%