The matlab files enclosed in this toolbox can be used to tabulate gain functions 
for log-spectral magnitude MMSE estimators under an assumed Generalized- 
Gamma model for the clean speech magnitude DFT coefficients.
 
For the theory behind these estimators and constraints on the parameters we refer to the article
 
[1] R.C.Hendriks, R.Heusdens and J.Jensen 
     "Log-spectral magnitude MMSE estimators under super-Gaussian densities", Interspeech, 2009.
 


 
Short description of the 2 main m-files (see the headers of the files for more info):
 
For an assumed Generalized-Gamma prior density of the magnitude DFT coefficients with gamma=2 and specific 
nu parameter the m-file [G1]=TabulateGainGamma2logmmse(Rprior,Rpost,nu)
tabulates the gain function  for the  log-spectral magnitude DFT coefficients, 
 For mathematical expressions 
of the gain functions for these estimators see [1].
The range of a priori and a posteriori SNRs is
-40 to 40 dB in 1 dB steps. Each row of the gain matrices is for a different a priori SNR, while a posteriori
SNR varies along columns.
 
Given the tabulated gain function, a vector of gain values for pairs of a priori and a posteriori SNRs can be
selected using the m-file
[gains]=lookup_gain_in_table(G,a_post,a_priori,a_post_range,a_priori_range);
where a_post and a_priori are vectors with the a posteriori and a priori SNRs respectively.
The vectors a_post and a_priori should have equal lengths. The parameters a_post_range and a_priori_range indicate
the ranges in dBs used in the gain table G.

The scriptfile log_mmse_supergaus gives an example of the usage of the aforementioned m-files.
 
Implementations of the special functions are based on 
S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996) with implementations available 
online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
 
The implementations of these special functions in the toolbox have been adapted with respect to the original implementations 
such that they can handle vector arguments as well.
 
Copyright 2009: Delft University of Technology, Information and
Communication Theory Group. The software is free for non-commercial use.
This program comes WITHOUT ANY WARRANTY.
 
June, 2009
R. C. Hendriks

 