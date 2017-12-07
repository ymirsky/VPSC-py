from numpy import *
from simulator import *
from LogEncryption import *

# Simulation parameters
fc = 1.9*1000*1000#carrier freq [kHz]
fs = 8.08*fc
N=1024
n_frames = 10000  # number of frames to simulate
delays = array([0, 0.5e-06, 1.0e-06, 0.2e-06])  # array([0,0.2e-06,0.5e-06])# #MPP delays in sec
gains = array([0, -9, -12, -6])  # array([0,-6,-9])# #MPP signal gains in DB
SNRs = [300] #arange(0,25) #calc_SNR_from_EbN0(0.0000667,1/fs,EbN0s)  # db
delay_factor = [0]#arange(0,5,.1)
relative_velocity = 1.4#m/s

print("Initializing Simulation")
BERs = zeros((len(delay_factor),5))
for i in range(len(delay_factor)):
    BER_Baseline, BER_Scrambler, BER_LOG, BER_RSA, BER_VPSC, symbs_o_mpp, symbs_b_mpp, symbs_scr_mpp, symbs_log_mpp, symbs_rsa_mpp, symbs_v_mpp, ErrL_Baseline, ErrL_Scrambler, ErrL_RSA, ErrL_LOG, ErrL_VPSC = sim(fc,fs,N,n_frames,SNRs[i],delays,gains,relative_velocity)
    BERs[i,0]=BER_Baseline
    BERs[i,1]=BER_Scrambler
    BERs[i,2]=BER_LOG
    BERs[i,3]=BER_VPSC
    BERs[i,4]=BER_RSA
print("Simulation Complete")

# ######## PLOTS ############
from matplotlib.pyplot import *
f, axarr = subplots(ncols=4, sharex=True,sharey=True)
f.text(0.5, 0.04, 'In-phase Amplitude', ha='center')
f.text(0.04, 0.5, 'Quadrature Amplitude', va='center', rotation='vertical')
axarr[0].scatter(imag(symbs_b_mpp[logical_not(ErrL_Baseline)]),real(symbs_b_mpp[logical_not(ErrL_Baseline)]),alpha=.2,color="gray")
axarr[0].scatter(imag(symbs_b_mpp[ErrL_Baseline]),real(symbs_b_mpp[ErrL_Baseline]),alpha=.2,color="red")
axarr[0].scatter(imag(symbs_o_mpp),real(symbs_o_mpp),alpha=.2,color="black",marker="+",s=150)
axarr[0].set_title('No Encryption. BER: '+str(BER_Baseline))
axarr[1].scatter(imag(symbs_scr_mpp[logical_not(ErrL_Scrambler)]),real(symbs_scr_mpp[logical_not(ErrL_Scrambler)]),alpha=.2,color="gray")
axarr[1].scatter(imag(symbs_scr_mpp[ErrL_Scrambler]),real(symbs_scr_mpp[ErrL_Scrambler]),alpha=.2,color="red")
axarr[1].scatter(imag(symbs_o_mpp),real(symbs_o_mpp),alpha=.2,color="black",marker="+",s=150)
axarr[1].set_title('ALM. BER: '+str(BER_Scrambler))
axarr[2].scatter(imag(symbs_log_mpp[logical_not(ErrL_LOG)]),real(symbs_log_mpp[logical_not(ErrL_LOG)]),alpha=.2,color="gray")
axarr[2].scatter(imag(symbs_log_mpp[ErrL_LOG]),real(symbs_log_mpp[ErrL_LOG]),alpha=.2,color="red")
axarr[2].scatter(imag(symbs_o_mpp),real(symbs_o_mpp),alpha=.2,color="black",marker="+",s=150)
axarr[2].set_title('FCS. BER: '+str(BER_Scrambler))
axarr[3].scatter(imag(symbs_v_mpp[logical_not(ErrL_VPSC)]),real(symbs_v_mpp[logical_not(ErrL_VPSC)]),alpha=.2,color="gray")
axarr[3].scatter(imag(symbs_v_mpp[ErrL_VPSC]),real(symbs_v_mpp[ErrL_VPSC]),alpha=.2,color="red")
axarr[3].scatter(imag(symbs_o_mpp),real(symbs_o_mpp),alpha=.2,color="black",marker="+",s=150)
axarr[3].set_title('VPSC. BER: '+str(BER_VPSC))
axarr[0].set_xlim([-4, 4])
axarr[0].set_ylim([-4, 4])
show()