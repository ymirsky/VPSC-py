from numpy import *
from numpy.fft import *
from VPSC import *
from FFTscrambler import *
from LogEncryption import *
from RSAencryption import *
import warnings
warnings.filterwarnings('ignore') #this is here because the ALM and RSA encryption methods have many overflow issues

#bits is a vector of four '0' or '1' values
#fc is carrier frequency. eg: fc = 19.53125;
#Ts is the sampling interval. eg: Ts = 1 / 1000
def mod_QAM16_frame(bits,fc,Ts,frame_size,method=1):
    #calculate timestamps of samples
    #t = arange(0,Ts*(frame_size - 1),Ts)
    t = arange(0, Ts * (frame_size), Ts)

    #determine qam16 cart x-y values
    if bits[0]==0 and bits[1]==0:
        x = -3
    elif bits[0]==0 and bits[1]==1:
        x = -1
    elif bits[0]==1 and bits[1]==1:
        x = 1
    else:
        x = 3

    if bits[2]==0 and bits[3]==0:
        y = -3
    elif bits[2]==0 and bits[3]==1:
        y = -1
    elif bits[2]==1 and bits[3]==1:
        y = 1
    else:
        y = 3

    symbol = x + 1j*y

    if method == 1:
        F = zeros(int(ceil(frame_size/2))+1,dtype=complex)
        fs = 1/Ts
        carrier_indx = int(ceil(fc / (fs / frame_size)))
        F[carrier_indx] = symbol*(frame_size/2)
        frame = irfft(F)
    else:
    #frame = sin(2 * pi * fc * t + angle(symbol) + pi * .77) * abs(symbol)
        frame = sin(2 * pi * fc * t + angle(symbol)) * abs(symbol)
    return frame


# fc is carrier frequency. eg: fc = 19.53125;
# Ts is the sampling interval. eg: Ts = 1 / 1000
def demod_QAM16_frame(frame,fc, Ts):
    symbols = array([(-3, -3),(-3,-1),(-3,1),(-3,3),
               (-1,-3),(-1,-1),(-1,1),(-1,3),
               (1,-3),(1,-1),(1,1),(1,3),
               (3,-3),(3,-1),(3,1),(3,3)])

    values = array([(0,0,0,0),(0,0,0,1),(0,0,1,1),(0,0,1,0),
              (0,1,0,0),(0,1,0,1),(0,1,1,1),(0,1,1,0),
              (1,1,0,0),(1,1,0,1),(1,1,1,1),(1,1,1,0),
              (1,0,0,0),(1,0,0,1),(1,0,1,1),(1,0,1,0)])

    fs = 1/Ts
    carrier_indx = int(ceil(fc / (fs / len(frame))))
    symb = rfft(frame)[carrier_indx]/( len(frame)/2)
    #symb = abs(symb)*exp(1j*(angle(symb)+pi/2))
    dist2symbols = sqrt((real(symb) - symbols[:, 0] )**2 + (imag(symb) - symbols[:,1])**2)
    bits = values[argmin(dist2symbols), :]
    return bits

def ber(binvec1,binvec2):
    Errors = sum(logical_xor(binvec1,binvec2))
    BER = Errors/len(binvec1)
    return (BER,Errors)

def MPP(samples, delays, gains, Ts, relative_velocity = 0): #relative velocity is in m/s
    n_paths = len(delays)
    n_samples = len(samples)
    maxDelay = max(delays)
    binDelay = int(floor(maxDelay/Ts))
    #normalize the gains
    amps = 10**(array(gains)/10)
    #amps_norm = amps/sum(amps)
    amps_norm = amps / amps[0] #normalized relative to LOS

    M = zeros((n_samples + binDelay, n_paths))
    for i in range(n_paths):
        bindel = int(floor(delays[i]/Ts))
        M[bindel:(n_samples+bindel),i] = samples*amps_norm[i]#*sqrt(amplification)
    return sum(M,axis=1)[0:n_samples]

def doppler_shift_samples(samples,frame_size,fs,relative_velocity):#X is the rfft of the signal
    n_frames = int(floor(len(samples)/frame_size))
    for i in range(n_frames):
        frame = samples[i * frame_size:(i + 1) * frame_size]
        X = rfft(frame)
        freq_bins = arange(0,fs/2,(fs/2)/len(X))
        d_shift = 1 + relative_velocity/3e8  #relative velocity is in m/s
        freq_bins_shifted = freq_bins*d_shift
        samples[i * frame_size:(i + 1) * frame_size] = irfft(interp(freq_bins_shifted, freq_bins, X))
    return samples

def sim(fc,fs,N,n_frames,SNR,delays,gains,relative_velocity):
    #scenario config
    Ts = 1/fs#1 / 1000 #sample interval [sec]

    #sim config
    SNRlin = 10**(SNR/10)
    encBand = array([fc - .7 * 1000 * 1000, fc + .7 * 1000 * 1000])  # 0.7 Mhz from either side of the carrier freq
    freqBins = arange(int(ceil(encBand[0] / ((1 / Ts) / N))),
                      int(ceil(encBand[1] / ((1 / Ts) / N))))  # arange(0,int(ceil(N/2+1))) #array([126,127,128])
    Es = sqrt(3 ** 2 + 3 ** 2) * len(freqBins)  # mean(mod_QAM16_frame(array([1,0,1,0]),fc,Ts,N)**2)
    noise_eng = sqrt(Es/SNRlin)
    relative_velocity = relative_velocity*1000 #m/s  :I multiply by 1.4 m/s is walking speed, I multiply by 1000 to give the relative doppler affect as if we are working in 1900Mhz as ooposed to 1.9Mhz

    ####### INIT ############
    #freqBins = array([int(ceil(fc/ ((1/Ts) / N)))])
    cipher_vpsc = VPSC(N, 5, (noise_eng)**(1/5), freqBins)
    cipher_scram = FScrambler(N,freqBins)
    cipher_log = LogEncryption(N,5,5)
    cipher_rsa = RSAencryption(N,5)
    Data_TX = random.randint(0,2,(n_frames,4))
    Data_RX_VPSC = zeros((n_frames,4))
    Data_RX_Scrambler = zeros((n_frames,4))
    Data_RX_LOG = zeros((n_frames,4))
    Data_RX_RSA = zeros((n_frames,4))
    Data_RX_Baseline = zeros((n_frames,4))

    S = zeros(n_frames*N)
    ####### Gen Timeseries ###########
    print("Generating Modulated Signal")
    for i in range(n_frames):
        data_tx = Data_TX[i]
        # Modulate QAM16 frame
        frame_m_tx = mod_QAM16_frame(data_tx, fc, Ts, N)
        S[i*N:(i+1)*N] = frame_m_tx


    S_baseline = zeros(n_frames*N)
    S_vpsc = zeros(n_frames*N)
    S_scram = zeros(n_frames*N)
    S_log = zeros(n_frames*N)
    S_rsa = zeros(n_frames*N)

    ####### SEND ###########
    print("Encrypting Signal... VPSC, FCS, ALM, and RSA")
    for i in range(n_frames):
        if i % 1000 == 0:
            print("Symbols Processed: "+str(i))
        data_tx = Data_TX[i]

        # Modulate QAM16 frame
        frame_m_tx = mod_QAM16_frame(data_tx,fc,Ts,N)

        #Encrypt with VPSC
        #print("Encrypting with VPSC")
        frame_c_tx_vpsc = cipher_vpsc.encrypt(frame_m_tx, i)
        #Encrypt with Scrambler
        #print("Encrypting with FCS")
        frame_c_tx_scram = cipher_scram.encrypt(frame_m_tx, i)
        #Encrypt with log
        #print("Encrypting with ALM")
        frame_c_tx_log = cipher_log.encrypt(frame_m_tx, i)
        #Encrypt with RSA
        #print("Encrypting with RSA")
        frame_c_tx_rsa = cipher_rsa.encrypt(frame_m_tx)

        #Pass through channel
        S_baseline[i * N:(i + 1) * N] = frame_m_tx
        S_vpsc[i * N:(i + 1) * N] = frame_c_tx_vpsc
        S_scram[i * N:(i + 1) * N] = frame_c_tx_scram
        S_log[i * N:(i + 1) * N] = frame_c_tx_log
        S_rsa[i * N:(i + 1) * N] = frame_c_tx_rsa


    # S_baseline_mpp = S_baseline
    # S_vpsc_mpp = S_vpsc
    # S_scram_mpp = S_scram
    # S_log_mpp = S_log
    # S_rsa_mpp = S_rsa
    ####### Process MPP ###########
    print("Processing MPP")
    S_baseline_mpp = MPP(S_baseline,delays,gains,Ts)
    S_vpsc_mpp = MPP(S_vpsc,delays,gains,Ts)
    S_scram_mpp = MPP(S_scram,delays,gains,Ts)
    S_log_mpp = MPP(S_log,delays,gains,Ts)
    S_rsa_mpp = MPP(S_rsa,delays,gains,Ts)

    ####### Process Doppler ###########
    print("Processing Doppler")
    if relative_velocity != 0:
        S_baseline_mpp = doppler_shift_samples(S_baseline_mpp,N,fs,relative_velocity)
        S_vpsc_mpp = doppler_shift_samples(S_vpsc_mpp,N,fs,relative_velocity)
        S_scram_mpp = doppler_shift_samples(S_scram_mpp,N,fs,relative_velocity)
        S_log_mpp = doppler_shift_samples(S_log_mpp, N, fs, relative_velocity)
        S_rsa_mpp = doppler_shift_samples(S_rsa_mpp, N, fs, relative_velocity)

    ####### Process AGWN ###########
    print("Processing AWGN")
    awgn = random.randn(n_frames*N)*noise_eng
    S_baseline_mpp += awgn
    S_vpsc_mpp += awgn
    S_scram_mpp += awgn
    S_log_mpp += awgn
    S_rsa_mpp += awgn

    ####### Receive ###########
    print("Decrypting Signal...  VPSC, FCS, ALM, and RSA")
    S_vpsc_decrypted = zeros(n_frames*N)
    S_scram_decrypted = zeros(n_frames*N)
    S_log_decrypted = zeros(n_frames*N)
    S_rsa_decrypted = zeros(n_frames*N)
    S_vpsc_decrypted_noMPP = zeros(n_frames*N)
    for i in range(n_frames):
        if i % 1000 == 0:
            print("Symbols Processed: "+str(i))
        frame_m_rx_baseline = S_baseline_mpp[i * N:(i + 1) * N]
        frame_c_rx_vpsc = S_vpsc_mpp[i * N:(i + 1) * N]
        frame_c_rx_scram = S_scram_mpp[i * N:(i + 1) * N]
        frame_c_rx_log = S_log_mpp[i * N:(i + 1) * N]
        frame_c_rx_rsa = S_rsa_mpp[i * N:(i + 1) * N]
        frame_c_rx_vpsc_noMPP = S_vpsc[i * N:(i + 1) * N]

        #Decrypt with VPSC
        #print("Decrypting with VPSC")
        frame_m_rx_vpsc_noMPP = cipher_vpsc.decrypt(frame_c_rx_vpsc_noMPP, i)
        frame_m_rx_vpsc = cipher_vpsc.decrypt(frame_c_rx_vpsc, i)
        S_vpsc_decrypted[i * N:(i + 1) * N] = frame_m_rx_vpsc
        S_vpsc_decrypted_noMPP[i * N:(i + 1) * N] = frame_m_rx_vpsc_noMPP

        #Decrypt with Scramber
        #print("Decrypting with FCS")
        frame_m_rx_scram = cipher_scram.decrypt(frame_c_rx_scram, i)
        S_scram_decrypted[i * N:(i + 1) * N] = frame_m_rx_scram

        #Decrypt with Log
        #print("Decrypting with ALM")
        frame_m_rx_log = cipher_log.decrypt(frame_c_rx_log, i)
        S_log_decrypted[i * N:(i + 1) * N] = frame_m_rx_log

        #Decrypt with RSA
        #print("Decrypting with RSA")
        frame_m_rx_rsa = cipher_rsa.decrypt(frame_c_rx_rsa, i)
        S_rsa_decrypted[i * N:(i + 1) * N] = frame_m_rx_rsa

        # Demodulate QAM16 frame
        #print("Demodulating Singal")
        data_rx_vpsc = demod_QAM16_frame(frame_m_rx_vpsc,fc,Ts)
        data_rx_scram = demod_QAM16_frame(frame_m_rx_scram,fc,Ts)
        data_rx_log = demod_QAM16_frame(frame_m_rx_log,fc,Ts)
        data_rx_rsa = demod_QAM16_frame(frame_m_rx_rsa,fc,Ts)
        data_rx_baseline = demod_QAM16_frame(frame_m_rx_baseline,fc,Ts)

        Data_RX_VPSC[i,:] = data_rx_vpsc
        Data_RX_Baseline[i, :] = data_rx_baseline
        Data_RX_Scrambler[i, :] = data_rx_scram
        Data_RX_RSA[i, :] = data_rx_rsa
        Data_RX_LOG[i, :] = data_rx_log

    ####### END ###########
    # collect results

    print("Collecting Results...")
    BER_VPSC, Err_VPSC = ber(Data_TX.ravel(),Data_RX_VPSC.ravel())
    ErrL_VPSC = zeros(n_frames,dtype=bool)
    for i in range(n_frames):
        if sum(logical_xor(Data_TX[i,:],Data_RX_VPSC[i,:]))!=0:
            ErrL_VPSC[i] = True
    BER_Baseline, Err_Baseline = ber(Data_TX.ravel(),Data_RX_Baseline.ravel())
    ErrL_Baseline = zeros(n_frames,dtype=bool)
    for i in range(n_frames):
        if sum(logical_xor(Data_TX[i,:],Data_RX_Baseline[i,:]))!=0:
            ErrL_Baseline[i] = True
    BER_Scrambler, Err_Scrambler = ber(Data_TX.ravel(),Data_RX_Scrambler.ravel())
    ErrL_Scrambler = zeros(n_frames,dtype=bool)
    for i in range(n_frames):
        if sum(logical_xor(Data_TX[i,:],Data_RX_Scrambler[i,:]))!=0:
            ErrL_Scrambler[i] = True
    BER_LOG, Err_LOG = ber(Data_TX.ravel(),Data_RX_LOG.ravel())
    ErrL_LOG = zeros(n_frames,dtype=bool)
    for i in range(n_frames):
        if sum(logical_xor(Data_TX[i,:],Data_RX_LOG[i,:]))!=0:
            ErrL_LOG[i] = True
    BER_RSA, Err_RSA = ber(Data_TX.ravel(),Data_RX_RSA.ravel())
    ErrL_RSA = zeros(n_frames,dtype=bool)
    for i in range(n_frames):
        if sum(logical_xor(Data_TX[i,:],Data_RX_RSA[i,:]))!=0:
            ErrL_RSA[i] = True

    print("No Enc: the BER was: "+str(BER_Baseline)+ " with "+str(Err_Baseline)+" errors.")
    print("Freq Scramber: the BER was: "+str(BER_Scrambler)+ " with "+str(Err_Scrambler)+" errors.")
    print("Log Enc: the BER was: "+str(BER_LOG)+ " with "+str(Err_LOG)+" errors.")
    print("RSA Enc: the BER was: "+str(BER_RSA)+ " with "+str(Err_RSA)+" errors.")
    print("VPSC: the BER was: "+str(BER_VPSC)+ " with "+str(Err_VPSC)+" errors.")


    #collect symbols:
    symbs_o_mpp = zeros(n_frames, dtype=complex)
    symbs_v_mpp = zeros(n_frames, dtype=complex)
    symbs_v_NOmpp = zeros(n_frames, dtype=complex)
    symbs_b_mpp = zeros(n_frames, dtype=complex)
    symbs_scr_mpp = zeros(n_frames, dtype=complex)
    symbs_log_mpp = zeros(n_frames, dtype=complex)
    symbs_rsa_mpp = zeros(n_frames, dtype=complex)
    for i in range(n_frames):
        symbs_o_mpp[i] = rfft(S_baseline[i * N:(i + 1) * N])[int(ceil(fc / ((1 / Ts) / N)))] / (N / 2)
        symbs_v_mpp[i] = rfft(S_vpsc_decrypted[i * N:(i + 1) * N])[int(ceil(fc / ((1 / Ts) / N)))] / (N / 2)
        symbs_v_NOmpp[i] = rfft(S_vpsc_decrypted_noMPP[i * N:(i + 1) * N])[int(ceil(fc / ((1 / Ts) / N)))] / (N / 2)
        symbs_b_mpp[i] = rfft(S_baseline_mpp[i * N:(i + 1) * N])[int(ceil(fc / ((1 / Ts) / N)))] / (N / 2)
        symbs_scr_mpp[i] = rfft(S_scram_decrypted[i * N:(i + 1) * N])[int(ceil(fc / ((1 / Ts) / N)))] / (N / 2)
        symbs_log_mpp[i] = rfft(S_log_decrypted[i * N:(i + 1) * N])[int(ceil(fc / ((1 / Ts) / N)))] / (N / 2)
        symbs_rsa_mpp[i] = rfft(S_rsa_decrypted[i * N:(i + 1) * N])[int(ceil(fc / ((1 / Ts) / N)))] / (N / 2)

    return (BER_Baseline,BER_Scrambler,BER_LOG,BER_RSA,BER_VPSC,symbs_o_mpp,symbs_b_mpp,symbs_scr_mpp,symbs_log_mpp,symbs_rsa_mpp,symbs_v_mpp,ErrL_Baseline,ErrL_Scrambler,ErrL_RSA,ErrL_LOG,ErrL_VPSC)

def calc_SNR_from_EbN0(Tsymb,Ts,target_ebn0):
    return 10 * log10(0.5 * 0.0000667 / Ts) - target_ebn0







