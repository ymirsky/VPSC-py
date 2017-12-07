from numpy import *
from numpy.fft import *

class VPSC:
    def __init__(self,frame_size,max_freq_energy,expected_noise_energy=0,targetBand=None):
        self.N = frame_size
        self.Phi = max_freq_energy
        self.preBuf = sqrt(expected_noise_energy)
        self.curKey = None
        self.targetBand = targetBand # the bins to encrypt

    def encrypt(self,frame, keySeed=None):
        self.curKey = genKey(self.N,self.Phi+2*self.preBuf,keySeed)
        return enc(frame, self.curKey, self.Phi, preBuf=self.preBuf, binRange=self.targetBand)

    def decrypt(self,frame, keySeed=None):
        if keySeed is not None:
            self.curKey = genKey(self.N, self.Phi+2*self.preBuf, keySeed)
        return dec(frame, self.curKey, self.Phi, preBuf=self.preBuf, binRange=self.targetBand)

def pol2complex(rho, phi):
    x = rho * cos(phi)
    y = rho * sin(phi)
    return x+1j*y

def enc(M,K, Phi, preBuf = 0, binRange=None):
    if binRange is None:
        binRange = arange(0, int(ceil(len(M) / 2)) + 1)  # the bins to encrypt
    Phi *= len(M)/2
    preBuf *= len(M)/2

    #transform to freq plane
    Fm = rfft(M)

    #decontruct: magnitudes and phases
    Mm = abs(Fm)
    Ma = angle(Fm)

    #add key
    Cm = Mm
    Cm[binRange] += preBuf #boost origional singal
    Cm[binRange] = (Cm[binRange] + K[0][binRange]) % (Phi + 2*preBuf) #encrypt boosted signal
    Cm[binRange] += preBuf #add carrier energy
    Ca = Ma
    Ca[binRange] = Ma[binRange] + K[1][binRange]

    #recontruct:
    Fc = pol2complex(Cm,Ca)
    return irfft(Fc)

def dec(C,K, Phi, preBuf=0, binRange = None):
    if binRange is None:
        binRange = arange(0, int(ceil(len(C) / 2)) + 1)  # the bins to encrypt
    Phi *= len(C)/2
    preBuf *= len(C)/2

    #transform to freq plane
    Fc = rfft(C)

    #decontruct: magnitudes and phases
    Cm = abs(Fc)
    Ca = angle(Fc)

    #remove impossible magnitudes
    CmR = Cm[binRange]
    CmR[CmR>(Phi + 3*preBuf)] = Phi + 3*preBuf
    Cm[binRange] = CmR

    #remove carrier energy
    Cm[binRange] -= preBuf
    Cm[Cm<0]=0

    #decrypt
    Mm = Cm
    Mm[binRange] = ((Cm[binRange] - K[0][binRange]) % (Phi + 2 * preBuf))

    #remove energy boost
    Mm[binRange] -= preBuf
    Mm[Mm<0]=0

    Mm[-1] = Mm[-2]

    #decrypt phase
    Ma = Ca
    Ma[binRange] = Ca[binRange] - K[1][binRange]

    #recontruct:
    Fm = pol2complex(Mm,Ma)
    return irfft(Fm)

def genKey(n, Phi,seed=None):
    if seed is not None:
        random.seed(seed)
    Phi *= n/2
    keylen = int(ceil(n / 2))+1
    Km = random.rand(keylen) * Phi
    Ka = zeros(keylen)
    Ka[1:] = random.rand(keylen-1) * 2 * pi
    return (Km,Ka)

