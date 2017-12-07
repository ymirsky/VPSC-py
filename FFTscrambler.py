from numpy import *
from numpy.fft import *

# This class impliments the frequency scrambling signal cipher described in:
# Matsunaga, Akira, Keiichiro Koga, and Michihisa Ohkawa. "An analog speech scrambling system using the FFT technique with high-level security." IEEE Journal on Selected Areas in Communications 7.4 (1989): 540-547.

class FScrambler:
    def __init__(self,frame_size,targetBand=None):
        self.N = frame_size
        self.curKey = None
        if targetBand is None:
            targetBand = arange(0,self.N/2).astype(int)
        self.binRange = targetBand # the bins to encrypt

    def encrypt(self,frame, keySeed=None):
        #init
        self.curKey = self.genKey(keySeed)
        #prep frame
        X = rfft(frame)
        X[self.binRange] = X[self.binRange][self.curKey]
        return irfft(X)

    def decrypt(self,frame, keySeed=None):
        if keySeed is not None:
            self.curKey = self.genKey(keySeed)
        #prep frame
        X = rfft(frame)
        Xm = zeros(len(self.binRange),dtype=complex)
        Xm[self.curKey] = X[self.binRange]
        X[self.binRange] = Xm
        return irfft(X)

    def genKey(self,seed=None):
        if seed is not None:
            random.seed(seed)
        keylen = len(self.binRange)
        K = random.permutation(keylen)
        return K

