from numpy import *
from numpy.fft import *

# This class impliments the log-based physical signal cipher described in:
# Khalil, M. I. "Real-Time Encryption/Decryption of Audio Signal." International Journal of Computer Network and Information Security 8.2 (2016): 25.

class LogEncryption:
    def __init__(self, frame_size, max_keyval, raising = 10):
        self.N = frame_size
        self.max_keyval = max_keyval
        self.raising = raising
        self.curKey = None

    def encrypt(self, frame, keySeed=None):
        self.curKey = self.genKey(keySeed)
        raisedM = frame + self.raising
        return (log10(self.curKey [0] * raisedM) / log10(self.curKey [1]))

    def decrypt(self, frame, keySeed=None):
        if keySeed is not None:
            self.curKey = self.genKey(keySeed)
        # decrypt with logA and logB keys
        decM = (self.curKey[1] ** frame) / self.curKey[0]
        return decM - self.raising

    def genKey(self, seed=None):
        if seed is not None:
            random.seed(seed)
        #max_value *= n / 2
        keylen = self.N#int(ceil(n / 2)) + 1
        Ka = random.rand(keylen) * self.max_keyval + 0.000001
        Kb = random.rand(keylen) * self.max_keyval + 0.000001
        return (Ka, Kb)
