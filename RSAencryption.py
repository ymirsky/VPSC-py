from numpy import *

# This class impliments the RSA physical signal cipher described in:
# Khalil, M. I. "Real-Time Encryption/Decryption of Audio Signal." International Journal of Computer Network and Information Security 8.2 (2016): 25.
# Warning: this implimentation is not secure since it uses small primes which are hard-coded. This implimentation is to be used for demonstration purposes only.

class RSAencryption:
    def __init__(self, frame_size, minValue=5):
        self.N = frame_size
        self.e = 17#3#17.0 # HARD CODED KEYS FOR TESTIGN ONLY
        self.d = 413#27#413.0 # HARD CODED KEYS FOR TESTIGN ONLY
        self.n = 3233#55#3233.0 # HARD CODED KEYS FOR TESTIGN ONLY
        self.minValue = minValue

    def encrypt(self, frame, keySeed=None):
        #F = ((frame+self.minValue)*1000).astype(int64)
        out = zeros(len(frame),dtype=double)
        for i in range(len(frame)):
            out[i] = int((frame[i]+self.minValue)*300)**self.e % self.n
        return out
        #return mod(round_(10 * (1 + frame)) ** self.e, self.n)

    def decrypt(self, frame, keySeed=None):
        #F = frame.astype(int64)
        out = zeros(len(frame),dtype=double)
        for i in range(len(frame)):
            out[i] = double(int(frame[i])**self.d % self.n)/300 - self.minValue
        return out#(F**self.d % self.n).astype(double)/1000 - self.minValue
        #return (mod(frame ** self.d, self.n)) / 10 - 1
