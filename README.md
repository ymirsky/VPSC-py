# Overview
In this repository you will find a Python implementation of the Vernam Pysical Signal Cipher (VPSC); a method for encrypting waveform signals in a manner which preserves the bandwidth of the origional singal and can *theoretically* provide perfect secrecy (depending on the manner in which the symetric keys are created). From,

An Encrytpion System for Securing Physical Signals
By *Yisroel Mirsky, Benjamin Fedidat, and Yoram Haddad*
Published in the SecureComm 2020 Proceedings (citation information below)

# Abstract
Secure communication is a necessity. However, common practice is to apply encryption
only to the upper layers of the communication protocol stack. This type of
implementation exposes network information to eavesdroppers and malicious attackers,
as the lower layers travel across public channels.
Some of this exposed information includes the channel's data rate, channel type,
protocol type, frame header and routing information. A solution to this problem is
to encrypt the physical layer (the lowest layer) of the protocol stack, thereby
securing all subsequent layers. In order for this method to be practical in a modern
network, the encrypted signal must use the same amount of bandwidth as the original
signal, and the encryption process must be fast enough to avoid undesirable delays.
Furthermore, the method must also deal with the issue of noise mitigation and synchronization.

The Vernam Physical Signal cipher (VPSC) is a novel method for encrypting any analog signal. The VPSC can be used as a physical-layer signal-cipher, capable of encrypting any analog waveform signal while maintaining the properties mentioned above. The VPSC can theoretically provide perfect security, but in practicality, it is as secure as the cryptographic number generator which generates its keys. 
We also offer a method for noise mitigation which are essential for a practical implementation.
To the best of our knowledge, the VPSC is the only physical signal cipher which provides high secrecy, preserves bandwidth, and is robust against noise.

# The Code
In this repo you will find a class which impliments the VPSC, along with several other class files which impliment other physical signal ciphers (used as a baseline for performance). These implimentations are for demonstration purposes only.

The VPSC class's constructor has the following arguments:
* frame_size: the FFT frame sizes (e.g., 256, 512,...)
* max_freq_energy: the phi parameter (largest frequency magnitude expected in the input signal)
* expected_noise_energy=0: the noise energy (per frequency) expected to be added to the singal after transmission (before reception)
* targetBand=None: the indexes to the bins int the FFT which should be encrypted. If 'None' then the entire spectrum is encryted

To use an instace of the VPSC, use the encrypt and decrypt functions. You can specify the key (seed) manually, or alow it to increment automatically by default. 

The snippet below shows how to encrypt and decrypt a sine wave using the VPSC:
```
import numpy as np
from VPSC import *

# Generate sample sine wave:
N = 1024 #number of samples
t = np.arange(0,.01*N,.01)
x = np.sin(2 * pi * t)

#Set VPSC Parameters
frame_size = N # (number of samples per frame)
max_freq_energy = 1.2

#Init VPSC
cipher = VPSC(frame_size, max_freq_energy)

#Encrypt the signal x
c = cipher.encrypt(x)

#decrypt the signal x
x_m = cipher.decrypt(c)
```

# The Test Code
The repo also includes a simulator which tests each of the methods in encrypting an OFMDA wireless channel using a QAM-16 modualtion, in the presence of multipath propagation, doppler effect, and additivee white gaussain noise.

To run the simulation, simply run the script entitled 'testscript.py'

```
python testscript.py
```

# Citations
If you use this code in any way, please cite our paper:
```
@InProceedings{mirsky2020vpsc,
author="Yisroel Mirsky, Benjamin Fedidat, Yoram Haddad",
title="An Encryption System for Securing Physical Signals",
booktitle="Security and Privacy in Communication Networks (SecureComm)",
year="2020",
publisher="Springer International Publishing"}
```

Yisroel Mirsky
yisroel@post.bgu.ac.il

