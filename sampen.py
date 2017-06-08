
import os
import ctypes
import warnings
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt

def init_libsampen():
    libdir = os.path.dirname(__file__)
    libfile = os.path.join(libdir, "libsampen.so")
    lib = ctypes.CDLL(libfile)

    roptr = ndpointer(float, flags=('C','A'))
    sizeptr = ndpointer(ctypes.c_size_t, flags=('C','A'))
    rwptr = ndpointer(float, flags=('C','A','W'))
    charptr = ndpointer(ctypes.c_char, flags=('C','A'))
    intptr = ndpointer(ctypes.c_int, flags=('C','A'))

    lib.sampen.restype = None
    lib.sampen.argtypes = [rwptr, ctypes.c_int, ctypes.c_double, ctypes.c_ulong, rwptr]
    lib.sampen2.restype = None 
    lib.sampen2.argtypes = [rwptr, ctypes.c_int, ctypes.c_double, ctypes.c_ulong, rwptr, rwptr]
    lib.normalize.restype = None
    lib.normalize.argtypes = [rwptr, ctypes.c_ulong]
    lib.readdata.restype = roptr
    lib.readdata.argtypes = [charptr, intptr]
    lib.SampleEntropy.restype = ctypes.c_double
    lib.SampleEntropy.argtypes = [rwptr, ctypes.c_ulong, ctypes.c_int, ctypes.c_double, ctypes.c_int]
    lib.MSE.restype = None
    lib.MSE.argtypes = [rwptr, ctypes.c_ulong, ctypes.c_int, ctypes.c_double, ctypes.c_int, rwptr]
    lib.MSESD.restype = None
    lib.MSESD.argtypes = [rwptr, ctypes.c_ulong, ctypes.c_int, ctypes.c_double, ctypes.c_int, rwptr, rwptr]
    lib.CoarseGraining.restype = None
    lib.CoarseGraining.argtypes = [rwptr, ctypes.c_ulong, ctypes.c_int, roptr]

    return lib;

libsampen = init_libsampen()

def sampen(data, m=2, r=0.2, norm=True, sd=True):
    N = len(data)
    data = np.array(data)##make sure we are dealing with NumPy array
    if data.ndim != 1:
        raise ValueError("Data must be 1D array")
    data = np.require(data, float, ('C', 'A'))
    
    ## Output arrays
    SampEn = np.empty(m+1, dtype=float)
    SampEn = np.require(SampEn, float, ('C', 'A'))
    if sd:
        SampEnSD = np.empty(m+1, dtype=float)
        SampEnSD = np.require(SampEnSD, float, ('C', 'A'))

    ## Call libsampen.so routine
    if norm:
        libsampen.normalize(data, N)
    if sd:
        libsampen.sampen2(data, m, r, N, SampEn, SampEnSD)
    else:# faster code without estimating SD
        libsampen.sampen(data, m, r, N, SampEn)

    if sd:
        return (SampEn, SampEnSD)
    else:
        return SampEn

def mse(data, m=2, r=0.2, maxscale=10, sd=False):
    N = len(data)
    data = np.array(data)## make sure we deal with NumPy array
    if data.ndim != 1:
        raise ValueError("Data must be 1D array")
    data = np.require(data, float, ('C', 'A'))
    SampEn = np.empty(maxscale, dtype=float)
    SampEn = np.require(SampEn, float, ('C', 'A'))
    if sd:
        SampEnSD = np.empty(maxscale, dtype=float)
        SampEnSD = np.require(SampEnSD, float, ('C', 'A'))

    ## Call C routine
    if sd:
        libsampen.MSESD(data, N, m, r, maxscale, SampEn, SampEnSD)
        return (SampEn, SampEnSD);
    else:
        libsampen.MSE(data, N, m, r, maxscale, SampEn)
        return SampEn;
    

# np.random.seed(10)
# d = np.random.randn(10**5)
# d = np.require(d, float, ('C', 'A'))
# d1 = np.loadtxt('../sampentest.txt')
# d2 = np.loadtxt('../test.txt')

# SE = mse(d)
# plt.plot(np.arange(10)+1,SE,'o')
# SE = mse(d1)
# plt.plot(SE,'x')
# SE = mse(d2)
# plt.plot(SE,'>')
# plt.show()

# print(libsampen.SampleEntropy(d, len(d), 2, 0.2, 1))
# print(sampen(d, 2, 0.2, norm=False))
# SE, SESD = mse(d, sd=True)
# plt.errorbar(np.arange(10)+1,SE, yerr=SESD)
# plt.show()

# print("SampEn(sd=False) =", sampen(d,sd=False))
# print("SampEn =", sampen(d))
# S2,Sd2 = sampen(d1)
# print("SampEn =", S2,"SampEnSD =",Sd2)

# SE = np.empty(1, dtype=float)
# libsampen.SampleEntropy(d1,len(d1),2,0.2,1,SE)
# print("SampleEntropy =",SE)

# SES = np.empty(10, dtype=float)
# libsampen.MSE(d,len(d),2,0.2,10,SES)
# plt.plot(SES,'o')

# libsampen.normalize(d, len(d))
# dc = np.empty(len(d), dtype=float)
# sampEn = np.empty(3, dtype=float)
# for i in range(10):
#     libsampen.CoarseGraining(dc,len(d),i+1,d)
#     libsampen.sampen(dc, 2, 0.2, len(d)//(i+1), sampEn)
#     ##SES[i] = libsampen.SampleEntropy(dc,len(d),2,0.2,i+1)
#     SES[i] = sampEn[-1]

# plt.plot(SES,'x')
# plt.show()
