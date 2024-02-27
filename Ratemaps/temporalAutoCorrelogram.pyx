from cython.view cimport array as cvarray
import numpy as np
cimport numpy as np

cdef extern from "temporalAutoCorr.c":
    void  binnedISI(double *spikeTs, long *ISICounts, int spikeTsLength)
    
def getBinnedISI(np.ndarray spikeTimestamps, int spikeTimestampsLength):
    ISI = np.array([0]*101)
    cdef int spikeTsLength = spikeTimestampsLength
    cdef double[:] spikeTs = spikeTimestamps.reshape((spikeTimestamps.size,))
    cdef long[:] ISICounts = ISI.reshape((ISI.size,))
    
    binnedISI(&spikeTs[0], &ISICounts[0], spikeTsLength)
    
    return ISICounts
