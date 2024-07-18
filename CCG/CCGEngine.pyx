# cython: language_level=3
import numpy as np
cimport numpy as np
from libc.math cimport floor

def CCGEngine(np.ndarray[np.float64_t, ndim=1] times, 
              np.ndarray[np.int32_t, ndim=1] ids, 
              double binSize, 
              int halfBins):
    cdef int nTimes = times.shape[0]
    cdef int nIDs = ids.shape[0]
    cdef int nBins = 2 * halfBins + 1
    cdef double maxDt = binSize * (halfBins + 0.5)
    cdef int i, size, index, current, previous, next, id1, id2, bin
    cdef double t1, t2
    
    if nTimes != nIDs:
        raise ValueError("Times and IDs have different lengths.")
    
    # Number of different IDs
    cdef int maxID = 0
    for i in range(nTimes):
        if ids[i] > maxID:
            maxID = ids[i]
    
    # Output parameter
    size = maxID * maxID * nBins
    cdef np.ndarray[np.float64_t, ndim=1] ccg = np.zeros(size, dtype=np.float64)
    
    # Loop through events
    for current in range(nTimes):
        id1 = ids[current]
        t1 = times[current]
        
        # Previous events
        for previous in range(current-1, -1, -1):
            id2 = ids[previous]
            t2 = times[previous]
            
            # Make sure the event falls in the time window
            if t1 - t2 > maxDt:
                break
            
            # Add an event pair in the CCG
            bin = halfBins + int(floor(0.5 + (t2 - t1) / binSize))
            index = nBins * maxID * (id1 - 1) + nBins * (id2 - 1) + bin
            if index < 0 or index >= size:
                raise IndexError("Index out of bounds")
            ccg[index] += 1
        
        # Next events
        for next in range(current+1, nTimes):
            id2 = ids[next]
            t2 = times[next]
            
            # Make sure the event falls in the time window
            if t2 - t1 >= maxDt:
                break
            
            # Add an event pair in the CCG
            bin = halfBins + int(floor(0.5 + (t2 - t1) / binSize))
            index = nBins * maxID * (id1 - 1) + nBins * (id2 - 1) + bin
            if index < 0 or index >= size:
                raise IndexError("Index out of bounds")
            ccg[index] += 1
    
    return ccg
