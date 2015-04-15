#!/usr/bin/env python

# todo: obv make more general

from forcebalance.nifty import lp_load
import numpy as np

x = lp_load('npt_result.p')
al = x[-7]
t = np.arange(len(al))
np.savetxt('al.dat', zip(t,al))
