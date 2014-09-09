#!/usr/bin/env python

import numpy as np
import sys

def avg_data(d_f):
	data = np.loadtxt(d_f)
	vs = []
	print len(data)
	for i in range(len(data)):
		vs.append(data[i][1])
	avg = np.average(vs)
	return avg	

avg = avg_data(sys.argv[1])
print avg
