#!/usr/bin/env python

import pickle

pkl = open("forcebalance.p","rb")

data = pickle.load(pkl)

ff, mvals, h, agrad = data

ff.make(mvals)
