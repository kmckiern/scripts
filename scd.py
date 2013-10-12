#!/usr/bin/env python

""" 
This script calculates S_CD(t), assuming a perfectly tetrahedral geometry, from a Gromacs index (input 1) and trajectory (input 2) file for a lipid bilayer.
Last revision: Oct 10 2013
"""

import sys, string
import numpy as np
import math

"""
Parsing/recombining index and trajectory data in a convenient way.
"""
# Parses a data file and extracts a list of demarcation elements.
def c_ls(d_f,line_start):
	data = open(d_f).readlines()
	c_ls = []
	for line in data:
		if line[0] == line_start:
			l, d, r = string.split(line)
			c_ls.append(d)
	return c_ls

# Reads in a trajectory file and stores the coordinates of the lipid tail atoms as a function of snapshot time.
def trj_info(ndx_f,trj_f):
	snap_ts = []
	trj_qs = {}
	tail_atms = c_ls(ndx_f,'[')
	c_0 = tail_atms[0]
	coords = {c: [] for c in tail_atms}
	trj = open(trj_f).readlines()
	c_ndx = 0
	t_ndx = 0
	for line in trj:
		l_a = line.split()
		if len(l_a) == 8:
			snap_t = l_a[-1]
			snap_ts.append(snap_t)
			trj_qs[snap_t] = {}
			if len(snap_ts) > 1:
				# Add the coordinate info from the previous block for coordinate data.
				trj_qs[snap_ts[t_ndx]] = coords
				# Reset coordinate info for new snapshot.
				coords = {c: [] for c in tail_atms}
				t_ndx = t_ndx+1
		while len(l_a) == 6:
			if l_a[1] == tail_atms[c_ndx]:
				atm_nm = l_a[1]
				x_c = l_a[-3]
				y_c = l_a[-2]
				z_c = l_a[-1]
				coords[atm_nm].append([x_c,y_c,z_c])
				if c_ndx < len(tail_atms)-1:
					c_ndx = c_ndx+1
				else:
					c_ndx = 0
			break
	trj_qs[snap_ts[t_ndx]] = coords
	return trj_qs, snap_ts, c_0

"""
Define some vector operations.
"""
# For a given atom, creates a dictionary which specifies its two neighbors.
def neighbors(n):
	c_list = c_ls(n,'[')
	cs = len(c_list)
	n_ls = []
	for c in range(1,cs-1):
		n_ls.append([c_list[c-1], c_list[c], c_list[c+1]])
	return n_ls

# Finds the modulus of a vector.
def modulus(v):
<<<<<<< HEAD
	modulus = np.sqrt(np.dot(v,v))
=======
	modulus = np.absolute(np.sqrt(np.dot(v,v)))
>>>>>>> d023e608b2f20cd1d0210cd7a3c67c42bd070067
	return modulus

# Normalizes a vector.
def unit(v):
	u = v / modulus(v)
	return u

# Finds the angle between two vectors.
def angle(v1,v2):
	dp = np.dot(v1,v2)
	ang = np.arccos(dp / (modulus(v1) * modulus(v2)))
	return ang

# Finds the distance between two vectors.
def distance(v1,v2):
	v1 = np.array(v1)
	v2 = np.array(v2)
	dis = np.sqrt((v2[2]-v1[2])**2 + (v2[1]-v1[1])**2 + (v2[0]-v1[0])**2)
	return dis

# Rotates a 3d vector counterclockwise about an arbitrary 3d axis.
def rotation(v,axis_v,theta):
	u = unit(axis_v)

	# Axis components.
	ux = u[0]
	uy = u[1]
	uz = u[2]
	# Cross product matrix of u.
	cpm = np.array([[0, -uz, uy],[uz, 0, -ux],[-uy, ux, 0]])
	# Tensor product of u with u.
	uutpm = np.array([[ux**2, ux*uy, ux*uz],[ux*uy, uy**2, uy*uz],[ux*uz, uy*uz, uz**2]])
	# Identity matrix.
	im = np.array([[1,0,0],[0,1,0],[0,0,1]])
	# Rodrigues' rotation formula.
	R = im*np.cos(theta) + np.sin(theta)*cpm + (1-np.cos(theta))*uutpm

	rotated_v = np.dot(R,np.transpose(v))
	return rotated_v

"""
Calculates deuterium locations via some rotations.  Finds the angle between the C-D axis and the bilayer normal.  Calculates the order parameter from this angle and stores it in a dictionary indexed by snapshot time and carbon position.
"""
def s_cd(ndx_f,trj_f):
	n_ls = neighbors(ndx_f)

	trj_prop = trj_info(ndx_f,trj_f)
	trj_qs = trj_prop[0]
	snap_ts = trj_prop[1]
	c_0 = trj_prop[2]
	n_snaps = len(snap_ts)

	td_scd = {}

	phi = 0.9557
	z = np.array([0,0,1])

<<<<<<< HEAD
	for snap_t in range(0,n_snaps):
=======
	for snap_t in range(n_snaps-10,n_snaps):
>>>>>>> d023e608b2f20cd1d0210cd7a3c67c42bd070067
		tail_length = len(trj_qs[snap_ts[snap_t]])
		s_cd = {}

		for c in range(len(n_ls)):
			scd1s = []
			scd2s = []
			tail_c_qs = trj_qs[snap_ts[snap_t]]
			n_qs = len(tail_c_qs[n_ls[c][1]])
			nbrs = n_ls[c]

			for l in range(n_qs):
				# Find the coordinates of the carbon of interest and its neighbors.
				cq_prior = np.array(map(float, tail_c_qs[nbrs[0]][l]))
				cq = np.array(map(float, tail_c_qs[nbrs[1]][l]))
				cq_post = np.array(map(float, tail_c_qs[nbrs[2]][l]))

				# Find the vector pointing from the carbon of interest to the carbon before it.
				v_prior = cq - cq_prior
				# Find the vector pointing from the carbon of interest to the carbon after it.
				v_post = cq_post - cq

				# Define the plane containing v_prior and v_post by taking the cross product.
				nrml_v = np.cross(v_prior,v_post)

				# Find the angle between v_prior and v_post.
				angl = angle(v_prior,v_post)

<<<<<<< HEAD
				# Rotate v_post about the nrml_v to bisect the v plane, and then invert to bisect the c-d plane.
				cd_bisector = (-.5)*rotation(v_post,nrml_v,-1*angl/2)

				# Rotate v_post so that it is perpendicular to the c_v normal vector.
				beta = ((np.pi - angl)/2)
				vv_plane_x = rotation(v_post,nrml_v,beta)
				# Rotate v_post_to_cd_plane to the cd_axes.
				cd_axis1 = rotation(cd_bisector,vv_plane_x,phi)
				cd_axis2 = rotation(cd_bisector,vv_plane_x,-1*phi)
=======
				# Rotate v_post into the C-D plane.
				alpha = -1*(np.pi - (angl/2))
				v_post_to_cd_plane = rotation(v_post,nrml_v,alpha)
				# Rotate v_post so that it is perpendicular to the c_v normal vector.
				beta = -1*((np.pi - angl)/2)
				vv_plane_x = rotation(v_post,nrml_v,beta)
				# Rotate v_post_to_cd_plane to the cd_axes.
				cd_axis1 = rotation(v_post_to_cd_plane,vv_plane_x,phi)
				cd_axis2 = rotation(v_post_to_cd_plane,vv_plane_x,-1*phi)
>>>>>>> d023e608b2f20cd1d0210cd7a3c67c42bd070067

				# Find angle between CD axis and the global normal vector of the bilayer.
				beta1 = angle(cd_axis1,z)
				beta2 = angle(cd_axis2,z)

				# Find S_CD using beta.
				scd1 = np.absolute((3*np.cos(beta1)**2 - 1)/2)
				scd2 = np.absolute((3*np.cos(beta2)**2 - 1)/2)
				scd1s.append(scd1)
				scd2s.append(scd2)
				
			s_cd.update({n_ls[c][1] : [np.average(scd1s),np.average(scd2s)]})
		td_scd.update({snap_t : s_cd})
	return td_scd

def main():
	n = sys.argv[1] 
	t = sys.argv[2]

	scd = s_cd(n,t)
	print scd

if __name__ == "__main__":
	main()
