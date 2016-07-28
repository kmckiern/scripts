import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

df = pd.read_pickle('forces/PROJ9712_RUN3_CLONE18_310.pkl')
df = df.loc[df['resName'] == 'K+']

qs = np.array([np.array(i) for i in df['q']])
x, y, z = qs.T

f_nb = np.array([np.array(i) for i in df['NonbondedForce']])
u, v, w = f_nb.T

fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.quiver(x, y, z, u, v, w, length=1)
ax = fig.gca()
ax.quiver(x, z, u, w)
plt.savefig('field_xz.png')
