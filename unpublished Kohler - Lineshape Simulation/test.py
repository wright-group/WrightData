import numpy as np

import matplotlib.pyplot as plt

import WrightTools as wt

fig, gs = wt.artists.create_figure(cols=[1, 1])

ax = plt.subplot(gs[0, 0])
p = r'measured\dpr 0.5 TO all w1 w2 d1 d2 data.hdf5'
d = wt.kit.read_h5(p)
mine = d['arr']
arr = d['arr']
xi = d['w1']
yi = d['w2']
zi = np.sqrt(arr[:, :, 10, 10])
X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
ax.pcolor(X, Y, Z, cmap=wt.artists.colormaps['default'])
ax.set_xlim(xi.min(), xi.max())
ax.set_ylim(yi.min(), yi.max())
ax.grid()

ax = plt.subplot(gs[0, 1])
p = r'C:\Users\blais\Desktop\TrEE run original\dpr 0.5 TO all w1 w2 d1 d2 arr.npz'
d = np.load(p)
dan = d['arr']
arr = d['arr']
xi = d['w1']
yi = d['w2']
zi = np.sqrt(arr[:, :, 10, 10])
X, Y, Z = wt.artists.pcolor_helper(xi, yi, zi)
ax.pcolor(X, Y, Z, cmap=wt.artists.colormaps['default'])
ax.set_xlim(xi.min(), xi.max())
ax.set_ylim(yi.min(), yi.max())
ax.grid()

plt.savefig('out.png')

