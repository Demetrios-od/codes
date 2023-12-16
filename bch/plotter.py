#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


folder = './res/'
fig = plt.figure()
ax = fig.add_subplot(111)

def add_line(fname, legend, fmt='', line_prop={}, xcol_name='SNR', ycol_name='BER'):
	data = pd.read_csv(folder + fname, sep='\s+', comment='#')
	line_prop.update({'data': data, 'label': legend})
	ax.plot(xcol_name, ycol_name, fmt, **line_prop)


# line params:
#     linestyle = ls: '-', '--', '-.', ':', ' '
#     linewidth = lw
#     color = c: b, g, r, c, m, y, k, w
#     gapcolor
#     marker: '.', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X'
#     markersize = ms
#     markeredgewidth = mew
#     markeredgecolor = mec
#     markerfacecolor = mfc
#     markerfacecoloralt = mfcalt
#     fillstyle: 'full', 'left', 'right', 'bottom', 'top', 'none'
#     antialiased = aa: bool
#     solid_capstyle: 'butt', 'projecting', 'round'
#     dash_capstyle:  'butt', 'projecting', 'round'
#     solid_joinstyle: 'miter', 'round', 'bevel'
#     dash_joinstyle:  'miter', 'round', 'bevel'
#     pickradius
#     drawstyle = ds: 'default', 'steps', 'steps-pre', 'steps-mid', 'steps-post'
#     markevery: None or int or (int, int) or slice or list[int] or float or (float, float) or list[bool]


add_line('BCH256-Chase-p4.txt', 'Chase p4', 'b')
add_line('BCH256-Chase-p5.txt', 'Chase p5', 'r')

ax.set(
	title = 'Performance BCH-256',
	yscale = 'log',
	xlabel = 'SNR, dB',
	ylabel = 'BER',
	xlim = [7.6, 10],
	ylim = [1e-7, 0.01],
)
ax.grid(
	which = 'major',
	color = [0.7, 0.7, 0.7],
	linestyle = '-',
)
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_major_locator(ticker.FixedLocator([10**i for i in range(-2,-16,-1)]))
ax.legend()
fig.set_figwidth(10)
fig.set_figheight(9)
plt.show()
