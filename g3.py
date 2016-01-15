#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.integrate import simps
import argparse
import math
import os

parser = argparse.ArgumentParser()
parser.add_argument('inputs', nargs='+')
parser.add_argument('--suffix', '-s', nargs='?', default='', const='')
parser.add_argument('--output', '-o', nargs='?', default=None, const=None)
parser.add_argument('--xlog', nargs='?', default=False, const=True)
parser.add_argument('--y_rescale_power', nargs=1, default='1')
parser.add_argument('--xlabel', nargs=1, default=None)
parser.add_argument('--ylabel', nargs=1, default=None)
parser.add_argument('--xlims', nargs=2, default=['auto', 'auto'])
parser.add_argument('--ylims', nargs=2, default=['auto', 'auto'])
parser.add_argument('--n_columns', nargs=1, default=['1'])
parser.add_argument('--linestyle', nargs=1, default=['-'])
parser.add_argument('--markerstyle', nargs=1, default=['None'])
parser.add_argument('--markersize', nargs=1, default=['10'])
parser.add_argument('--single_colour', nargs='?', default=None, const=None)
parser.add_argument('--aspect', nargs=1, default=None)
parser.add_argument('--fig_size', nargs=2, default=['auto', 'auto'])
parser.add_argument('--annotation_id', nargs=1, default=None)
parser.add_argument('--save_info', nargs='?', default=None, const=1)

args = parser.parse_args(sys.argv[1:])
#fig=plt.figure()
#q = fig.add_subplot(111, aspect=5)
#fig, ax = plt.subplots(1)
x_dim = 21.0
y_dim = 15.0
print(args.aspect)
if (args.aspect):
    y_dim = float(args.aspect[0]) * x_dim
fig = plt.figure(figsize=(x_dim / 2.5, y_dim / 2.5)) # INCHES TO CM THANKS AMERICA
ax = fig.add_subplot(111)

labelsize = 20


if (args.save_info and args.output):
    f = open(args.output + ".info", 'w+')
    r = ''
    toggle = False
    for a in sys.argv:
        if (toggle):
            r += '"' + a.replace('$', '\$') + '"' + ' '
            toggle = False
        else:
            r += a + ' '
        if ('label' in a):
            toggle = True
    f.write(r)


data_files = [x + args.suffix for x in args.inputs]

n=len(data_files)+1
color=iter(plt.cm.cool(np.linspace(0,1,n)))
c=next(color)

def rescale_x(x):
  return float(x.strip(':'))
def rescale_y(y):
    scale_factor = float(args.y_rescale_power[0])
    if (scale_factor != 1):
        return math.copysign(pow(abs(float(y)), scale_factor), float(y))
    return float(y)

lims = [0,0,0,0]
if (args.xlims[0] != 'auto'):
  lims[0] = float(args.xlims[0])
if (args.xlims[1] != 'auto'):
  lims[1] = float(args.xlims[1])
if (args.ylims[0] != 'auto'):
  lims[2] = float(args.ylims[0])
if (args.ylims[1] != 'auto'):
  lims[3] = float(args.ylims[1])

def print_file(name):
    f = open(name)
    data_raw = []

    rep = '.' + (name.split('.'))[-1]
    length_scale_filename = name.replace(rep, '') + ".length"
    if (args.suffix != ''):
      length_scale_filename = length_scale_filename.replace(args.suffix, '')
    rescale_factor = 1
    if (os.path.isfile(length_scale_filename)):
        length_scale_file = open(length_scale_filename)
        l = length_scale_file.readlines()
        lj = map(float, l[6:9])
        print(lj)
        rescale_sf2 = False 
        if (rescale_sf2 and '.sf2' in name):
          rescale_factor = lj[1] / lj[0]



    for line in f.readlines():
        sp = line.split()
        temp = int(args.n_columns[0]) + 1
        sp = sp[0:temp]
        rescaled = []
        rescaled.append(rescale_factor * rescale_x(sp[0]))
        for y in sp[1:]:
            rescaled.append(rescale_y(y) / rescale_factor)
        data_raw.append(rescaled)
    f.close()

    #with open(name, 'r') as infile:
        #data = np.array([map(float, line.split()) for line in infile])

    data = np.array(data_raw)
    data = data.T
    if (args.single_colour):
      col = args.single_colour
    else:
      col = next(color)

    print(args.linestyle[0])
    ax.plot(data[0], data[1], args.linestyle[0], marker=args.markerstyle[0], color=col, markersize=int(args.markersize[0]))
    #plt.plot(data[0], data[1], args.linestyle[0])

    if (args.xlims[0] == 'auto'): lims[0] = min(lims[0], min(data[0]))
    if (args.xlims[1] == 'auto'): lims[1] = max(lims[1], max(data[0,0], data[0,-1]))
    if (args.ylims[0] == 'auto'): lims[2] = min(lims[2], min(data[1]))
    if (args.ylims[1] == 'auto'): lims[3] = max(lims[3], max(data[1]))
    print(lims)
    int1 = simps(data[1], data[0])
    print(int1)



for f in data_files:
    print_file(f)

#identify what kind of graph we have if possible
ident = None
#is it a structure factor?
for i in [0, 2, 3]:
    test = 'sf' + str(i)
    if (test in args.suffix or test in args.inputs[0]):
        ident = test

ylabel=""
if (args.ylabel):
    ylabel = args.ylabel[0]
else:
    ylabel = ""
    if (ident == "sf3"):
        ylabel = "$4\pi P_{hl}(r)r^2 / a^{-1}_{\mathrm{ho}}$"
    if (ident == "sf2"):
        ylabel = "$4\pi P_{hh}(r)r^2 / a^{-1}_{\mathrm{ho}}$"
    if (ident == "sf0"):
        ylabel = "$ P_{\mathrm{hyper}}(R) / ( \sqrt{2} a^{-1}_{\mathrm{ho}})$"
plt.ylabel(ylabel, fontsize=labelsize)

xlabel=""
if (args.xlabel):
    xlabel = args.xlabel[0]
if (not args.xlabel):
    xlabel = ""
    if (ident == "sf3" or ident == "sf2"):
        xlabel = "$r / a_\mathrm{ho}$"
    if (ident == "sf0"):
        xlabel = "$R / \sqrt{2} a_\mathrm{ho}$"

plt.xlabel(xlabel, fontsize=labelsize)

#if (args.aspect): ax.set_aspect(int(args.aspect[0]), 'box-forced')
plt.axis(lims)

def special_1():
  fsiz = 18
  shrink = 0.03
  width=1.0
  headwidth=6
  ax.annotate('$L = 0$', xy=(2.159, 0.296), xytext=(2.2654, 0.4255), fontsize=fsiz,
                  arrowprops=dict(frac = 0.06, facecolor='black', shrink=shrink, width=width, headwidth=headwidth))
  ax.annotate('$L = 1$', xy=(2.206, 0.232), xytext=(2.45, 0.37), fontsize=fsiz,
                  arrowprops=dict(frac = 0.05, facecolor='black', shrink=shrink, width=width, headwidth=headwidth))
  ax.annotate('$L = 2$', xy=(2.429, 0.247), xytext=(2.65, 0.308), fontsize=fsiz,
                  arrowprops=dict(frac = 0.08, facecolor='black', shrink=shrink, width=width, headwidth=headwidth))
  ax.annotate('$L = 3$', xy=(2.635, 0.233), xytext=(2.888, 0.255), fontsize=fsiz,
                  arrowprops=dict(frac = 0.1, facecolor='black', shrink=shrink, width=width, headwidth=headwidth))
  ax.text(5.2, 0.41, '$\eta = 64$', fontsize = fsiz + 2)

def special_2():
  fsiz = 18
  shrink = 0.03
  width=1.0
  headwidth=6
  ax.annotate('$L = 0$', xy=(2.575, 0.210), xytext=(2.7, 0.25), fontsize=fsiz,
                  arrowprops=dict(frac = 0.08, facecolor='black', shrink=shrink, width=width, headwidth=headwidth))
  ax.annotate('$L = 1, \, 2,\, 3$', xy=(1.76, 0.38), xytext=(1.87, 0.44), fontsize=fsiz,
                  arrowprops=dict(frac = 0.1, facecolor='black', shrink=shrink, width=width, headwidth=headwidth))
  ax.text(4.3, 0.45, '$\eta = 1$', fontsize = fsiz + 2)


print args.annotation_id
if (args.annotation_id):
  if (args.annotation_id[0] == 'fig_L_64'):
    special_1()
  if (args.annotation_id[0] == 'fig_L_1'):
    special_2()
#q.axis(lims)
#ax = plt.gca()
#ax.ticklabel_format(useOffset=False)
plt.ticklabel_format(useOffset=False)
#plt.axes().set_aspect(float(args.aspect))
if (not args.output):
  plt.show()
else:
  plt.savefig(args.output, bbox_inches='tight')

