#!/usr/bin/python3.5

import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
import sys

lines = []

def main():

    file_name = sys.argv[1]
    print(file_name)
    save_file_name = sys.argv[2]
    print(save_file_name)
    num_modes = 0
    num_dims = 0
    with open(file_name, 'r') as file:
        num_modes, num_dims = next(file).split()
        for i in range(int(num_modes)):
            amp, freq, damp, l_idx, r_idx = next(file).split()
            lines.append((int(l_idx), int(r_idx)))
    num_modes = int(num_modes)
    num_dims = int(num_dims)

    names = ['f', 'real', 'found']
    for i in range(num_modes):
        names.append('found_%s' % str(i+1))
    names.append('gauss')
    print(names)
    data = np.genfromtxt(file_name, skip_header=num_modes+1, delimiter=';', names=names)

    fig = plt.figure()

    ax = fig.add_subplot(111)
    for rect in lines:
        ax.add_patch(patches.Rectangle((rect[0], 0.0), rect[1]-rect[0], 0.01, alpha=0.1))
    ax.plot(data['f'], data['found'], color='b', label='found data')
    for i in range(num_modes):
        ax.plot(data['f'], data['found_%s' % str(i+1)], label='found data %s' % str(i+1))
    ax.plot(data['f'], data['gauss'], label='gauss')
    legend = ax.legend(loc='upper right', shadow=True)

    fig.savefig(save_file_name, bbox_inches='tight')
    fig.savefig(save_file_name)


    interactive_legend().show()

def interactive_legend(ax=None):
    if ax is None:
        ax = plt.gca()
    if ax.legend_ is None:
        ax.legend()

    return InteractiveLegend(ax.legend_)

class InteractiveLegend(object):
    def __init__(self, legend):
        self.legend = legend
        self.fig = legend.axes.figure

        self.lookup_artist, self.lookup_handle = self._build_lookups(legend)
        self._setup_connections()

        self.update()

    def _setup_connections(self):
        for artist in self.legend.texts + self.legend.legendHandles:
            artist.set_picker(10) # 10 points tolerance

        self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)

    def _build_lookups(self, legend):
        labels = [t.get_text() for t in legend.texts]
        handles = legend.legendHandles
        label2handle = dict(zip(labels, handles))
        handle2text = dict(zip(handles, legend.texts))

        lookup_artist = {}
        lookup_handle = {}
        for artist in legend.axes.get_children():
            if artist.get_label() in labels:
                handle = label2handle[artist.get_label()]
                lookup_handle[artist] = handle
                lookup_artist[handle] = artist
                lookup_artist[handle2text[handle]] = artist

        lookup_handle.update(zip(handles, handles))
        lookup_handle.update(zip(legend.texts, handles))

        return lookup_artist, lookup_handle

    def on_pick(self, event):
        handle = event.artist
        if handle in self.lookup_artist:
            artist = self.lookup_artist[handle]
            artist.set_visible(not artist.get_visible())
            self.update()

    def on_click(self, event):
        if event.button == 3:
            visible = False
        elif event.button == 2:
            visible = True
        else:
            return

        for artist in self.lookup_artist.values():
            artist.set_visible(visible)
        self.update()

    def update(self):
        for artist in self.lookup_artist.values():
            handle = self.lookup_handle[artist]
            if artist.get_visible():
                handle.set_visible(True)
            else:
                handle.set_visible(False)
        self.fig.canvas.draw()

    def show(self):
        plt.show()

if __name__ == "__main__":
    main()
