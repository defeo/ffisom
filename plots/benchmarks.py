#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import re, math
from functools import reduce

def prime_powers(bound):
    '''
    Very naive sieve outputting the lists of primes and prime powers up to bound.
    '''
    
    primes = list(range(2, bound + 1))
    ppowers = list(range(2, bound + 1))
    for p in primes:
        for n in range(2*p, bound+1, p):
            try:
                primes.remove(n)
            except:
                try:
                    ppowers.remove(n)
                except:
                    pass
    return primes, ppowers

def parse_timings(datasets=['04/{:0>2}'.format(i) for i in range(1,12)] +
                      ['05/{:0>2}'.format(i) for i in [1,2,3,4,5,6,7,8,10]]):
    '''
    Parse timing data from files. Returns a pandas DataFrame.
    '''
    def read_file(f):
        d = pd.read_table(f, sep=r'[,\)]? \(?', comment='#', engine='python',
                            names=['prime', 'degree', 'rains_aux', 'allombert_aux',
                                       't_rains', 't_conic_rains', 't_ellrains',
                                       't_pari', 't_allom_linalg', 't_allom_modcomp',
                                       't_allom_cofactor', 't_allom_iterfrob'],
                            na_values='0')
        d['dataset'] = f
        return d
    d = pd.concat([read_file('benchdata/%s.dat' % run) for run in datasets])
    return d

def plot_rains(d, size=(10,10)):
    '''
    Plot cyclotomic Rains' vs conic Rains' vs elliptic Rains' algorithm.
    '''
    
    # Figure settings
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    ax.loglog(basex=2, basey=2)
    ax.set_xlabel('degree $r$')
    ax.set_ylabel('seconds')

    # Plot cyclotomic Rains', one line per auxiliary extension degree up to 9
    df = d.groupby([d.rains_aux, d.degree]).median().loc[:9]
    ax.set_xlim(df.index.min()[1], df.index.max()[1])
    ax.set_ylim(df.t_rains.min(), df.t_rains.max())
    for key, g in df.groupby(level=0):
        ax.plot(g.index.get_level_values('degree'), g.t_rains, label=key, color='blue')
    # Plot conic Rains'
    df = d[~d.t_conic_rains.isnull()]
    df = df.groupby(df.degree).median()
    ax.plot(df.index, df.t_conic_rains, 'g')
    # Plot elliptic Rains', scatter style
    ax.plot(d.degree, d.t_ellrains, 'r.', ms=2)

    # Legend
    cyclo = mlines.Line2D([], [], color='blue', label="Cyclotomic Rains'")
    conic = mlines.Line2D([], [], color='green', label="Conic Rains'")
    ell = mlines.Line2D([], [], color='red', label="Elliptic Rains'")
    ax.legend(handles=[cyclo, conic, ell], loc=2)

    return fig

def plot_allombert_lowaux(d, size=(10,10)):
    '''
    Plot various implementations of Allombert's algorithm for small auxiliary degree.
    '''
    
    # Figure settings
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    ax.set_xlabel('degree $r$')
    ax.set_ylabel('seconds')

    # Plot this columns, with legend and z-order
    cols = [
        ('Case 1', 't_allom_modcomp', 3),
        ('Case 2', 't_allom_cofactor', 4),
        ('Case 2 (variant)', 't_allom_linalg', 2),
        ('Case 3', 't_allom_iterfrob', 0),
        ('PARI/gp', 't_pari', 1)
    ]
    # Only plot if aux degree is <= 10 and all algorithms were run
    df = d[(d.allombert_aux <= 10) & reduce(lambda x,y: x & ~d[y[1]].isnull(), cols, True)]
    
    ax.set_xlim(df.degree.min(), df.degree.max()+5)
    ax.set_ylim(0, df.t_pari.max())

    for name, c, z in cols:
        data = df[c].values
        # scatter plot
        scat, = ax.plot(df.degree, data, '.', alpha=0.4, zorder=z)
        # degree 2 linear regression
        linreg = np.polyfit(df.degree, data, 2)
        x = np.linspace(100, df.degree.max(), 100)
        ax.plot(x, np.polyval(linreg, x), color=scat.get_color(), label=name, zorder=z)
        
    ax.legend(loc=2)

    return fig

def plot_allombert_anyaux(d, size=(10,10)):
    '''
    Plot various implementations of Allombert's algorithm with respect to auxiliary degree.
    '''
    
    # Figure settings
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'order of $q$ mod $r$')
    ax.set_ylabel(r'ratio $\times 10^4$')

    # Add a reference column for scaling timings with respect to degree^2
    df = pd.concat([d, pd.Series(d.degree**2 * 10**-4, name='deg_sq')], axis=1)
    
    def plot_algo(col, label, marker='.', *args, **kwds):
        dd = df[~df[col].isnull()]
        data = dd[col] / dd.deg_sq
        # scatter plot
        scat,  = ax.plot(dd.allombert_aux, data, *args, ls='none',
                             marker=marker, alpha=0.4, ms=2, **kwds)
        # degree 2 linear regression
        linreg = np.polyfit(dd.allombert_aux, data, 2)
        x = np.linspace(1, dd.allombert_aux.max(), 100)
        ax.plot(x, np.polyval(linreg, x), *args, lw=1,
                    color=scat.get_color(), label=label, **kwds)

    # Plot the algorithms
    plot_algo('t_allom_modcomp', label="Case 1", zorder=1)
    plot_algo('t_allom_cofactor', label="Case 2", zorder=0)
    plot_algo('t_allom_linalg', label="Case 2 (variant)", zorder=2)
    plot_algo('t_allom_iterfrob', label="Case 3", zorder=3)
    plot_algo('t_pari', marker='x', label="PARI/gp", zorder=4)

    # Add a vertical line for the reference
    ax.plot([0, df.allombert_aux.max()+10], [1, 1], 'k--', label='$r^2$')
    
    ax.set_xlim(0, df.allombert_aux.max()+10)
    ax.set_ylim(0, 2*(df.t_allom_iterfrob / df.deg_sq).max())

    # Legend
    ax.legend(loc=1)
    
    return fig

def plot_all(d, size=(10,10)):
    '''
    Comparison of all algorithms by degree.
    '''
    
    # Figure settings
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    ax.loglog(basex=2, basey=2)
    ax.set_xlabel('degree $r$')
    ax.set_ylabel('seconds')

    # Plot Allombert's Case 2 variant
    # Group auxiliary degree by successive powers of 8
    b = 8
    partitions = [(b**i, b**(i+1)) for i in range(int(math.log(d.allombert_aux.max(), b) + 1))]
    for l, h in partitions:
        df = d[(d.allombert_aux >= l) & (d.allombert_aux < h)
                   & ~d.t_allom_cofactor.isnull()].groupby('degree').t_allom_cofactor
        # plot median time as a line
        median = df.median()
        pl, = ax.plot(median.index, median, alpha=0.3,
                          label=r"Allombert's (2) $s\in[%d,%d]$" % (l, h-1))
        # and min/max as an area
        ax.fill_between(median.index, df.min(), df.max(), alpha=0.1,
                            facecolor=pl.get_color(), edgecolor='none')

    # Plot cyclotomic Rains', auxiliary degree 1
    df = d[d.rains_aux==1]
    df = df.groupby(df.degree).median()
    ax.plot(df.index, df.t_rains, lw=2, color='black', label="Cyclotmic Rains' $s=1$")
    # Plot conic Rains'
    df = d[~d.t_conic_rains.isnull()]
    df = df.groupby(df.degree).median()
    ax.plot(df.index, df.t_conic_rains, lw=2, color='purple', label="Conic Rains'")
    # Plot elliptic Rains', scatter style
    df = d.sort_values('degree')
    ax.plot(df.degree, df.t_ellrains, '.', ms=1, label="Elliptic Rains'")
    
    # Legend
    ax.legend(loc=2)
    
    return fig

def pdf_plots(prefix='bench-'):
    '''
    Plot benchmarks.
    '''
    d = parse_timings()
    
    # Take out non-prime-powers to limit noise
    primes, ppowers = prime_powers(d.degree.max())
    d = d[d.degree.isin(ppowers) & (d.prime > 103) & (d.prime < 2**13)]

    # Global figure configurations
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', serif='Computer Modern')
    plt.rc('legend', fontsize=10)

    # Save pdf plots
    plot = plot_rains(d, (4,3))
    plot.savefig(prefix + 'rains.pdf', bbox_inches='tight')

    plot = plot_allombert_lowaux(d, (5,4))
    plot.savefig(prefix + 'allombert-lowaux.pdf', bbox_inches='tight')
    
    plot = plot_allombert_anyaux(d, (5,4))
    plot.savefig(prefix + 'allombert-anyaux.pdf', bbox_inches='tight')

    allombert = plot_all(d, (5,4))
    allombert.savefig(prefix + 'all.pdf', bbox_inches='tight')

    
if __name__ == '__main__':
    pdf_plots()
