#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib
if not matplotlib.get_backend():
    matplotlib.use('Agg')
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

def parse_timings(datasets=['12/{:0>2}'.format(i) for i in filter(lambda x: x != 36 and x != 39,
                                                                      range(1,41))]):
    '''
    Parse timing data from files. Returns a pandas DataFrame.
    '''
    def read_file(f):
        d = pd.read_table(f, sep=r'[,\)]? \(?', comment='#', engine='python',
                            names=['prime', 'degree', 'rains_aux', 'kummer_aux',
                                       't_magma', 't_cyclo_rains', 't_conic_rains',
                                       't_elliptic_rains', 't_pari', 't_kummer_cyclo_linalg',
                                       't_kummer_linalg_only', 't_kummer_linalg',
                                       't_kummer_modcomp', 't_kummer_cofactor',
                                       't_kummer_iterfrob', 't_kummer_mpe'],
                            na_values='0')
        d['dataset'] = f
        return d
    d = pd.concat([read_file('benchdata/%s.dat' % run) for run in datasets])
    return d

def plot_char(d, size=(10,10)):
    '''
    Plot something as p grows
    '''
    
    # Figure settings
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    ax.loglog(basex=2, basey=2)
    ax.set_xlabel('Characteristic')
    ax.set_ylabel('seconds')

    df = d[(d.degree <= 200) & (d.kummer_aux != 0)]
    data = df.t_pari / df.kummer_aux
    ax.plot(df.prime, data, '.')
    data = df.t_kummer_cofactor / df.kummer_aux
    ax.plot(df.prime, data, '.')

    return fig
    

def plot_magma(d, size=(10,10)):
    '''
    Plot our implementation of cyclotomic Rains' vs Magma
    '''
    
    # Figure settings
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    ax.loglog(basex=2, basey=2)
    ax.set_xlabel("Cyclotomic Rains' (seconds)")
    ax.set_ylabel('Magma (ratio)')

    df = d[~d.t_magma.isnull() & ~d.t_cyclo_rains.isnull()]
    data = df.t_cyclo_rains
    ax.plot(df.t_cyclo_rains, data, '.', ms=1, alpha=0.5)

    return fig

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
    ax.set_ylim(df.t_cyclo_rains.min(), df.t_cyclo_rains.max())
    for key, g in df.groupby(level=0):
        ax.plot(g.index.get_level_values('degree'), g.t_cyclo_rains, label=key, color='m')
    # Plot conic Rains'
    df = d[~d.t_conic_rains.isnull()]
    df = df.groupby(df.degree).t_conic_rains.median()
    ax.plot(df.index, df, color='y')
    # Plot elliptic Rains'
    df = d[~d.t_elliptic_rains.isnull()]
    df = df.groupby(df.degree).t_elliptic_rains.median()
    ax.plot(df.index, df, color='k')

    # Legend
    cyclo = mlines.Line2D([], [], color='m', label="Cyclotomic Rains'")
    conic = mlines.Line2D([], [], color='y', label="Conic Rains'")
    ell = mlines.Line2D([], [], color='k', label="Elliptic Rains'")
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
        ('Divide \& conquer', 't_kummer_modcomp', 4),
        ('Automorphism eval.', 't_kummer_cofactor', 5),
        ('Multipoint eval.', 't_kummer_mpe', 0),
        ('Multipoint eval. (var)', 't_kummer_iterfrob', 2),
        ('PARI/GP', 't_pari', 1),
        ('Allombert (rev)', 't_kummer_linalg_only', 3),
#        ('cyclo linalg', 't_kummer_cyclo_linalg', 6),
#        ('linalg', 't_kummer_linalg', 7),
    ]
    # Only plot if aux degree is <= 10 and all algorithms were run
    df = d[(d.kummer_aux <= 10) & reduce(lambda x,y: x & ~d[y[1]].isnull(), cols, True)]
    
    ax.set_xlim(df.degree.min(), df.degree.max()+5)
    ax.set_ylim(0, df.t_kummer_iterfrob.max())

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
    ax.set_ylabel(r'ratio')

    # Add a reference column for scaling timings with respect to degree^2
    df = pd.concat([d, pd.Series(d.degree**2 * 10**-4, name='deg_sq')], axis=1)
    
    def plot_algo(col, label, marker='.', *args, **kwds):
        dd = df[~df[col].isnull()]
        data = dd[col] / dd.deg_sq
        # scatter plot
        scat,  = ax.plot(dd.kummer_aux, data, *args, ls='none',
                             marker=marker, alpha=0.4, ms=2, **kwds)
        # degree 2 linear regression
        linreg = np.polyfit(dd.kummer_aux, data, 2)
        x = np.linspace(1, dd.kummer_aux.max(), 100)
        ax.plot(x, np.polyval(linreg, x), *args, lw=1,
                    color=scat.get_color(), label=label, **kwds)

    # Plot the algorithms
    plot_algo('t_kummer_modcomp', label="Divide \& conquer", zorder=2)
    plot_algo('t_kummer_cofactor', label="Automorphism eval.", zorder=0)
    plot_algo('t_kummer_mpe', label="Multipoint eval.", zorder=4)
    plot_algo('t_kummer_iterfrob', label="Multipoint eval. (var)", zorder=1)
    plot_algo('t_pari', label="PARI/GP", zorder=5)
    plot_algo('t_kummer_linalg_only', label="Allombert (rev)", zorder=3)
#    plot_algo('t_kummer_cyclo_linalg', label="cyclo linalg", zorder=6)
#    plot_algo('t_kummer_linalg', label="linalg", zorder=7)

    # Add a vertical line for the reference
    ax.plot([0, df.kummer_aux.max()+10], [1, 1], 'k--', label='$r^2$')
    
    ax.set_xlim(0, df.kummer_aux.max()+10)
    ax.set_ylim(0, (df.t_kummer_modcomp / df.deg_sq).max())

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
    partitions = [(b**i, b**(i+1)) for i in range(int(math.log(d.kummer_aux.max(), b) + 1))]
    for l, h in partitions:
        df = d[(d.kummer_aux >= l) & (d.kummer_aux < h)
                   & ~d.t_kummer_cofactor.isnull()].groupby('degree').t_kummer_cofactor
        # plot median time as a line
        median = df.median()
        pl, = ax.plot(median.index, median, alpha=0.3,
                          label=r"Allombert (AE) $s\in[%d,%d]$" % (l, h-1))
        # and min/max as an area
        ax.fill_between(median.index, df.min(), df.max(), alpha=0.1,
                            facecolor=pl.get_color(), edgecolor='none')

    # Plot cyclotomic Rains', auxiliary degree 1
    df = d[d.rains_aux==1]
    df = df.groupby(df.degree).median()
    ax.plot(df.index, df.t_cyclo_rains, lw=2,  label="Cyclotmic Rains' $s=1$")
    # Plot conic Rains'
    df = d[~d.t_conic_rains.isnull()]
    df = df.groupby(df.degree).median()
    ax.plot(df.index, df.t_conic_rains, lw=2, label="Conic Rains'")
    # Plot elliptic Rains'
    df = d[~d.t_elliptic_rains.isnull()]
    df = df.groupby(df.degree).median()
    ax.plot(df.index, df.t_elliptic_rains, lw=2, label="Elliptic Rains'")
    
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
    d = d[d.degree.isin(ppowers) & (d.prime > 100) & (d.prime < 2**20)]

    # Global figure configurations
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', serif='Computer Modern')
    plt.rc('legend', fontsize=10)

    # Save pdf plots
    plot = plot_magma(d, (5,4))
    plot.savefig(prefix + 'magma.pdf', bbox_inches='tight')
    
    plot = plot_rains(d, (4,3))
    plot.savefig(prefix + 'rains.pdf', bbox_inches='tight')

    plot = plot_allombert_lowaux(d, (5,4))
    plot.savefig(prefix + 'allombert-lowaux.pdf', bbox_inches='tight')
    
    plot = plot_allombert_anyaux(d, (5,4))
    plot.savefig(prefix + 'allombert-anyaux.pdf', bbox_inches='tight')

    allombert = plot_all(d[d.degree < 2**10], (6,5))
    allombert.savefig(prefix + 'all.pdf', bbox_inches='tight')

    
if __name__ == '__main__':
    pdf_plots()
