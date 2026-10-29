"""Saved-output visualization; importing this module never starts a solver."""
from pathlib import Path
import numpy as np

def plot_run(directory):
    p=Path(directory)
    with open(p/'history.csv') as f:columns=f.readline().strip().split(',')
    if 'radius_rsun' in columns:  # Ball history (quasistatic.history_row): other columns and profile keys
        raise ValueError('qstars plot does not support quasi-static runs; see history.csv and profile_*.npz')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    h=np.genfromtxt(p/'history.csv',delimiter=',',names=True)
    h=np.atleast_1d(h)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.grid':True,'grid.alpha':.15,'figure.dpi':130})
    fig,ax=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
    t=h['time_years']
    ax[0,0].plot(t,h['bh_mass_msun']-h['bh_mass_msun'][0],color='#146a79');ax[0,0].set(ylabel='BH mass gained [solar masses]')
    ax[0,1].plot(t,h['radius_cm']/1e14,color='#784bb2');ax[0,1].set(ylabel='Outer radius [10¹⁴ cm]')
    ax[0,2].plot(t,h['luminosity_erg_s']/1e42,label='Boundary outflow',color='#df8832')
    ax[0,2].plot(t,h['engine_luminosity_erg_s']/1e42,label='BH deposition',color='#146a79')
    ax[0,2].set(ylabel='Luminosity [10⁴² erg/s]');ax[0,2].legend(fontsize=8)
    ax[1,0].plot(t,h['energy_residual_relative'],label='Energy / initial U')
    ax[1,0].plot(t,h['mass_residual_relative'],label='Mass / initial mass');ax[1,0].legend(fontsize=8);ax[1,0].set(ylabel='Signed conservation residual')
    for k,label in [('supplied_mass','Delivered'),('wind_mass','Escaping wind'),('accreted_rest_mass','BH rest mass')]:
        ax[1,1].plot(t,h['ledger_'+k]/1.98847e33,label=label)
    ax[1,1].set(ylabel='Cumulative transported mass [solar masses]');ax[1,1].legend(fontsize=8)
    for k,label in [('central_generated','BH source'),('radiated_energy','Radiated'),('supplied_energy','Supplied (signed)'),('wind_energy','Wind (signed)')]:
        ax[1,2].plot(t,h['ledger_'+k]/1e48,label=label)
    ax[1,2].set(ylabel='Cumulative energy [10⁴⁸ erg]');ax[1,2].legend(fontsize=8)
    for a in ax.flat:a.set_xlabel('Time [yr]')
    fig.suptitle('QStarS · thermal evolution and conservation',fontsize=16)
    fig.savefig(p/'evolution.png');plt.close(fig)
    profiles=sorted(p.glob('profile_*.npz'))
    if profiles:
        fig,axes=plt.subplots(1,3,figsize=(14,4.5),constrained_layout=True)
        chosen=sorted(set(np.linspace(0,len(profiles)-1,min(5,len(profiles)),dtype=int)))
        for i in chosen:
            with np.load(profiles[i]) as d:
                label=f"{float(d['time_seconds'])/31557600:.4g} yr"
                axes[0].loglog(d['r_cm'],d['density_cgs'],label=label)
                axes[1].loglog(d['r_cm'],d['temperature_K'])
                axes[2].plot(d['r_cm'],d['entropy_cgs']/1e9)
        for a,title in zip(axes,['Density [g cm⁻³]','Temperature [K]','Specific entropy [10⁹ erg g⁻¹ K⁻¹]']):
            a.set(xlabel='Radius [cm]',ylabel=title,xscale='log')
        axes[0].legend(fontsize=8);fig.suptitle('Resolved radial profiles',fontsize=15)
        fig.savefig(p/'profiles.png');plt.close(fig)
    if len(profiles)>2:
        data=[dict(np.load(f)) for f in profiles]
        r=np.geomspace(min(d['r_cm'][0] for d in data),max(d['r_cm'][-1] for d in data),128)
        times=np.array([d['time_seconds']/31557600 for d in data])
        field=np.array([np.interp(np.log(r),np.log(d['r_cm']),np.log10(d['temperature_K']),left=np.nan,right=np.nan) for d in data])
        fig,ax=plt.subplots(figsize=(9,5),constrained_layout=True)
        im=ax.pcolormesh(times,r,field.T,shading='nearest',cmap='magma')
        ax.set(xlabel='Time [yr]',ylabel='Radius [cm]',yscale='log',title='Resolved envelope temperature')
        fig.colorbar(im,ax=ax,label='log₁₀ T [K]');fig.savefig(p/'temperature_spacetime.png');plt.close(fig)
    return [str(x) for x in p.glob('*.png')]
