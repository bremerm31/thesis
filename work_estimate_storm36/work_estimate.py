import numpy as np
from scipy.io import netcdf

import matplotlib.pyplot as plt
from matplotlib import ticker

#constants
#radius of the earth
R=6738206.4

#reference latitude and longitudes
lat0=29
long0=265.5
R0=R*np.cos(np.pi/180. * lat0)

#acceleration due to gravity
g=9.81

#CFL constant
CFL=1/3

def get_inscribed_radius(element_,x_,y_):
    element = np.array(element_[:], dtype=int)-1
    x = np.array(x_[:]) #longidute
    y = np.array(y_[:]) #latitude

    #CPP projection convert (lat,long) to meters
    x = R0 * ( x - long0) * np.pi/180.
    y = R * y * np.pi/180.

    length = np.empty(element.shape)
    for i in range(3):
        length[:,i] = np.hypot( x[element[:,(i+1)%3]] - x[element[:,i]],
                                y[element[:,(i+1)%3]] - y[element[:,i]] )

    #semi-perimeter
    s = 0.5 * np.sum(length,axis=1)

    return np.sqrt( (s-length[:,0]) * (s-length[:,1]) * (s-length[:,2])/ s)

def get_lambda(idx, element_np, zeta, depth_np, u_vel, v_vel):
    c = np.empty(element_np.shape[0])

    ze = np.array(zeta[idx,:])
    ze[ze<-888888]=np.nan
    h_node = ze + depth_np
    h_node[np.isnan(h_node)] = 0.
    h_node[h_node<0.] = 0.

    c = np.sqrt(g*np.max(h_node[element_np],axis=1))

    vel_node = np.hypot(u_vel[idx,:], v_vel[idx,:])

    return c + np.max(vel_node[element_np],axis=1)


if __name__=="__main__":
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rc('font', size=10)

    fort63 = netcdf.netcdf_file('Storm36.fort.63.nc','r')
    fort64 = netcdf.netcdf_file('Storm36.fort.64.nc','r')

    r = get_inscribed_radius(fort63.variables['element'],
                             fort63.variables['x'],
                             fort63.variables['y'])

    min_r = np.min(r)
    max_r = np.max(r)
    print('Min radius: {:f} [m]'.format(min_r))
    print('Max radius: {:f} [m]'.format(max_r))
    bins_r = min_r * np.logspace(0,11, num=20, base=2)
    plt.figure(figsize=(3.125,2.5))
    plt.hist(r, bins=bins_r)
    plt.gca().set_xscale("log")
    plt.xlabel('Inscribed radius, $h$ (in $\mathrm{m}$)')
    plt.ylabel('\# of elements')
    plt.savefig('h_distribution.pdf', bbox_inches="tight")
    plt.show()

    #compute dt
    element_np = np.array(fort63.variables['element'][:], dtype=int)-1
    depth_np = np.array(fort63.variables['depth'][:])

    n_elements = element_np.shape[0]

    #do 2 passes over the data
    #1st: compute dt_min
    '''dt_min = np.finfo(float).max

    for idx in range(fort63.variables["time"].shape[0]):
        c = get_lambda(idx, element_np,
                       fort63.variables["zeta"], depth_np,
                       fort64.variables["u-vel"],
                       fort64.variables["v-vel"])

        dt_min_loc = np.min(r[c>0]/c[c>0])
        print("  {:d}  {:f}".format(idx, dt_min_loc))

        dt_min = min(dt_min,
                     dt_min_loc)

    #round dt min down to the nearest 10th
    dt_min = np.floor(dt_min*10.)/10.
    '''
    dt_min=0.13
    print("dt_min: {:f} [s]".format(dt_min))

    work_sync = fort63.variables["time"][-1]/dt_min*n_elements

    t_series = np.array(fort63.variables["time"][:])

    #2nd pass: compute dt_lts
    work = np.empty((t_series.size, n_elements//1000+1))
    wet_frac = np.empty(t_series.shape)

    timestepping_snapshots=set([1,76])
    bins_snap = dt_min * np.array([2**ex for ex in range(int(np.ceil(np.log2(3600/dt_min)))+1)])


    work_lts = 0
    for idx,t in enumerate(t_series):
        c = get_lambda(idx, element_np,
                       fort63.variables["zeta"], depth_np,
                       fort64.variables["u-vel"],
                       fort64.variables["v-vel"])

        wet_frac[idx] = np.sum(c>0)/n_elements

        t_prev = t_series[idx-1] if idx>0 else 0

        dt_loc = CFL * r/c
        dt_loc[c<0] = t - t_prev

        dt_loc = dt_min * 2**(np.floor(np.log2(dt_loc/dt_min)))
        dt_loc = np.minimum(dt_loc, float(t - t_prev))

        work[idx,:] = (t-t_prev)/dt_loc[::1000]
        if idx==0 or idx==t_series.size-1:
            work_lts += 0.5 * np.sum((t - t_prev)/dt_loc)
        else:
            work_lts += np.sum((t-t_prev)/dt_loc)

        if idx in timestepping_snapshots:
            sanity_check,_=np.histogram(dt_loc, bins=bins_snap)
            print(" Elements in histogram: {:d} (should be {:d})".format(np.sum(sanity_check), n_elements))
            print(" {:.1f}% of elements at smallest level".format(sanity_check[0]/np.sum(sanity_check)*100))
            print(" {:.1f}% of elements are dry".format(sanity_check[-1]/np.sum(sanity_check)*100))
            print(sanity_check)

            percentiles = [25, 50, 75]
            dt_pct = np.percentile(dt_loc, percentiles, interpolation='higher')
            for jdx, pct in enumerate(percentiles):
                print( " {:.2f}%: {:f}".format(pct, dt_pct[jdx]))

            plt.figure(figsize=(2.75,2.5))
            plt.hist(dt_loc, bins=bins_snap)
            plt.gca().set_xscale("log")
            plt.gca().set_ylim([0,1.4e6])
            plt.xlabel('Timestep size (in $\mathrm{s}$)')
            plt.ylabel('\# of elements')
            plt.savefig('t_distribution'+str(int(t))+'.pdf', bbox_inches="tight")
            plt.show()

        if idx%10==0:
            print(" finished iteration: {:d}\r".format(idx))

    srtd_indices = np.argsort(np.sum(work, axis=0))

    fig, ax = plt.subplots(figsize=(6.25,3))
    tv, xv = np.meshgrid(t_series, np.arange(n_elements//1000+1))
    CS = ax.contourf(tv, xv, work[:,srtd_indices].T,)
    cbar = fig.colorbar(CS)
    cbar.ax.set_ylabel('Timesteps to simulate 1 hour')
    plt.xlabel('Simulation time (in $\mathrm{s}$)')
    plt.ylabel('Elements')
    plt.savefig('work_contour.pdf', bbox_inches="tight")
    plt.show()

    print("Synchronous updates required: {:f}".format(work_sync))
    print("LTS updates required:         {:f}".format(work_lts))
    print("Speed-up work:                {:1.2f}".format(work_sync/work_lts))

    plt.figure(figsize=(3.125,2.5))
    plt.plot(t_series, wet_frac)
    plt.xlabel('Simulation time (in $\mathrm{s}$)')
    plt.ylabel('Fraction of wet elements')
    plt.savefig("wet_fraction.pdf", bbox_inches="tight")
    plt.show()
