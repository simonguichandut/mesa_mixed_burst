# The last inlist (fallback) starts by removing the outer envelope
# which is being ejected (though some part of it may have negative velocity). This script
# finds the outermost point with acceleration 0 (or other methods) and then modifies the 
# line remove_initial_surface_by_density in the inlist

from utils import *

def remove_inversions(x,y):
    # remove inversions in x (since x needs to be monotonically increasing for spline fits)
    # 
    pass

def find_k(data, method, special=None):

    r = data.R
    rho = data.d
    T = data.T
    v = data.v
    dr = abs(data.R[1:]-data.R[:-1])
    dr = np.insert(dr, 0, dr[0]) # just duplicate first element so that array is same length as others

    ycol = -integrate.cumtrapz(rho,r,initial=0)

    # rhop = (rho[2:]-rho[:-2])/(r[2:]-r[:-2])
    # rhopp = (rho[2:]+rho[:-2]-2*rho[1:-1])/((r[2:]-r[:-2])/2)**2 
    # Tp = (T[2:]-T[:-2])/(r[2:]-r[:-2])
    # Tpp = (T[2:]+T[:-2]-2*T[1:-1])/((r[2:]-r[:-2])/2)**2 
    # vp = (v[2:]-v[:-2])/(r[2:]-r[:-2])
    # vpp = (v[2:]+v[:-2]-2*v[1:-1])/((r[2:]-r[:-2])/2)**2 


    ## Different ways to find where to cut the data

    def zero_acceleration():
        # Where dv/dr first becomes zero, starting from outside
        a = (v[2:]-v[:-2])/(r[2:]-r[:-2])
        for i,ai in enumerate(a):
            if ai<0: break
        return i+1

    def stepsize_jump(threshold=10):
        # Where the stepsize increases suddenly (measured by the instantaneous d(stepsize)/dr)

        # 10 is arbitrary, found by graphing and looking at ddr. Does not actually work all the time for varying grids
        if type(threshold) == str:
            threshold = eval(threshold)

        ddr = (dr[2:]-dr[:-2])/(r[2:]-r[:-2])
        for i,ddri in enumerate(ddr):
            if ddri>threshold:
                break
        if i==len(ddr)-1:
            print("Did not find a stepsize drop larger than the threshold. Expecting a crash on next list indexing")
        return i+1
        # return i + 25  
        # D1 didn't work with +10 but does with +20! Not going deep enough somehow breaks things. Finds L<0, then many retries, then L and Teff both have dropped a lot. Probably something in the atm calculation..

    def sonic():
        # cut at sonic point
        X,Y = data.h1, data.he4
        Z = 1-X-Y
        mu = 1/(2*X + 0.75*Y + 0.5*Z)
        cs2 = kB*T/mu/mp
        isonic = np.argmin(abs(v**2 - cs2))
        return isonic

    def power_law_rho_r(threshold=10):
        # Detect where the power law index dln(rho)/dln(r) exceeds threshold. We expect n->2 at infinity for winds
        # Need to smooth the density points to get rid of noise
        xx,yy = np.log(r[::-1]), np.log(rho[::-1])
        sp = interpolate.UnivariateSpline(xx,yy, s=0.5) # s=0.5 smoothing factor seems to work well
        n = -sp.derivative()(xx)
        for i,ni in enumerate(n[::-1]):
            if ni>threshold:
                break
        if i==len(n)-1:
            print("Did not find a power law index larger than the threshold. Expecting a crash on next list indexing")
        return i+1

    def power_law_T_rho(threshold=0.1):
        # Detect where the power law index dln(T)/dln(rho) drops below threshold. We expect n=1/3 in hydrostatic, radiation dominated regions
        xx,yy = np.log(rho), np.log(T)
        print(np.diff(xx))
        sp = interpolate.UnivariateSpline(xx,yy, s=0.5) # s=0.5 smoothing factor seems to work well
        n = sp.derivative()(xx)
        for i,ni in enumerate(n):
            if ni<threshold:
                break
        if i==len(n)-1:
            print("Did not find a power law index larger than the threshold. Expecting a crash on next list indexing")
        return i+1


    def other_method():
        pass

    methods = {func.__name__ : func for func in (zero_acceleration,stepsize_jump,sonic,power_law_rho_r,power_law_T_rho,other_method)}
    if special is None:
        return methods[method]()
    else:
        return methods[method](special)


## Examine graphically
def plot(data,k_remove,method):
    fig,axes = plt.subplots(3,4,figsize=(13,8))

    # substract_rns = False
    substract_rns = True

    for ax in axes[-1]:
        if substract_rns:
            ax.set_xlabel("r-R (cm)")
        else:
            ax.set_xlabel("r (cm)")

    # T plot as function of rho
    axes[-1][1].set_xlabel(r"$\rho$")
    
    for i,s in enumerate((r'$\rho$',r'$T$',r'$v$',r'$\Delta r$')):
        axes[0][i].set_ylabel(r"%s"%s)
        axes[1][i].set_ylabel(r"$d$%s/$dr$"%s)
        axes[2][i].set_ylabel(r"$d^2$%s/$dr^2$"%s)

    r,T,rho,v,L = data.R,data.T,data.d,data.v,data.L
    dr = abs(r[1:]-r[:-1])
    dr = np.insert(dr, 0, dr[0]) 
    
    if substract_rns:
        r -= 12e5

    # photosphere
    term1 = sigmarad*T**4
    term2 = L / (4*np.pi*r**2)
    iph = np.argmin(abs(term1/term2-1))

    # sonic point
    X,Y = data.h1, data.he4
    Z = 1-X-Y
    mu = 1/(2*X + 0.75*Y + 0.5*Z)
    cs2 = kB*T/mu/mp
    isonic = np.argmin(abs(v**2 - cs2))


    for i,var in enumerate((rho,T,v,dr)):

        x = rho if i==1 else r

        # first and second derivatives
        varp = (var[2:]-var[:-2])/(r[2:]-r[:-2])
        varpp = (var[2:]+var[:-2]-2*var[1:-1])/((x[2:]-x[:-2])/2)**2 
       
        axes[0][i].plot(x, var, 'k.')
        axes[1][i].plot(x[1:-1], varp, 'k.')
        axes[2][i].plot(x[1:-1], varpp, 'k.')

        axes[0][i].plot(x[iph], var[iph], 'go', label='photosphere')
        axes[1][i].plot(x[iph], varp[iph], 'go')
        axes[2][i].plot(x[iph], varpp[iph], 'go')

        axes[0][i].plot(x[isonic], var[isonic], 'bx', label='sonic point')
        axes[1][i].plot(x[isonic], varp[isonic], 'bx')
        axes[2][i].plot(x[isonic], varpp[isonic], 'bx')
        
        for j in range(3):
            axes[j][i].axvspan(xmin=x[k_remove], xmax=x[0], color='r', alpha=0.2)
            # axes[j][i].axvline(x[k_remove], ls='--', color='r')
            axes[j][i].set_xlim(auto=True)
            axes[j][i].set_xscale('log')

    
    axes[0][0].set_yscale('log')
    axes[0][1].set_yscale('log')

    axes[0][0].legend()

    ax_highlight = None
    if method == 'zero_acceleration':
        ax_highlight = axes[1,2]
    elif method == 'sonic':
        ax_highlight = axes[0,2]
    elif method == 'stepsize_jump':
        ax_highlight = axes[1,3]
    elif method == 'power_law_rho_r':
        ax_highlight = axes[1,0]

        # in that case also plot the power law index
        figb, axb = plt.subplots(1,1)
        axb.set_ylim([0,30])
        axb.set_xlabel('r (cm)')
        axb.set_ylabel(r'$n=d\log\rho/d\log r$')
        if substract_rns:
            r += 12e5
        xx,yy = np.log10(r[::-1]), np.log10(rho[::-1])
        n = -np.diff(yy)/np.diff(xx)
        sp = interpolate.UnivariateSpline(xx,yy, s=0.5) # s=0.5 smoothing factor seems to work well
        n_smooth = -sp.derivative()(xx)
        
        xx = 10**(xx[::-1])
        yy,n,n_smooth = 10**yy[::-1],n[::-1],n_smooth[::-1]

        axb.plot(xx[1:], n, '-', color='gray', alpha=0.5, lw=0.3, label='instantaneous derivative')
        axb.plot(xx, n_smooth, 'k-', label=r'from smoothed $\rho(r)$')
        axb.plot(xx[iph],n_smooth[iph],'go', ms=5)
        axb.plot(xx[isonic],n_smooth[isonic],'bx', ms=5)
        axb.set_xscale('log')
        axb.legend()
        # axb.axvline(xx[k_remove], ls='--', color='r', lw=0.7)
        axb.axvspan(xx[k_remove],xx[0], color='r', alpha=0.2)

    elif method == 'power_law_T_rho':
        ax_highlight = axes[1,1]

        figb, axb = plt.subplots(1,1)
        axb.set_ylim([0,30])
        axb.set_xlabel(r'$\rho$ (g cm$^{-3}$)')
        axb.set_ylabel(r'$n=d\log T/d\log\rho$')
        xx,yy = np.log10(rho), np.log10(T)
        n = np.diff(yy)/np.diff(xx)
        sp = interpolate.UnivariateSpline(xx,yy, s=0.5) # s=0.5 smoothing factor seems to work well
        n_smooth = sp.derivative()(xx)
        
        xx,yy = 10**xx, 10**yy
        axb.plot(xx[1:], n, '-', color='gray', alpha=0.5, lw=0.3, label='instantaneous derivative')
        axb.plot(xx, n_smooth, 'k-', label=r'from smoothed $\rho(r)$')
        axb.plot(xx[iph],n_smooth[iph],'go', ms=5)
        axb.plot(xx[isonic],n_smooth[isonic],'bx', ms=5)
        axb.set_xscale('log')
        axb.legend()
        # axb.axvline(xx[k_remove], ls='--', color='r', lw=0.7)
        axb.axvspan(xx[k_remove],xx[0], color='r', alpha=0.2)

    if ax_highlight is not None:
        for spine in ax_highlight.spines.values():
            spine.set_edgecolor('red')

    plt.tight_layout()
    plt.show()


def get_k_remove_info(data, k_remove):
    k = k_remove 
    ycol = -integrate.cumtrapz(data.d,data.R,initial=0)
    print(f"Removing at index {k}, r={data.R[k]/1e5:.5f} km, rho={data.d[k]:.3e} g/cm3, y={ycol[k]:.3e} g/cm2")

    m = np.cumsum(data.dq)*data.xmstar

    # another way to calculate the mass
    # m2 = -integrate.cumtrapz(data.d*4*np.pi*data.R**2, data.R, initial=0)
    # print(max(abs(m2[1:]-m[1:])/m[1:]))
    # print(m[k],m2[k])
    r1,r2 = data.R[k],data.R[0] 

    print(f'Mass being removed is {m[k]:.3e} g, out of a total {m[-1]:3e} g for the envelope (m_remove/M = {m[k]/m[-1]:.1e})')
    print(f'(estimating the mass from the first density point and assuming a r^-2 power law gives {4*np.pi*data.d[k]*r1**2*(r2-r1):.2e} g)')
    print(f'(" " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " r^-3 " " " " " " " " {4*np.pi*data.d[k]*r1**3*np.log(r2/r1):.2e} g)')
    print(f'X={data.h1[k]:.2f} there ({data.h1[0]:.2f} at the surface)\n')


    print(f'Optical depth estimates')
    tau = -integrate.cumtrapz(data.d*0.2*(1+data.h1), data.R, initial=2/3)
    print(f'Integrated density directly (opacity is scattering only) gives tau={tau[k]:.3e}')

    tau_avg_0 = 0.2*(1+data.h1[0]) * (r2-r1) * m[k]/(4/3*np.pi*(r2**3-r1**3))
    print(f'Given r2={r2/1e5:.2f} km, the approximate optical depth (kappa0*(r2-r1)*(M/V)) gives tau={tau_avg_0:.2e}')
    print(f'tau* = kappa*rho*r1 = {0.2*(1+data.h1[0])*data.d[k]*data.R[k]:.2e}')
    print(f'tau* with corrections for not going to infty : {0.2*(1+data.h1[0])*data.d[k] * data.R[k]**2*(1/data.R[k] - 1/data.R[0]):.2e}')
    print(f'tau from power law r^-3 (basically tau*/2) :   {0.2*(1+data.h1[0])*data.d[k]/2 * data.R[k]**3*(1/data.R[k]**2 - 1/data.R[0]**2):.2e}\n')

    v1,v2 = data.v[k],data.v[0]
    print(f'Given (v1,v2)=({v1:.2e} cm/s), {v2:.2e} cm/s)')
    



# Edit inlist file
def change_file(k_remove):
    filename = "inlist_folder/inlist8_fallback"
    tag = "remove_initial_surface_at_cell_k"
    new_file_contents = ""
    with open(filename,'r') as f:
        for line in f:
            if tag in line:
                new_file_contents += "   remove_initial_surface_at_cell_k = " + str(k_remove) + "\n"
            else:
                new_file_contents += line

    with open(filename,'w') as f:
        f.write(new_file_contents)



# Command line call
parser = argparse.ArgumentParser(description="find index k at which to cut the data given some method (default zero acceleration)")
parser.add_argument('-i', '--input', type=str, default='models/ns_env_ejected.mod')
parser.add_argument('-m', '--method', type=str, default='zero_acceleration')
parser.add_argument('-k','--index', type=int, default=None, help="give index k directly (probably using -p to plot)")
parser.add_argument('-p', '--plot', action='store_true', help="plot to show location of k_remove")
args = parser.parse_args()

if __name__ == "__main__":

    ## Load data
    if not os.path.exists(args.input):
        sys.exit("Could not load mod file " + args.input + " (you should be running this from a run dir)")

    data = mr.MesaData(args.input)
    print("Loaded mod file %s\n"%args.input) 

    if args.index is not None:
        print("Specified index where to remove")
        k_remove = args.index
        method = None
    else:
        print("Removing with method: ", args.method)

        if ":" not in args.method:
            method = args.method
            # k_remove = methods[method]()
            k_remove = find_k(data, method)
        else:
            method, special = args.method.split(":")
            # k_remove = methods[method](special)
            k_remove = find_k(data, method, special)

    get_k_remove_info(data, k_remove)

    if args.plot:
        plot(data,k_remove,method)
    else:
        try:
            change_file(k_remove)
        except Exception as e:
            print("could not change inlist, because:")
            print(e)