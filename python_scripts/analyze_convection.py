from utils import *

def get_burn_mix_matrix_from_history(hist, lgy_arr):
    # assuming history MesaData object with burn_regions and mixing_regions columns (same number)
    # will interpolate onto the column depth array lgy_arr

    # since coordinates are M/Mstar, need to calculate column depth knowing
    # mass and radius of the original star.
    Mcenter = 1.3999999999859998 # Msun
    Rcenter = 1.2e6
    def y(m,q):
        # Assuming radius to be ~cst, ok in hydrostatic phase
        # m in Msun units
        return (m-Mcenter)*Msun*(1-q)/(4*np.pi*Rcenter**2)

    # How many models, zones
    Nmod = len(hist.model_number)
    n=1
    while True:
        if not hist.in_data("burn_type_%d"%n):
            break
        n+=1
    N_zones = n-1

    # Initialize
    burn = np.zeros((Nmod, len(lgy_arr)))   
    mix = np.copy(burn)

    # Go through the file, figuring out the zones and interpolating onto the column depth array
    for i in range(Nmod):

        # current total star mass
        m = hist.star_mass[i]
        ybot = y(m,0) # total column depth of envelope

        ## Burn data
        
        # column depth bounds and burning types
        ybounds = [ybot] + [y(m, hist.data('burn_qtop_%d'%(n+1))[i]) for n in range(N_zones)]
        burn_raw = [hist.data('burn_type_%d'%(n+1))[i] for n in range(N_zones)]

        # Remove invalid zones
        N_burn_zones = N_zones
        if -9999 in burn_raw:
            j = burn_raw.index(-9999)
            burn_raw = burn_raw[:j]
            ybounds = ybounds[:j+1]
            N_burn_zones = len(burn_raw)
        
        lgybounds = np.log10([max(yi,1) for yi in ybounds])  # replace y=0(lgy=-inf) by y=1(lgy=0)

        # store values near (eps) both edges of the zones, instead of just the edge or center value, for more accurate 'nearest' interpolation

        # example visualization
        # raw data (|:zone edges, x:zone centers)
        # |    x    |  x  |         x         |    x    |
        # new data (o:points saved in lgy)
        # o    x  o-|-oxo-|-o       x      o--|--o x    o
        # the edges of the domain are saved directly
        # the distance (-) between | and o is based on the smallest neighbor zone

        # every burn value gets duplicated
        burn_raw = np.repeat(burn_raw, 2)
        lgy = [lgybounds[0],]
        dlgy = np.diff(lgybounds)
        eps = 1e-5
        for n in range(N_burn_zones - 1):
            if dlgy[n] < dlgy[n+1]:
                lgy.append(lgybounds[n+1] - eps*dlgy[n])
                lgy.append(lgybounds[n+1] + eps*dlgy[n])
            else:
                lgy.append(lgybounds[n+1] - eps*dlgy[n+1])
                lgy.append(lgybounds[n+1] + eps*dlgy[n+1])
        lgy.append(lgybounds[-1])

        # Use interp1d 'nearest' on the grid with values stored neared zone edges
        burn[i,:] = interpolate.interp1d(lgy, burn_raw, kind='nearest', bounds_error=False, fill_value=0)(lgy_arr)


        ## Mixing data
        
        # column depth bounds and mixing types
        ybounds = [ybot] + [y(m, hist.data('mix_qtop_%d'%(n+1))[i]) for n in range(N_zones)]
        mix_raw = [hist.data('mix_type_%d'%(n+1))[i] for n in range(N_zones)]

        # Remove invalid zones
        N_mix_zones = N_zones
        if -1 in mix_raw:
            j = mix_raw.index(-1)
            mix_raw = mix_raw[:j]
            ybounds = ybounds[:j+1]
            N_mix_zones = len(mix_raw)
        
        lgybounds = np.log10([max(yi,1) for yi in ybounds])

        # There seems to be an issue when there are too many mixing zones and not enough columns to be able to report them
        # Mesa will (might?) output the last zone (mix_top_20=1.0) as some convective type, when really the last zone should always (?) be non-convective
        # mix_raw[-1] = 0

        mix_raw = np.repeat(mix_raw, 2)
        lgy = [lgybounds[0],]
        dlgy = np.diff(lgybounds)
        for n in range(N_mix_zones - 1):
            if dlgy[n] < dlgy[n+1]:
                lgy.append(lgybounds[n+1] - eps*dlgy[n])
                lgy.append(lgybounds[n+1] + eps*dlgy[n])
            else:
                lgy.append(lgybounds[n+1] - eps*dlgy[n+1])
                lgy.append(lgybounds[n+1] + eps*dlgy[n+1])
        lgy.append(lgybounds[-1])

        # Use interp1d 'nearest' on the grid with values stored neared zone edges
        mix[i,:] = interpolate.interp1d(lgy, mix_raw, kind='nearest', bounds_error=False, fill_value=0)(lgy_arr)

    return burn,mix




def plot_kipp_from_history(fig, ax, cbax, hist, lgy_arr, xaxis="star_age", show_burn=True, show_mix=True, show_luminosity=True, show_num_zones=False, show_time_steps=False):

    ax.set_ylabel(r"log column depth (g cm$^{-2}$)")
    yy = lgy_arr

    if xaxis == "star_age":
        ax.set_xlabel("t (s)")
        xx = hist.star_age*yr
    elif xaxis == "model_number":
        ax.set_xlabel("Model number")
        xx = hist.model_number

    X,Y = np.meshgrid(xx,yy)

    burn,mix = get_burn_mix_matrix_from_history(hist, lgy_arr)

    ## Burning
    if show_burn:
        vmin,vmax = np.min(burn), np.max(burn)
        if vmin<0 and vmax>0:
            vmax = np.maximum(vmax, np.abs(vmin))
            vmin = -vmax

        if xaxis == 'model_number':
           levels = np.linspace(vmin,vmax,20)
           im = ax.contourf(xx, yy, burn.T, cmap='bwr', levels=levels)

        elif xaxis == 'star_age':
            # NonUniformImage is weird about non-increasing y-data. Flipping data seems to work
            im = mpl.image.NonUniformImage(ax, cmap='bwr', interpolation='nearest', extent=(xx[0],xx[-1],yy[0],yy[-1]))
            im.set_data(xx, yy[::-1], burn.T[::-1])
            im.set_clim(vmin,vmax)
            ax.add_image(im)
            ax.set_xlim([xx[0],xx[-1]])
            ax.set_ylim([yy[-1],yy[0]]) # flipped

            # show model (x000) numbers at the top
            last_mod = hist.model_number[-1]
            N,T = [],[]
            n = 1000
            while n < last_mod:
                i = list(hist.model_number).index(n)
                N.append(n)
                T.append(xx[i])
                n+=1000
            N.append(last_mod)
            T.append(xx[-1])
            ax.tick_params(top=False)
            ax_top = ax.secondary_xaxis('top')
            ax_top.set_xticks(T)
            ax_top.set_xticklabels(("\#%d"%n for n in N))

        # cbaxes = fig.add_axes([0.95,0.1,0.03,0.8])
        if cbax != None:
            cbar = fig.colorbar(im,cax=cbax,pad=0.5)
            # cbar.set_label(r"log$_{10}\;\epsilon_{\rm{nuc}}$ (erg g$^{-1}$ s$^{-1}$)")  
            # its actually sign(eps_nuc-eps_neu) * log10(max(1, |eps_nuc))
            cbax.set_title(r"$\pm\log\left\vert\epsilon_{\rm{nuc}}-\epsilon_\nu\right\vert$"
                            "\n"
                           r"(erg g$^{-1}$ s$^{-1}$)")


    ## Mixing
    if show_mix:

        # Mixing types
        mix_types = np.unique(mix)[1:] # ignore 0 no-mixing

        # change hatch size
        mpl.rcParams['hatch.linewidth'] = 0.5

        # custom legend for mixing types
        l,b,w,h = ax.get_position().bounds
        hax_left = l+0.01
        hax_bot = b+h+0.01
        hax_w = w/10
        hax_h = h/20
        haxes = []
        for i in range(len(mix_types)):
            hax = fig.add_axes([hax_left, hax_bot, hax_w, hax_h])
            hax.xaxis.set_visible(False)
            hax.yaxis.set_visible(False)
            haxes.append(hax)
            # next one is moved to the right
            hax_left += hax_w

        mixing_types = {1:('conv', 'g', '//////'), 
                        2:('over', 'm', '******'), 
                        3:('semi', 'c', 'xxxxxxxxx'), 
                        4:('th',   'y', '+++++++'), 
                        5:('rot',  'k', '\\\\\\'), 
                        6:('rt',   'k', '\\\\\\'), 
                        7:('min',  'k', '\\\\\\'), 
                        8:('left', 'k', '\\\\\\')} 
     
        for i,mi in enumerate(mix_types):
            mask = mix.T!=mi
            Xm = np.ma.MaskedArray(X, mask)
            Ym = np.ma.MaskedArray(Y, mask)
            Zm = np.ma.MaskedArray(mix.T, mask)
            img = ax.contourf(Xm,Ym,Zm,[mi-1,mi+1], colors='none', hatches=[mixing_types[mi][2], None])
     
            for coll in img.collections:
                coll.set_edgecolor(mixing_types[mi][1])
                coll.set_linewidth(0.05)

            hax = haxes[i]
            hax.set_title(mixing_types[mi][0], pad=1)
            hax.add_patch(mpl.patches.Rectangle((0, 0), 2, 2, fill=False, hatch=mixing_types[mi][2], edgecolor=mixing_types[mi][1]))

    ax.invert_yaxis()
    axes = [ax, cbax]
    if show_mix:
        axes += haxes

    ## Luminosity, zone number, timestep curves
    if show_luminosity or show_num_zones or show_time_steps:
        ax.tick_params(right=False)
        axb = ax.twinx()
        ylabel = ""
        if show_luminosity:
            lgL_max = max(hist.log_L)
            axb.plot(xx, hist.log_L - lgL_max, 'k-', lw=0.5, alpha=0.7, label=r"$L_{\rm{max}}$ = %.2e erg s$^{-1}$"%10**(lgL_max+np.log10(Lsun)))
            ylabel += r"log $L/L_{\rm{max}}$"
        if show_num_zones:
            num_zones_max = max(hist.num_zones)
            axb.plot(xx, np.log10(hist.num_zones/num_zones_max), 'k--', lw=0.5, alpha=0.7, label=r"$N_{z,\rm{max}}$ = %d"%num_zones_max)
            if ylabel != "": ylabel += " ; "
            ylabel += r"log $N_z/N_{z,\rm{max}}$"
        if show_time_steps:
            #lgdt_max = max(hist.log_dt)
            lgdt_min = min(hist.log_dt)
            # axb.plot(xx, hist.log_dt - lgdt_max, 'k:', lw=0.5, alpha=0.7, label=r"$dt_{\rm{max}}$ = %.2e s"%10**(lgdt_max+np.log10(yr)))
            axb.plot(xx, lgdt_min -  hist.log_dt, 'k:', lw=0.5, alpha=0.7, label=r"$dt_{\rm{min}}$ = %.2e s"%10**(lgdt_min+np.log10(yr)))
            if ylabel != "": ylabel += " ; "
            ylabel += r"log $dt_{\rm{min}}/dt$"

        axb.set_ylim(top=1)
        axb.set_ylabel(ylabel)
        leg = axb.legend(loc=2,ncol=1,shadow=False,framealpha=0.5,facecolor="gray")
        for line in leg.get_lines():
            line.set_linewidth(1.0)

    return axes


# # testing kipp plot
# fig,ax = plt.subplots(1,1)
# hist = mr.MesaData("runs/G1/histories/history_5_flash.data")
# lgy = np.linspace(9,3,500)
# plot_kipp_from_history(fig,ax,cbax=None,hist=hist,lgy_arr=lgy,xaxis='model_number',show_burn=True,show_mix=False)
# plt.show()
            

def plot_composition_from_model(fig, ax, mod):

    ax.set_xlabel(r"log column depth (g cm$^{-2}$)")
    ax.set_ylabel(r"log X")

    # column depth is sum of dq/r^2
    column = integrate.cumtrapz(mod.dq/mod.R**2) * mod.xmstar/(4*np.pi)
    lgy = np.log10(column)

    linestyles = ['-',':','-.','--']
    cnt = 0
    for iso in get_isotope_list():
        if mod.in_data(iso):
            X = mod.data(iso)
            if max(X) > 1e-4:
                X[X==0] = 1e-99
                lgX = np.log10(mod.data(iso))
                ax.plot(lgy, lgX[1:], ls=linestyles[cnt//10], label=latexify_iso(iso))
                cnt+=1

    box,trans = (1.1,1),ax.transAxes
    ax.legend(frameon=False, ncol=2, bbox_to_anchor=box, bbox_transform=trans, loc='upper left')





def plot_kipp_comp(history_file, mod_file, kipp_xaxis='star_age'):

    # Kippenhan plot from history file, composition profile from mod file

    hist = mr.MesaData(history_file)
    lgy = np.linspace(9,3,500)
    
    fig = plt.figure(figsize=(13,5))
    gs = mpl.gridspec.GridSpec(1,3, width_ratios=(0.05,1,1), wspace=0.5)

    ax1 = fig.add_subplot(gs[0,1])
    cbax = fig.add_subplot(gs[0,0])
    axes = plot_kipp_from_history(fig, ax1, cbax, hist, lgy, show_luminosity=True, show_time_steps=True, show_num_zones=True, xaxis=kipp_xaxis)
    # plot_kipp_from_history(fig, ax1, cbax, hist, lgy, xaxis='model_number', show_luminosity=True, show_time_steps=True, show_num_zones=True)
    cbax.yaxis.set_ticks_position('left')
    cbax.yaxis.set_label_position('left')

    if kipp_xaxis == 'star_age':
        ax1.set_xlim([7.5,hist.star_age[-1]*yr])

    # Bring colorbar closer
    l1,b1,_,_ = ax1.get_position().bounds
    l2,b2,w2,h2 = cbax.get_position().bounds
    cbax.set_position([l2+0.5*(l1-l2), b2, w2, h2])

    if mod_file != None:
        mod = mr.MesaData(mod_file)
        ax2 = fig.add_subplot(gs[0,2])
        plot_composition_from_model(fig, ax2, mod)
        ax2.set_xlim([3,10])
        ax2.set_ylim([-4,0.01])
        ax2.set_yticks([0,-1,-2,-3,-4])
        axes += [ax2]

    # Bring all plots towards the left
    for ax in axes:
        l,b,w,h = ax.get_position().bounds
        ax.set_position([l-0.1,b,w,h])

    return fig


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("give run_dir")

    else:
        run_dir = sys.argv[1]
        if run_dir[-1] != '/': run_dir+='/'

        ## Kipp for flash, composition at end of flash (some fraction of Eddington)

        history_file = run_dir + "histories/history_5_flash.data"
        if not os.path.exists(history_file):
            sys.exit("Can't find history file ", history_file)

        # mod_file = run_dir + "models/ns_env_Edd.mod" 
        # if not os.path.exists(mod_file):
        #     sys.exit("Can't find mod file ", mod_file)

        print("Making kippenhan/composition plot for ", run_dir)

        #fig = plot_kipp_final_composition("runs/G1/histories/history_5_flash.data", "runs/G1/models/ns_env_Edd.mod")
        # fig = plot_kipp_comp(history_file, mod_file)
        fig = plot_kipp_comp(history_file, mod_file=None, kipp_xaxis='model_number')
        fig.savefig(run_dir + "kipp_comp_Edd.pdf")