from utils import *

# name, color, hatch style (symbols repeated to make them appear smaller in the plot)
mixing_styles = {1:('conv', 'g', '//////'), 
                2:('over', 'm', '******'), 
                3:('semi', 'c', 'xxxxxxxxx'), 
                4:('th',   'y', '+++++++'), 
                5:('rot',  'k', '\\\\\\'), 
                6:('rt',   'k', '\\\\\\'), 
                7:('min',  'k', '\\\\\\'), 
                8:('left', 'k', '\\\\\\')} 

# Since coordinates in history files are M/Mstar (q), need to calculate column depth knowing
# mass and radius of the original star.
Mcenter = 1.4 # Msun
Rcenter = 1.2e6 #cm 
def y(m,q):
    # Assuming radius to be ~cst, ok in hydrostatic phase
    # m in Msun units
    return (m-Mcenter)*Msun*(1-q)/(4*np.pi*Rcenter**2)


def get_burn_mix_matrix_from_history_fixed_lgy(hist, lgy_arr):
    # assuming history MesaData object with burn_regions and mixing_regions columns (same number)
    # will interpolate onto the column depth array lgy_arr

    # How many models, zones
    Nmod = len(hist.model_number)
    n=1
    while True:
        if not hist.in_data("burn_type_%d"%n):
            break
        n+=1
    N_zones = n-1

    if hist.in_data("N_burn_regions"):
        if max(hist.N_burn_regions) > N_zones:
            print("Warning: max number of burning (%d) regions exceeds %d"%(max(hist.N_burn_regions), N_zones))
        if max(hist.N_mix_regions) > N_zones:
            print("Warning: max number of mixing (%d) regions exceeds %d"%(max(hist.N_mix_regions), N_zones))

    # Initialize
    burn = np.zeros((len(lgy_arr), Nmod))   
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

        #smallest_mixing_zone = min(np.abs(np.diff(lgybounds)))
        #if abs(smallest_mixing_zone) < abs(lgy_arr[1]-lgy_arr[0]):
        #    print(smallest_mixing_zone, lgy_arr[1]-lgy_arr[0])
        #    print("Warning")

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
        eps = 1e-10
        for n in range(N_burn_zones - 1):
            if dlgy[n] < dlgy[n+1]:
                lgy.append(lgybounds[n+1] - eps*dlgy[n])
                lgy.append(lgybounds[n+1] + eps*dlgy[n])
            else:
                lgy.append(lgybounds[n+1] - eps*dlgy[n+1])
                lgy.append(lgybounds[n+1] + eps*dlgy[n+1])
        lgy.append(lgybounds[-1])

        # Use interp1d 'nearest' on the grid with values stored neared zone edges
        burn[:,i] = interpolate.interp1d(lgy, burn_raw, kind='nearest', bounds_error=False, fill_value=0)(lgy_arr)


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
        mix[:,i] = interpolate.interp1d(lgy, mix_raw, kind='nearest', bounds_error=False, fill_value=0)(lgy_arr)

    return burn,mix

def get_burn_mix_matrix_from_history_dynamic_lgy(hist):
    # In this case the lgy array will be different for every model. Points are taken
    # at the edges of zones so that all zones all automatically resolved (even the small ones)

    Mcenter = 1.4 # Msun
    Rcenter = 1.2e6 #cm 
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

    if hist.in_data("N_burn_regions"):
        if max(hist.N_burn_regions) > N_zones:
            print("Warning: max number of burning (%d) regions exceeds %d"%(max(hist.N_burn_regions), N_zones))
        if max(hist.N_mix_regions) > N_zones:
            print("Warning: max number of mixing (%d) regions exceeds %d"%(max(hist.N_mix_regions), N_zones))

    # Initialize
    lgy_mat_burn = np.zeros((2*N_zones, Nmod))
    lgy_mat_mix = np.copy(lgy_mat_burn)
    burn = np.copy(lgy_mat_burn)
    mix = np.copy(lgy_mat_burn)

    # Go through the file, figuring out the zones and writing into the matrices
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

        #smallest_mixing_zone = min(np.abs(np.diff(lgybounds)))
        #if abs(smallest_mixing_zone) < abs(lgy_arr[1]-lgy_arr[0]):
        #    print(smallest_mixing_zone, lgy_arr[1]-lgy_arr[0])
        #    print("Warning")

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
        eps = 1e-10
        for n in range(N_burn_zones - 1):
            if dlgy[n] < dlgy[n+1]:
                lgy.append(lgybounds[n+1] - eps*dlgy[n])
                lgy.append(lgybounds[n+1] + eps*dlgy[n])
            else:
                lgy.append(lgybounds[n+1] - eps*dlgy[n+1])
                lgy.append(lgybounds[n+1] + eps*dlgy[n+1])
        lgy.append(lgybounds[-1])

        n = len(lgy)
        #lgy_mat_burn[i, :n] = lgy  
        lgy_mat_burn[:n, i] = lgy  
        # In the end for plotting the points have to be ordered monotonically, can't just leave 0
        # So for the extra points (if this model has less than N_zones), we just give negative coordinates which appear on the plot (due to set_ylim)
        #lgy_mat_burn[i, n:] = np.arange(-1, -(2*N_zones-n+1), -1)
        lgy_mat_burn[n:, i] = np.arange(-1, -(2*N_zones-n+1), -1)

        #burn[i, :n] = burn_raw  # the non-used points stay =0
        burn[:n, i] = burn_raw  # the non-used points stay =0


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

        n = len(lgy)
        #lgy_mat_mix[i, :n] = lgy  
        #lgy_mat_mix[i, n:] = np.arange(-1, -(2*N_zones-n+1), -1)
        #mix[i, :n] = mix_raw  # the non-used points stay =0
        lgy_mat_mix[:n, i] = lgy  
        lgy_mat_mix[n:, i] = np.arange(-1, -(2*N_zones-n+1), -1)
        mix[:n, i] = mix_raw  # the non-used points stay =0


    return lgy_mat_burn,burn,lgy_mat_mix,mix


def plot_kippenhan(fig, ax, cbax, hist, lgy, xaxis="star_age", show_burn=True, show_mix=True, show_luminosity=True, show_num_zones=False, show_time_steps=False):

    ax.set_ylabel(r"log column depth (g cm$^{-2}$)")

    if xaxis == "star_age":
        ax.set_xlabel("t (s)")
        xx = hist.star_age*yr
    elif xaxis == "model_number":
        ax.set_xlabel("Model number")
        xx = hist.model_number

    # If lgy is given, the same is used for all models, both burning and mixing. we interpolate values on that array.
    if lgy is not None:
        yy = lgy
        X,Y = np.meshgrid(xx,yy)
        Yburn = Ymix = Y
        burn,mix = get_burn_mix_matrix_from_history_fixed_lgy(hist, lgy)

    # Otherwise, the column depth array varies for every model, different for burning and mixing zones
    else:
       lgy_mat_burn,burn,lgy_mat_mix,mix = get_burn_mix_matrix_from_history_dynamic_lgy(hist)
       Yburn = lgy_mat_burn
       Ymix = lgy_mat_mix
       X = np.outer(np.ones(Yburn.shape[0]), xx)

    ## Burning
    if show_burn:
        vmin,vmax = np.min(burn), np.max(burn)
        if vmin<0 and vmax>0:
            vmax = np.maximum(vmax, np.abs(vmin))
            vmin = -vmax
            # print(np.min(burn[burn>0]))
            # print(np.max(burn[burn<0]))

        # NonUniform is poorly documented. And it won't accept arrays for X and Y so can only do constant lgy array
        # im = mpl.image.NonUniformImage(ax, cmap='bwr', interpolation='nearest', extent=(xx[0],xx[-1],yy[0],yy[-1]))
        # # im.set_data(xx, yy[::-1], burn[::-1])
        # # print(X.shape,Yburn[::-1].shape,burn[::-1].shape)
        # im.set_data(X, Yburn[::-1], burn[::-1])
        # im.set_clim(vmin,vmax)
        # ax.add_image(im)
        # ax.set_xlim([xx[0],xx[-1]])
        # ax.set_ylim([yy[-1],yy[0]]) # flipped

        levels = np.linspace(vmin,vmax,50)
        # levels = np.sort(np.unique(burn[0,:]))
        # print(len(levels))
        im = ax.contourf(X, Yburn, burn, cmap='bwr', levels=levels)

        # Pcolor? Really slow and doesnt look good
        # cmap = plt.colormaps['bwr']
        # levels = mpl.ticker.MaxNLocator(nbins=20).tick_values(vmin, vmax)
        # norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        # im = ax.pcolor(X, Yburn, burn, cmap=cmap, norm=norm)

        if xaxis == 'star_age':
            # show model (x000) numbers at the top
            last_mod = hist.model_number[-1]
            N,T = [],[]
            n = 1000
            while n < last_mod:
                if n in hist.model_number:
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

        if cbax != None:
            cbar = fig.colorbar(im,cax=cbax,pad=0.5)
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
     
        for i,mi in enumerate(mix_types):
            mask = mix!=mi
            Xm = np.ma.MaskedArray(X, mask)
            Ym = np.ma.MaskedArray(Ymix, mask)
            Zm = np.ma.MaskedArray(mix, mask)

            img = ax.contourf(Xm,Ym,Zm,[mi-1,mi+1], colors='none', hatches=[mixing_styles[mi][2], None])
            # Set the hatch color https://github.com/matplotlib/matplotlib/issues/2789/#issuecomment-604599060
            for coll in img.collections:
                coll.set_edgecolor(mixing_styles[mi][1])
                coll.set_linewidth(0.05)
            
            hax = haxes[i]
            hax.set_title(mixing_styles[mi][0], pad=1)
            hax.add_patch(mpl.patches.Rectangle((0, 0), 2, 2, fill=False, hatch=mixing_styles[mi][2], edgecolor=mixing_styles[mi][1]))

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



def plot_mixing_zone_as_lines(fig, ax, hist, xaxis='star_age'):

    n=1
    while True:
        if not hist.in_data("mix_type_%d"%n):
            break
        n+=1
    N_zones = n-1

    mix_types = []

    for i in range(len(hist.model_number)):

        if xaxis=='star_age':
            x = hist.star_age[i]*yr
        elif xaxis=='model_number':
            x=hist.model_number[i]
                
        # current total star mass
        m = hist.star_mass[i]
        ybot = y(m,0) # total column depth of envelope

        ybounds = [ybot] + [y(m, hist.data('mix_qtop_%d'%(n+1))[i]) for n in range(N_zones)]
        lgybounds = np.log10([max(yi,1) for yi in ybounds])  # replace y=0(lgy=-inf) by y=1(lgy=0)
        mix_raw = [hist.data('mix_type_%d'%(n+1))[i] for n in range(N_zones)]

        for j,mj in enumerate(mix_raw):
            if mj!=0 and mj!=-1:
                ax.plot([x,x], [lgybounds[j],lgybounds[j+1]], '-', color=mixing_styles[mj][1], lw=0.5)

                if mj not in mix_types:
                    mix_types.append(mj)

    # custom legend for mixing types
    l,b,w,h = ax.get_position().bounds
    lax_left = l+0.01
    lax_bot = b+h+0.01
    lax_w = w/10
    lax_h = h/20
    laxes = []
    for i,mi in enumerate(mix_types):
        lax = fig.add_axes([lax_left, lax_bot, lax_w, lax_h])
        lax.xaxis.set_visible(False)
        lax.yaxis.set_visible(False)
        lax.set_title(mixing_styles[mi][0], pad=1)
        lax.add_patch(mpl.patches.Rectangle((0, 0), 2, 2, fill=False, hatch=mixing_styles[mi][2], edgecolor=mixing_styles[mi][1]))
        laxes.append(lax)

        # next one is moved to the right
        lax_left += lax_w

    return laxes
            

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





def plot_kipp_comp(history_file, mod_file, kipp_xaxis='star_age', t0=None, mixing_display='hatches'):

    # Kippenhan plot from history file, composition profile from mod file

    hist = mr.MesaData(history_file)
    lgy = np.linspace(10,3,500)
    
    fig = plt.figure(figsize=(13,5))
    gs = mpl.gridspec.GridSpec(1,3, width_ratios=(0.05,1,1), wspace=0.5)

    ax1 = fig.add_subplot(gs[0,1])
    cbax = fig.add_subplot(gs[0,0])

    if mixing_display == 'hatches':
        axes = plot_kippenhan(fig, ax1, cbax, hist, lgy, show_luminosity=True, show_time_steps=True, show_num_zones=True, xaxis=kipp_xaxis)
    
    elif mixing_display == 'lines':
        axes = plot_kippenhan(fig, ax1, cbax, hist, lgy, show_luminosity=True, show_time_steps=True, show_num_zones=True, xaxis=kipp_xaxis, show_mix=False)
        laxes = plot_mixing_zone_as_lines(fig, ax1, hist, xaxis=kipp_xaxis)
        axes.extend(laxes)

    cbax.yaxis.set_ticks_position('left')
    cbax.yaxis.set_label_position('left')

    ax1.set_ylim(10,3)

    if kipp_xaxis == 'star_age' and t0 != None:
        ax1.set_xlim([t0,hist.star_age[-1]*yr])

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


def test(run_dir,history_file, mod_file):
    pass




# Command line call
parser = argparse.ArgumentParser(description="Analyze convection")
parser.add_argument('-dir','--rundir', type=str, help='run directory', default=None)
parser.add_argument('-hist','--histfile', type=str, help='history file', default='histories/history_5_flash.data')
parser.add_argument('-mod','--modfile', type=str, help='mod file', default='models/ns_env_Edd.mod')
parser.add_argument('-t','--test', action='store_true', help='only run the test function')

parser.add_argument('-x','--xaxis', type=str, help='x-axis for kippenhan plot (star_age or model_number)', default='star_age')
parser.add_argument('-t0',type=float, help='minimum time in case x-axis=star_age', default=None)
parser.add_argument('-lines', action='store_true', help='show mixing regions as lines instead of hatches')

parser.add_argument('-s','--show', action='store_true', help='show plot dont save')
parser.add_argument('-o','--outfile', type=str,help='name of output file (will go in run dir)', default='kipp.pdf')

if __name__ == "__main__":
    args = parser.parse_args()

    if args.rundir == None:
        print("give run_dir")
        sys.exit()
    
    run_dir = args.rundir
    if run_dir[-1] != '/': run_dir += '/'

    history_file = run_dir + args.histfile
    if not os.path.exists(history_file):
        print("Can't find history file ", history_file) 
        if not args.test: sys.exit()

    mod_file = run_dir + args.modfile
    if not os.path.exists(mod_file):
        print("Can't find mod file ", mod_file) 
        #sys.exit()
        print("Will make kippenhan but not composition profile")
        mod_file = None
    
    if args.test:
        test(run_dir, history_file, mod_file)
        sys.exit()

    xaxis = args.xaxis
    if mod_file == None: 
        xaxis = 'model_number'

    mixing_display = 'hatches'
    if args.lines:
        mixing_display = 'lines'

    fig = plot_kipp_comp(history_file, mod_file, kipp_xaxis=xaxis, t0=args.t0, mixing_display=mixing_display)

    if not args.show:
        savefile = run_dir + args.outfile
        fig.savefig(savefile)
        print("Saved to ", savefile)
    else:
        plt.show()