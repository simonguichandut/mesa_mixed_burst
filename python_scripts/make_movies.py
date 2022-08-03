from utils import *
from matplotlib.animation import FuncAnimation,FFMpegWriter

# Movies from profiles

def make_nuc_movie(log_dir, movie_filename):
    # A 3-panel plot of temperature, nuclear reactions, and compositin, wrt column depth. animated using matplotlib funcanimation. The way I did this way by 
    # first figuring out how many lines are needed (and patches in the case of convective zones highlight), storing in dictionaries, initializing the plot elements,
    # and then updating their data in the update function. I'm not sure this was done correctly in the case of convective zones, the number of which 
    # can decrease... Regardless, it gives a good overview of the evolution of convection

    index = mr.MesaProfileIndex(log_dir+"profiles.index")
    num_profiles = len(index.profile_numbers)
    filename = lambda n: log_dir+"profile%d.data"%n
    print(f"Making movie using {num_profiles} log files in {log_dir}")

    def get_nuclear_reactions_list():
        data = mr.MesaData(filename(index.profile_numbers[0])) # load first profile 
        i0 = data.bulk_names.index("eps_nuc")
        iend = data.bulk_names.index("other")
        return data.bulk_names[i0+1:iend+1]

    fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(10,8),sharex=True)
    fig.subplots_adjust(hspace=0)
    ax3.set_xlabel(r'$y$ (g cm$^{-2}$)')
    ax1.set_ylabel(r'$T$ (K)')
    ax2.set_ylabel(r'$\epsilon_\mathrm{nuc}$ (erg g$^{-1}$ s$^{-1}$)')
    ax3.set_ylabel(r'X')
    ax1b = ax1.twinx()
    ax2b = ax2.twinx()
    ax3b = ax3.twinx()
    ax1b.set_ylabel(r'$L/L_\mathrm{Edd}$')
    ax2b.set_ylabel(r'$\epsilon_\mathrm{neu}$ (erg g$^{-1}$ s$^{-1}$)')
    ax3b.set_ylabel(r'$\kappa$')

    ax1b.axhline(1, xmin=0.15, xmax=1, color='r', ls='-', lw=0.5)

    line_T, = ax1.loglog([],[],'k-')
    line_L, = ax1b.loglog([],[],'k--',label=r'$L_\mathrm{tot}$', lw=0.8, alpha=0.3)
    line_Lrad, = ax1b.loglog([],[],'k:',label=r'$L_\mathrm{rad}$', lw=0.8)
    line_eps_nuc, = ax2.loglog([],[],'k-',label='total')
    line_eps_nuc_neu, = ax2b.loglog([],[],'k--', lw=0.8)
    line_kap, = ax3b.semilogy([],[],'k--', lw=0.8)
    lines_basic = [line_T,line_L,line_Lrad,line_eps_nuc,line_eps_nuc_neu,line_kap]

    # Get list of nuclear reaction categories that produce >1e-4 of the total nuclear
    # energy output (minus neutrinos) at some point, at some time
    # intialize a line for each in ax2, and save in a dict

    # Inside the same loop (so as to not re-load the data many times):
    # Get list of isotopes that have X>1e-3 at some grid point, at some time,
    # initialize a line for each in ax3, and save in a dict
    
    # Also do convective mixing types (store a list of zones in the data for each convective type)

    # Also check if (some) data surpasses default maximum values (so can be plotted on graph)
    Tmax = 2e9
    eps_max = 1e22
    eps_neu_max = 1e-10 # really don't know this a priori, will set bottom boudary a fixed # of orders of magnitude below the max

    linestyles = ['-',':','-.','--']#,'-','-','-']
    nuc_list = get_nuclear_reactions_list()
    lines_nuc = {}
    isotope_list = get_isotope_list()
    lines_iso = {}
    mixing_types = {'1':"convective",'2':"overshoot",'3':"semiconvective",'4':"thermohaline",'5':"rotational",'6':"rayleigh-taylor",'7':"minimumn mixing",'8':"anonymous mixing",'9':"leftover mixing"}
    mixing_colors = {'1':'b','2':'r','3':'g','4':'m','5':'gray','6':'y','7':'k','8':'k','9':'k'}
    patches_mixing = {}

    for i, prof_number in enumerate(index.profile_numbers):
        data = mr.MesaData(filename(prof_number))

        if max(data.T)>Tmax: Tmax=max(data.T)
        if max(data.eps_nuc)>eps_max: eps_max=max(data.eps_nuc)
        if max(data.eps_nuc_neu_total)>eps_neu_max: eps_neu_max=max(data.eps_nuc_neu_total)

        # nuclear reactions
        for nuc in nuc_list:
            if data.in_data(nuc) and max(data.bulk_data[nuc]/(data.eps_nuc+1e-5))>1e-4: # adding small number to avoid division by 0
                if nuc not in lines_nuc.keys():
                    linestyle = linestyles[len(lines_nuc)//10]
                    line, = ax2.loglog([],[],ls=linestyle,label=nuc)
                    lines_nuc[nuc] = line

        # composition
        for iso in isotope_list:
            if data.in_data(iso) and max(data.bulk_data[iso])>1e-3:
                if iso not in lines_iso.keys():
                    linestyle = linestyles[len(lines_iso)//10]
                    line, = ax3.loglog([],[],ls=linestyle,label=latexify_iso(iso))
                    lines_iso[iso] = line

        # mixing types
        mixing_zones = get_zones(data.conv_mixing_type)
        # mixing_zones.pop('0') # remove the 'no mixing' zones
        # Remove non-convective mixing ('0') and mixing types with no zones (dic value is empty list)
        mixing_zones = {key:val for key,val in mixing_zones.items() if key!='0' and len(val)>0}
        for mix in mixing_zones.keys():
            if mix not in patches_mixing.keys():
                patches_mixing[mix] = []
                new = True

            to_add = len(mixing_zones[mix]) - len(patches_mixing[mix])
            for i in range(to_add):
                if new and i==0:
                    patch = ax1.axvspan(0,0,alpha=0.4,fc=mixing_colors[mix],ec=None,label=mixing_types[mix]) 
                    new = False
                else:
                    patch = ax1.axvspan(0,0,alpha=0.4,fc=mixing_colors[mix],ec=None)

                patches_mixing[mix].append(patch)

    ## Remove mixing types that have no zones
    # patches_mixing = {key:val for key,val in patches_mixing.items() if len(val)>0}

    # Flatten dictionary into list of matplotlib patch objects into list for animation
    patches_list = [patch for mix in patches_mixing.keys() for patch in patches_mixing[mix]]

    print('\n%d nuclear reactions involved'%len(lines_nuc))
    print(list(lines_nuc.keys()))
    print('\n%d isotopes present'%len(lines_iso))
    print(list(lines_iso.keys()))
    print('\n%d mixing types present'%len(patches_mixing))
    print([f"{mixing_types[mix]} ({len(patches_mixing[mix])} patches)" for mix in patches_mixing])

    ax1b.legend(loc=4,frameon=False)

    box,trans = (1.1,1),ax1.transAxes
    ax1.legend(frameon=False, ncol=1, bbox_to_anchor=box, bbox_transform=trans, loc='upper left')
    box,trans = (1.1,1),ax2.transAxes
    ax2.legend(frameon=False, ncol=2, bbox_to_anchor=box, bbox_transform=trans, loc='upper left')
    box,trans = (1.1,1),ax3.transAxes
    ax3.legend(frameon=False, ncol=2, bbox_to_anchor=box, bbox_transform=trans, loc='upper left')


    def init():
        ax3.set_xlim([1e3,1e10])
        ax1.set_ylim([1e7,Tmax*2])
        ax1b.set_ylim([1e-4,1e1])
        ax2.set_ylim([10,eps_max*2])
        # ax2b.set_ylim([1e7,1e14])
        ax2b.set_ylim([eps_neu_max/1e3, eps_neu_max*2])
        ax3.set_ylim([1e-7,1.4])
        ax3b.set_ylim([1e-3,1])
        return lines_basic + list(lines_nuc.values()) + list(lines_iso.values()) + patches_list

    def update(frame):
        data = mr.MesaData(filename(frame))
        time = data.star_age*yr
        if time < 60:
            fig.suptitle(f"Model \#{data.model_number} --- t = {time:.3e} s --- dt = {data.time_step*yr:.1e} s")
        elif time < 3600:
            fig.suptitle(f"Model \#{data.model_number} --- t = {time/3600:.3e} h --- dt = {data.time_step*yr:.1e} s")
        elif time < day:
            fig.suptitle(f"Model \#{data.model_number} --- t = {time/day:.3e} days --- dt = {data.time_step*yr:.1e} s")
        else:
            fig.suptitle(f"Model \#{data.model_number} --- t = {time/yr:.3e} years --- dt = {data.time_step*yr:.1e} s")

        column = data.column_depth
        T = data.T 

        # L_over_LEdd = data.L_div_Ledd
        # L_over_LEdd = abs(data.L_div_Ledd)
        Ltot = data.luminosity
        Ledd = data.Ledd
        Lrad = data.Lrad
        # note: all luminosities in Lsun units
        L_over_LEdd = Ltot/Ledd
        # L_over_LEdd = abs(Ltot)/Ledd 
        Lrad_over_LEdd = Lrad/Ledd
        # Lrad_over_LEdd = abs(Lrad)/Ledd

        eps_nuc = data.eps_nuc
        eps_nuc_neu = data.eps_nuc_neu_total
        kap = data.opacity
        ybasic = [T,L_over_LEdd,Lrad_over_LEdd,eps_nuc,eps_nuc_neu,kap]

        for line,z in zip(lines_basic, ybasic):
            line.set_data(column, z)

        for nuc in lines_nuc.keys():
            lines_nuc[nuc].set_data(column, data.bulk_data[nuc])

        for iso in lines_iso.keys():
            lines_iso[iso].set_data(column, data.bulk_data[iso])

        mixing_zones = get_zones(data.conv_mixing_type)
        # mixing_zones.pop('0')
        mixing_zones = {key:val for key,val in mixing_zones.items() if key!='0' and len(val)>0}
        for mix in mixing_zones.keys():
            patch_counter = 0
            for zone in mixing_zones[mix]:
                a,b = column[zone[0]], column[zone[1]]
                patches_mixing[mix][patch_counter].set_xy([[a,0],[a,1],[b,1],[b,0],[a,0]])
                patch_counter += 1

            # the number of zones may have decreases, have to reset the remaining patches
            while patch_counter<len(patches_mixing[mix]):
                patches_mixing[mix][patch_counter].set_xy([[0,0],[0,1],[1,1],[1,0],[0,0]])
                patch_counter += 1

        # patches_list = [patch for patch in mix for mix in patches_mixing]
        patches_list = [patch for mix in patches_mixing.keys() for patch in patches_mixing[mix]]


        # for ax in (ax1,ax1b,ax2,ax3b):
        # for ax in (ax2b,):
        #     ax.relim()
        #     ax.autoscale_view(tight=True,scalex=False,scaley=True)

        plt.tight_layout()
        return lines_basic + list(lines_nuc.values()) + list(lines_iso.values()) + patches_list

        # fig.savefig("plots/%d.png"%frame, bbox_inches="tight")
        # plt.show()


    profile_list = index.profile_numbers
    anim = FuncAnimation(fig, update, frames=profile_list, init_func=init, blit=False)

    writervideo = FFMpegWriter(fps=30)
    anim.save(movie_filename, writer=writervideo)
    print("\nSaved to ",movie_filename)



def make_conv_movie(log_dir, movie_filename):
    # A closer inspection of convection zones. 

    g = G*1.4*Msun/12e5**2
    def H(T,mu):
        return kB*T/(mu*mp*g)

    index = mr.MesaProfileIndex(log_dir+"profiles.index")
    num_profiles = len(index.profile_numbers)
    filename = lambda n: log_dir+"profile%d.data"%n
    print(f"Making movie using {num_profiles} log files in {log_dir}")

    fig, ax = plt.subplots(1,1)
    axb = ax.twinx()
    ax.set_ylabel(r"log $\Delta r/H$")
    axb.set_ylabel(r"log $\epsilon_{\rm{nuc}}$")
    ax.set_xlabel(r"log $y$")
    ax.set_xlim(4,10)
    ax.set_ylim(-5,1.5)
    axb.set_ylim(14,20)

    colors = ['k','g','m','c','y']

    # Out of bounds points just for the legend
    ax.plot(0,0,'k.',label='gap')
    ax.plot(0,0,'g.',label='conv')
    ax.plot(0,0,'m.',label='over')
    ax.plot(0,0,'c.',label='semi')
    ax.plot(0,0,'y.',label='th')
    ax.legend(frameon=False,loc=3)

    def func_plot(prof_number):

        # Start by removing lines from previous plot
        # print(func_plot.lines)
        for line in func_plot.lines:
            line.remove()
        func_plot.lines = []

        # func_plot.lines = []

        data = mr.MesaData(filename(prof_number))

        # Burning data
        lgeps = [np.log10(e) if e>0 else 0 for e in data.eps_nuc]
        line = axb.plot(np.log10(data.column_depth), lgeps, 'r-', lw=0.5)
        func_plot.lines.append(line[0])

        # Convection zones
        mx = data.conv_mixing_type[::-1]
        r = data.R_cm[::-1]
        y = data.column_depth[::-1]
        T = data.T[::-1]
        mu = data.mu[::-1]

        mix_type = []
        mix_size = []
        mix_ytop = []
        cur_type = mx[0]
        ybot = y[0]
        rbot = r[0]
        Hbot = H(T[0],mu[0])
                
        for j in range(1,len(mx)-1):


            if mx[j] != cur_type: # change of mixing type

                mix_type.append(cur_type)

                ytop = y[j]
                mix_ytop.append(ytop)

                rtop = r[j]
                mix_size.append((rtop-rbot)/Hbot)

                xx1, xx2 = np.log10([ybot,ytop])
                yy = np.log10(mix_size[-1])
                col = colors[cur_type]
                marker = None if yy>-1 else 'o'  # cant see the lines if they are too short
                line = ax.plot([xx1,xx2], [yy,yy], color=col, ls='-', marker=marker, ms=2)
                func_plot.lines.append(line[0])

                cur_type = mx[j]
                rbot = rtop
                ybot = ytop
                Hbot = H(T[j],mu[j])

        Nzones = len([1 for x in mix_type if x!=0])
        fig.suptitle(f"Model \#{data.model_number} --- t = {data.star_age*yr:.1e} s --- log dt/s = {np.log10(data.time_step*yr):.1f} --- Nzones = {Nzones}", y=0.95)



    func_plot.lines = []

    profile_list = index.profile_numbers
    anim = FuncAnimation(fig, func_plot, frames=profile_list)
    writervideo = FFMpegWriter(fps=30)
    anim.save(movie_filename, writer=writervideo)
    print("\nSaved to ",movie_filename)



def make_wind_movie(log_dir, movie_filename):

    index = mr.MesaProfileIndex(log_dir+"profiles.index")
    filename = lambda n: log_dir+"profile%d.data"%n
    
    # profile_list = index.profile_numbers
    profile_list = []
    for n in index.profile_numbers:
        if os.path.exists(filename(n)):
            profile_list.append(n)
    num_profiles = len(profile_list)

    print(f"Making movie using {num_profiles} log files in {log_dir}")

    # Check if there is burning
    prof1 = mr.MesaData(filename(1))
    burning = False
    if max(prof1.eps_nuc)>0:
        burning = True
    
    # Panels: rho,T,L,mu,v. if burning: eps_nuc
    if burning:
        fig,axes = plt.subplots(6,1,figsize=(6,15),sharex=True)
    else:
        fig,axes = plt.subplots(5,1,figsize=(6,15),sharex=True) 
    fig.subplots_adjust(hspace=0.08)

    axes[-1].set_xlabel(r"$r$ (km)")
    axes[0].set_ylabel(r"$T$ (K)")
    axes[1].set_ylabel(r"$\rho$ (g cm$^{-3}$)")
    axes[2].set_ylabel(r"$L$ (erg s$^{-1}$)")
    axes[3].set_ylabel(r"$v$ (cm s$^{-1}$)")
    axes[4].set_ylabel(r"$\mu$")
    if burning: axes[5].set_ylabel(r"$\epsilon_{\rm nuc}$ (erg g$^{-1}$ s$^{-1}$)")

    def init_func():
        artists=[]
        for ax in axes:
            line, = ax.plot([],[],'k-')
            rph_pt, = ax.plot([],[],'k.',ms=7,mec='b')
            rs_pt, = ax.plot([],[],'kx',ms=7,mec='r')
            artists.extend((line,rph_pt,rs_pt))

            ax.set_xscale('log')
            ax.set_xlim([10,5e2])
            ax.set_yscale('log')
        
        axes[0].set_ylim([1e6,1e9])
        axes[1].set_ylim([1e-8,1e0])
        axes[2].set_ylim([1e37,1e39])
        axes[3].set_ylim([1e5,3e9])
        axes[4].set_yscale('linear')
        axes[4].set_ylim([0.5,2.1])
        axes[5].set_ylim([1e5,1e20])

        return artists

    artists=init_func()

    def func_plot(prof_number):
        prof = mr.MesaData(filename(prof_number))
        fig.suptitle(f"Model \#{prof.model_number} --- t-t$_0$ = {prof.star_age*yr:.1e} s --- dt = {prof.time_step*yr:.1e} s", y=0.9)

        r,T,rho,L,v,mu,eps = prof.R_cm,prof.T,prof.Rho,prof.luminosity,prof.velocity,prof.mu,prof.eps_nuc
        r /= 1e5
        L *= Lsun

        cs = np.sqrt(kB*T/mu/mp)
        isonic = np.argmin(abs(v-cs))
        iphot = np.argmin(abs(L/(4*np.pi*(r*1e5)**2) - sigmarad*T**4))

        ydata = [T,rho,L,v,mu]
        if burning: ydata.append(eps)

        for i,y in enumerate(ydata):
            artists[3*i].set_data(r, y)
            artists[3*i+1].set_data(r[iphot], y[iphot])
            artists[3*i+2].set_data(r[isonic], y[isonic])

        return artists


    anim = FuncAnimation(fig, func_plot, init_func=init_func, frames=profile_list, blit=True)
    writervideo = FFMpegWriter(fps=30)
    anim.save(movie_filename, writer=writervideo)
    print("\nSaved to ",movie_filename)
    

    

# Command line call
parser = argparse.ArgumentParser(description="Analyze convection")
parser.add_argument('-dir','--rundir', type=str, help='run directory', default=None)
parser.add_argument('-L','--logdir', type=str, help='log directory (withing run_dir/LOGS/', default=None)
parser.add_argument('-m','--movie', type=str, help='movie function (nuc, conv, wind)', default=None)
parser.add_argument('-o','--outfile', type=str,help='name of output file (will go in run_dir/movies)', default='movie.mp4')


if __name__ == "__main__":
    args = parser.parse_args()

    if args.rundir == None:
        print("give run_dir")
        sys.exit()

    if args.logdir == None:
        print("give log_dir")
        sys.exit()

    run_dir = args.rundir
    if run_dir[-1] != '/': run_dir += '/'
    log_dir = args.logdir
    if log_dir[-1] != '/': log_dir += '/'

    full_log_dir = run_dir + 'LOGS/' + log_dir
    movie_filename = run_dir + 'movies/' + args.outfile

    if args.movie in ('nuc','nuc_movie','nuclear'):
        make_nuc_movie(full_log_dir, movie_filename)

    elif args.movie in ('conv','conv_movie','convection'):
        make_conv_movie(full_log_dir, movie_filename)

    elif args.movie in ('wind','wind_movie'):
        make_wind_movie(full_log_dir, movie_filename)

    else:
        print("Unknown movie function ", args.movie)
        parser.print_help()

