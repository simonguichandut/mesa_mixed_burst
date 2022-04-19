from compare_models import *
from matplotlib.animation import FuncAnimation,FFMpegWriter
from os import listdir
from os.path import isfile, join

yr = 3.1536e7
Lsun = 3.85e33

def get_zones(L):
    # given a list L that looks like [2,2,5,5,5,1,3,3,3,6,8,1,1]
    # return all zones of list indices which correspond to the same number
    # return as a dictionary where the key is the number which repeats, the value is 
    # a list with all zones as (beginning,end). Empty list if no zones for that number
    # For this example, return {'2':[(0,1)]}
    dic = {str(i):[] for i in np.unique(L)}
    for i in range(1,len(L)):
        if L[i] == L[i-1]:
            s = str(L[i])
            if len(dic[s])>=1 and dic[s][-1][1] == i-1: # same zone as before so just update final index
                dic[s][-1][1] = i
            else: # create new zone
                dic[s].append([i-1,i])
    return dic

def make_nuc_movie(log_dir,movie_filename):

    if log_dir[-1]!='/': log_dir += '/' 
    # num_profiles = len([f for f in listdir(logdir) if isfile(join(logdir, f)) and ".data" in f and "history" not in f])
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

        if i==0: t0=data.star_age

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
        fig.suptitle(f"Model \#{data.model_number} --- Age = {data.star_age:.2e} yr --- dt = {data.time_step*yr:.1e} s --- t-t$_0$ = {(data.star_age-t0)*yr:.2e} s")
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


if __name__ == "__main__":
    if len(sys.argv):
        logdir, movie_filename = sys.argv[1:]
        make_nuc_movie(logdir, movie_filename) 