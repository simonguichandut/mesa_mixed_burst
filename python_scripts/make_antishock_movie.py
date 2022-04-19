from compare_models import *
from matplotlib.animation import FuncAnimation,FFMpegWriter
from os import listdir
from os.path import isfile, join

yr = 3.1536e7
Lsun = 3.85e33

def make_movie(log_dir,movie_filename):

    if log_dir[-1]!='/': log_dir += '/' 
    # num_profiles = len([f for f in listdir(logdir) if isfile(join(logdir, f)) and ".data" in f and "history" not in f])
    index = mr.MesaProfileIndex(log_dir+"profiles.index")
    num_profiles = len(index.profile_numbers)
    filename = lambda n: log_dir+"profile%d.data"%n
    print(f"Making movie using {num_profiles} log files in {log_dir}")

    # fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(12,8),sharex=True)
    fig,(ax1,ax3) = plt.subplots(2,1,figsize=(12,8),sharex=True)
    fig.subplots_adjust(hspace=0)
    ax3.set_xlabel(r'$r-R$ (cm)')

    ax1.set_ylabel(r'$\rho$ (g cm$^{-3}$)')
    # ax2.set_ylabel(r'$v$ (cm s$^{-1}$)')
    ax3.set_ylabel(r'X')
    ax1b = ax1.twinx()
    # ax2b = ax2.twinx()
    # ax3b = ax3.twinx()
    ax1b.set_ylabel(r'$\mu$',color='r')
    # ax2b.set_ylabel(r'$c_s$ (cm s$^{-1}$)',color='r')

    # ax1b.axhline(4/3, xmin=0.5, xmax=1, color='r', ls=':', lw=0.5)
    # ax1b.axhline(0.5, xmin=0.5, xmax=1, color='r', ls=':', lw=0.5)

    line_rho, = ax1.loglog([],[],'k-',marker='.',ms=2)
    line_mu, = ax1b.semilogx([],[],'r-', lw=0.8)
    # line_v, = ax2.loglog([],[],'k-',marker='.',ms=2) 
    # line_cs, = ax2.loglog([],[],'r-',lw=0.8) 
    # lines_basic = [line_rho,line_mu,line_v,line_cs]
    lines_basic = [line_rho,line_mu]

    linestyles = ['-',':','-.','--']#,'-','-','-']
    isotope_list = get_isotope_list()
    lines_iso = {}

    for i, prof_number in enumerate(index.profile_numbers):
        data = mr.MesaData(filename(prof_number))

        if i==0: t0=data.star_age

        # if max(data.T)>Tmax: Tmax=max(data.T)

        # composition
        for iso in isotope_list:
            if data.in_data(iso) and max(data.bulk_data[iso])>1e-3:
                if iso not in lines_iso.keys():
                    linestyle = linestyles[len(lines_iso)//10]
                    line, = ax3.loglog([],[],ls=linestyle,label=latexify_iso(iso),marker='.',ms=1)
                    lines_iso[iso] = line


    print('\n%d isotopes present'%len(lines_iso))
    print(list(lines_iso.keys()))

    box,trans = (1.1,1),ax3.transAxes
    ax3.legend(frameon=False, ncol=2, bbox_to_anchor=box, bbox_transform=trans, loc='upper left')

    def init():
        ax3.set_xlim([3e3,1000e5])
        ax1.set_ylim([5e-9,1e5])
        ax1b.set_ylim([0.5,2])
        # ax2.set_ylim([1,1e9])
        # ax2b.set_ylim([1,1e9])
        ax3.set_ylim([1e-2,1.4])
        return lines_basic + list(lines_iso.values())

    def update(frame):
        data = mr.MesaData(filename(frame))
        fig.suptitle(f"Model \#{data.model_number} --- Age = {data.star_age*yr:.2e} s --- dt = {data.time_step*yr:.1e} s")
        r = data.R_cm
        x = r - 12e5
        rho = data.Rho
        mu = data.mu
        # v = data.v
        ybasic = [rho,mu]

        for line,z in zip(lines_basic, ybasic):
            line.set_data(x, z)

        for iso in lines_iso.keys():
            lines_iso[iso].set_data(x, data.bulk_data[iso])

        # for ax in (ax1,ax1b,ax3,ax3b):
        # for ax in (ax3b,):
        #     ax.relim()
        #     ax.autoscale_view(tight=True,scalex=False,scaley=True)

        plt.tight_layout()
        return lines_basic + list(lines_iso.values())

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
        make_movie(logdir, movie_filename) 