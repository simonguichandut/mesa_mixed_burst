from sympy import latex
from utils import *

def make_Trho_plot(models, filename=None):

    fig,ax = plt.subplots(1,1)
    ax.set_xlabel(r'$\rho$ (g cm$^{-3}$)')
    ax.set_ylabel(r'$T$ (K)')
    ax2 = ax.twinx()
    ax2.set_ylabel(r'$L$ (erg cm$^{-2}$ s$^{-2}$)',color='r')

    iso_colors = {'none':'r'} # don't use red for composition
    linestyles = [':','--','-','-.']
    if len(models)>len(linestyles):
        print("not enough linestyles")

    for m,model in enumerate(models):
        data = mr.MesaData(model)
        most_abundant = get_most_abundant_isotope(data)

        # Plot composition in rho,T plot
        current_iso = most_abundant[0]
        iprev=0
        ranges=[]
        for i in range(len(data.d)):
            if most_abundant[i]!=current_iso:
                ranges.append([current_iso,iprev,i])
                iprev=i
                current_iso=most_abundant[i]
        ranges.append([current_iso,iprev,-1])
            
        for x in ranges:
            iso,a,b = x
            if iso not in iso_colors.keys():
                line = ax.loglog(data.d[a:b],data.T[a:b],lw=3,label=latexify_iso(iso))
                iso_colors[iso] = line[0].get_color()
            else:
                ax.loglog(data.d[a:b],data.T[a:b],lw=3,color=iso_colors[iso])

        # Legend on the top right for isotopes
        box,trans = (1,1),ax.transAxes
        ax.legend(frameon=False, ncol=2, bbox_to_anchor=box, bbox_transform=trans, loc='lower right')

        # Overlay thin black line (different linestyles) over the composition to indicate which model
        ax.loglog(data.d,data.T,lw=1.5,color='k',ls=linestyles[m])

        # Plot the luminosity with the same linestyles, add the labels including 
        # ymax and rmax
        ymax = -integrate.trapz(data.d,data.R)
        name = model.replace('/',' : ').replace('_','\_')
        lab = ('%s\n'r'$r_\mathrm{max}$=%.3e km \quad $\log y_\mathrm{max}$=%.1f'%(name,data.R[0]/1e5,np.log10(ymax)))
        ax2.loglog(data.d,data.L,color='r',ls=linestyles[m],label=lab)

        # Legend on the top left for models
        box,trans = (0,1),ax2.transAxes
        ax2.legend(frameon=False, ncol=1, bbox_to_anchor=box, bbox_transform=trans, loc='lower left')

    if filename is None:
        plt.tight_layout()
        plt.show()
    else:
        plt.savefig(filename,bbox_inches='tight')
        print('\nSaved to ',filename)


def make_Xy_plot(models, isotopes, filename=None):

    fig,ax = plt.subplots(1,1)
    ax.set_xlabel(r'$y$ (g cm$^{-2}$)')
    ax.set_ylabel(r'$X$')

    ax.set_xlim([1e3,1e10])
    ax.set_ylim([1e-4,1.2])

    # Colors for models
    # mod_colors = {}

    # Linestyles for isotopes
    linestyles = ['-',':','--','-.']
    if len(isotopes) > len(linestyles):
        print("WARNING : not enough linestyles (%d) for number of isotopes (%d)"%(len(linestyles),len(isotopes)))

    # # Black line outside of the plot for isotopes legend (top right)
    ax2=ax.twinx()
    ax2.yaxis.set_visible(False)
    for i,iso in enumerate(isotopes):
        line = ax2.plot([1e-10,1e-10],[1e-10,1e-10],color='k',ls=linestyles[i],label=latexify_iso(iso))
    box,trans = (1,1),ax2.transAxes
    leg2 = ax2.legend(frameon=False, ncol=2, bbox_to_anchor=box, bbox_transform=trans, loc='lower right')


    for m,model in enumerate(models):
        data = mr.MesaData(model)
        column = integrate.cumtrapz(data.dq/data.R**2) * data.xmstar/(4*np.pi)

        for i,iso in enumerate(isotopes):
            if i==0:
                ymax = column[-1]
                name = model.replace('/',' : ').replace('_','\_')
                lab = ('%s\n'r'$r_\mathrm{max}$=%.3e km \quad $\log y_\mathrm{max}$=%.1f'%(name,data.R[0]/1e5,np.log10(ymax)))
                line = ax.loglog(column, data.bulk_data[iso][1:], ls=linestyles[i], lw=2, label=lab)
                color = line[0].get_color()
            else:
                ax.loglog(column, data.bulk_data[iso][1:], ls=linestyles[i], color=color, lw=2)


    # Legend on the top left for models    
    box,trans = (0,1),ax.transAxes
    ax.legend(frameon=False, ncol=1, bbox_to_anchor=box, bbox_transform=trans, loc='lower left')

    if filename is None:
        plt.tight_layout()
        plt.show()
    else:
        plt.savefig(filename,bbox_inches='tight')
        print('\nSaved to ',filename)




# Command line call
parser = argparse.ArgumentParser(description="Compare MESA .mod files")
parser.add_argument('-f','--files', type=str, nargs='+', help='paths of mod files to be plotted', default=None)
parser.add_argument('-p','--plot', type=str, help='plot function (Trho, Xy)', default='Trho')
parser.add_argument('-s','--show', action='store_true', help="plt.show, don't save")
parser.add_argument('-o','--outfile', type=str,help='path and filename of output file', default='plot.png')

parser.add_argument('-iso','--isotopes', type=str, nargs='+', help='list of specific isotopes to include in Xy plot', default=['h1'])


if __name__ == "__main__":
    args = parser.parse_args()

    if args.files == None:
        print("No mod files given")
        sys.exit()
    elif True in [".mod" not in file for file in args.files]:
        print("Expecting all .mod files")
        sys.exit()

    filename = args.outfile
    if args.show:
        filename = None

    if args.plot in ('Trho','T_rho'):
        make_Trho_plot(args.files, filename)
    
    elif args.plot in ('X','Xy','composition'):
        make_Xy_plot(args.files, args.isotopes, filename)

    else:
        print("Unknown plot function", args.plot)
        parser.print_help()