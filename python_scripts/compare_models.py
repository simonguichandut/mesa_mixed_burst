from utils import *

def make_figure(models, filename=None):

    fig,ax = plt.subplots(1,1)
    #fig,ax = plt.subplots(1,1,figsize=(6,4))
    #fig,ax = plt.subplots(1,1,figsize=(6.9, 4.3))
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

        #ax.legend(loc=4,frameon=False)
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

        #ax2.legend(loc=2,frameon=False)
        box,trans = (0,1),ax2.transAxes
        ax2.legend(frameon=False, ncol=1, bbox_to_anchor=box, bbox_transform=trans, loc='lower left')

    if filename is None:
        plt.tight_layout()
        plt.show()
    else:
        plt.savefig(filename,bbox_inches='tight')
        print('saved to ',filename)

#make_figure(models = ('run_5e34_1e-9/ns_Edd.mod','run_5e34_1e-9/ns_Edd_AC.mod'),filename='pdf/Edd_models.pdf')

if __name__ == "__main__":
    if len(sys.argv)>=2:
        if '.mod' not in sys.argv[-1]: # means the last arg is the filename
            make_figure(models=sys.argv[1:-1],filename=sys.argv[-1])
        else:
            make_figure(models=sys.argv[1:])