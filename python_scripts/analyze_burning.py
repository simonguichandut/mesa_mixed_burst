from utils import *

def plot_nuc_categories_from_history(history_file, xaxis="star_age", t0=None):

    hist = mr.MesaData(history_file)

    fig,ax = plt.subplots(1,1)

    ax.set_ylabel(r"$\epsilon_{\rm{nuc}}$ (erg g$^{-1}$ s$^{-1}$)")
    if xaxis == "star_age":
        ax.set_xlabel("t (s)")
        xx = hist.star_age*yr
    elif xaxis == "model_number":
        ax.set_xlabel("Model number")
        xx = hist.model_number

    pp = 10**hist.pp * Lsun
    cno = 10**hist.cno * Lsun
    tria = 10**hist.tri_alfa * Lsun

    line_zero = np.zeros(len(pp))
    line_tria = np.log10(tria)
    line_pp = np.log10(tria+pp)
    line_cno = np.log10(tria+pp+cno)
    

    ax.fill_between(x=xx, y1=line_zero, y2=line_tria, color='gray', alpha=0.7)
    ax.fill_between(x=xx, y1=line_tria, y2=line_cno, color='b', alpha=0.7)
    ax.fill_between(x=xx, y1=line_cno, y2=line_pp, color='r', alpha=0.7)


    return fig
            


# Command line call
parser = argparse.ArgumentParser(description="Analyze convection")
parser.add_argument('-dir','--rundir', type=str, help='run directory', default=None)
parser.add_argument('-hist','--histfile', type=str, help='history file', default='histories/history_5_flash.data')

parser.add_argument('-x','--xaxis', type=str, help='x-axis for kippenhan plot (star_age or model_number)', default='star_age')
parser.add_argument('-t0',type=float, help='minimum time in case x-axis=star_age', default=None)

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

    fig = plot_nuc_categories_from_history(history_file, xaxis=args.xaxis, t0=args.t0)

    if not args.show:
        savefile = run_dir + args.outfile
        fig.savefig(savefile)
        print("Saved to ", savefile)
    else:
        plt.show()