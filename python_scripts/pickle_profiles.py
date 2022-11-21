from utils import *

# Point of this script is to save relevant profiles data into binary, so as to not have to save every profile locally.

# To indicate progress in loading the profiles
## https://stackoverflow.com/a/34482761
def progressbar(it, prefix="", size=60, out=sys.stdout): # Python3.6+
    count = len(it)
    def show(j):
        x = int(size*j/count)
        print(f"{prefix}[{u'█'*x}{('.'*(size-x))}] {j}/{count}", end='\r', file=out, flush=True)
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    print("\n", flush=True, file=out)


class Extended_Profile:

    def __init__(self, prof):
        self.prof = prof
        self.iph = self.iphot()

    def iphot(self):
        term1 = sigmarad*self.prof.T**4
        term2 = self.prof.luminosity*Lsun / (4*np.pi*self.prof.R_cm**2)
        ratio = term1/term2
        if min(ratio)<1 and max(ratio)>1: # the terms cross so there is a photosphere
            return np.argmin(abs(ratio-1))
        else:
            return None

    def Lph(self): # Luminosity at photosphere (Lsun)
        if self.iph is not None:
            return self.prof.luminosity[self.iph]
        else:
            return None

    def rph(self): # Photospheric radius (cm)
        if self.iph is not None:
            return self.prof.R_cm[self.iph]
        else:
            return None

    def tauph(self): # Optical depth at photosphere

        if self.iph is not None:
            return self.prof.tau[self.iph]
        else:
            return None

    def yitf(self): # highest column at which H1 is the dominant species
        most = get_most_abundant_isotope(self.prof)
        i_itf = len(most) - most[::-1].index("h1") - 1   # .index() returns the first instance, so we search in the reversed array
        return self.prof.column_depth[i_itf]

    def yh1(self): # highest column where H1 is still its constant (accretion) value
        h1 = self.prof.h1
        ih1 = list(h1).index(h1[h1<max(h1)-0.01][0])
        return self.prof.column_depth[ih1]

    def Egen(self): # energy generated (nuclear - neutrinos) (erg/g/s)
        return self.prof.eps_nuc - self.prof.non_nuc_neu

    def Mdot(self): # mass-loss rate (g/s)
        return 4*np.pi*self.prof.R_cm**2 * self.prof.Rho * self.prof.velocity

    def Mh1(self): # total column of hydrogen
        return np.sum(self.dm*self.h1)

    def Mhe4(self): # total column of heliuem
        return np.sum(self.dm*self.he4)

    def Mc12(self): # total column of carbon
        return np.sum(self.dm*self.c12)


    def _is_available(self, var):
        # Either in the profile itself or a method in this class
        return self.prof.in_header(var) or self.prof._any_version(var) or hasattr(self, var)

    def get(self, var):
        if hasattr(self, var):
            return getattr(self, var)()
        else:
            return self.prof.__getattr__(var) # gets called if var is not a method name of this class


def make_data_dict(log_dir, vars):

    index = mr.MesaProfileIndex(log_dir+'profiles.index')
    profname = lambda n: log_dir+"profile%d.data"%n

    dic = {}

    if 'model_number' not in vars:
        dic = {'model_number':[]} # included by default

    # Start getting the keys and intialize values as empty lists
    prof1 = mr.MesaData(profname(index.profile_numbers[0]))
    P1 = Extended_Profile(prof1)

    for var in vars:
        if P1._is_available(var):
            dic[var] = []
            print(var)
        else:
            print("\n %s not available either in profile or method. Will be ignored.\n"%var)

    for n in progressbar(index.profile_numbers, prefix="Reading profiles... "):
        if os.path.exists(profname(n)):
            prof = mr.MesaData(profname(n))
            P = Extended_Profile(prof)

            for key in dic.keys():
                    dic[key].append(P.get(key))

    return dic



parser = argparse.ArgumentParser(description = "Export profile data to binary")

parser.add_argument('-dir','--rundir', type=str, help='run directory', default=None)
parser.add_argument('-L','--logdir', type=str, help='log directory (withing run_dir/LOGS/', default=None)
parser.add_argument('-o','--outfile', type=str, help='name of output file (will go in run_dir/pickle)', default='data.pickle')

# Options for what to save
parser.add_argument('-vars','--variables', type=str, help='variables (comma-separated) which are included in the profile MesaData object. Saves the whole array')


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

    pickle_dir = run_dir + "pickle/"
    if not os.path.exists(pickle_dir):
        os.mkdir(pickle_dir)

    filename = pickle_dir + args.outfile

    # D = make_data_dict(full_log_dir, args.vars.split(','), args.eval.split(','))
    D = make_data_dict(full_log_dir, args.variables.split(','))

    # Pickle dict
    with open(filename, 'wb') as f:
        pickle.dump(D, f)

    print("Saved to ", filename)
