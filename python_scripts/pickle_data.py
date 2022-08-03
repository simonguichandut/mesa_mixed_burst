from utils import *

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


def pickle_burn_mix_data_from_profiles(log_dir, filename=None):

    index = mr.MesaProfileIndex(log_dir+'profiles.index')
    t,lgy,burn,mix = [],[],[],[]

    for i in progressbar(index.profile_numbers, prefix="Reading profiles... "):
        prof = mr.MesaData(log_dir + 'profile%d.data'%i)
        t.append(prof.star_age*yr)
        lgy.append(np.log10(prof.column_depth))
        burn.append(prof.eps_nuc - prof.eps_nuc_neu_total)
        mix.append(prof.conv_mixing_type)

    if filename == None:
        filename = log_dir+"burn_mix_data.pickle"

    with open(filename, 'wb') as f:
        pickle.dump([t,lgy,burn,mix], f)

    print("Saved to ", filename)


if __name__ == "__main__":
    log_dir = sys.argv[-1]
    if log_dir[-1] != '/': log_dir+='/'
    pickle_burn_mix_data_from_profiles(log_dir)
