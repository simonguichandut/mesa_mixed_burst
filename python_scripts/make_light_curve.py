import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation,FFMpegWriter
import py_mesa_reader as mr
import argparse
import pickle

yr = 3.1536e7
Lsun = 3.85e33
G = 6.6726e-8
kB = 1.380658e-16
c = 2.99792458e10
mp = 1.67e-24

def make_lightcurve(log_dirs,fig_filename,zoom=None):

    fig,ax = plt.subplots(1,1,figsize=(5,4))
    ax.set_xlabel('t (s)')
    ax.set_ylabel('L (1e38 erg s/s)')

    LEdd_He = 4*np.pi*G*1.4*2e33*c/0.2 
    LEdd_H = LEdd_He/2
    LEdd_solar = LEdd_He/1.7
    ax.axhline(LEdd_He/1e38, ls='--', color='k', lw=0.7)
    ax.axhline(LEdd_solar/1e38, ls='--', color='k', lw=0.7)
    ax.axhline(LEdd_H/1e38, ls='--', color='k', lw=0.7)

    t0 = 0

    for k,log_dir in enumerate(log_dirs):

        try:

            if log_dir[-1]!='/': log_dir += '/' 
            index = mr.MesaProfileIndex(log_dir+"profiles.index")

            if index.profile_numbers.ndim==0: # only 1 profile, needs reformating
                index.profile_numbers = [index.profile_numbers.item(0), ]
                
            num_profiles = len(index.profile_numbers)
            filename = lambda n: log_dir+"profile%d.data"%n
            if not args.quiet:
                print(f"Making light curve using {num_profiles} log files in {log_dir}")

            t,Louter,Lphot,LEdd = [],[],[],[]
            for i, prof_number in enumerate(index.profile_numbers):
                
                if not args.quiet:
                    print(i,end='\r')
                
                data = mr.MesaData(filename(prof_number))
                t.append(data.star_age*yr)  
                Louter.append(data.luminosity[0]*Lsun) 

                iphot = np.argmin(np.abs(data.tau - 1))
                Lphot.append(data.luminosity[iphot]*Lsun)

                LEdd.append(LEdd_He*0.2/data.opacity[iphot])

            t,Louter,Lphot,LEdd = np.array((t,Louter,Lphot,LEdd))
            t += t0
            ax.axvline(t0, ls=':', lw=0.7)

            if k==0:
                ax.plot(t,Lphot/1e38,'k.-',ms=2.5,label='L(tau=1)')
                ax.plot(t,Louter/1e38,'r-',ms=2.5,label='L(r=rmax)')
                ax.plot(t,LEdd/1e38,'b-',label='Ledd',lw=0.6)
            else:
                ax.plot(t,Lphot/1e38,'k.-',ms=2.5)
                ax.plot(t,Louter/1e38,'r-',ms=2.5)
                ax.plot(t,LEdd/1e38,'b-',lw=0.6)

            t0 = t[-1]

        except Exception as e:
            print("Could not load profiles for ", log_dir)
            print(e)

        
    ax.legend(loc=1,framealpha=0.5)
    # ax.text(t[-1],LEdd_He/1e38,r'$L_{\rm Edd,He}(X=0)$',ha='right',va='top')
    # ax.text(t[-1],LEdd_H/1e38,r'$L_{\rm Edd,H}(X=1)$',ha='right',va='bottom')
    # ax.text(t[-1],LEdd_He/1e38,r'$L_{\rm Edd,He}(X=0)$',ha='right',va='top')

    # maximum time
    # ax.set_xlim([-5,65])
    # ax.text(64,LEdd_H/1e38,r'$L_{\rm Edd,H}(X=1)$',ha='right',va='bottom')
    # ax.text(64,LEdd_solar/1e38,r'$L_{\rm Edd,solar}(X=0.7)$',ha='right',va='bottom')
    # ax.text(64,LEdd_He/1e38,r'$L_{\rm Edd,He}(X=0)$',ha='right',va='bottom')
    
    # Set time limits
    if zoom == None:
        ax.set_xlim([0,30]) # default
    else:
        ax.set_xlim(zoom)

    fig.savefig(fig_filename, bbox_inches='tight')
    print("\nSaved to ",fig_filename)

    # Save fig handle with pickle. Serves as a lightweight portable interactive plot
    with open("lightcurve.pickle","wb") as f:
        pickle.dump(fig, f)
    

# Command line call
parser = argparse.ArgumentParser(description="Make a lightcurve")
parser.add_argument('-L','--logdirs', type=str, nargs='+', help="Directories containing profiles")
parser.add_argument('-F','--filename', type=str, help="Filename to save to")
parser.add_argument('-q','--quiet', action='store_true', help="Don't print stuff")
parser.add_argument('-z','--zoom', type=float, nargs=2, help="xlims in seconds (2 arguments)", default=None)

if __name__ == "__main__":
    args = parser.parse_args()
    make_lightcurve(args.logdirs, args.filename, args.zoom) 