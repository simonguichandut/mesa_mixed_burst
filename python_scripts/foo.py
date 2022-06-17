from utils import *

# this is the profile for the first model (280) that exists for this run 
# right after hitting the hydrogen region. It's at this moment (or before 
# but I dont have profiles) where a bunch of convection zones appear. 
# lets have a close look at temperature, mu, .. , in and outside the convection
# zones at that moment
prof = mr.MesaData("runs/I8c/LOGS/5_flash/profile15.data")

fig,ax=plt.subplots(1,1)

lgy = np.log10(prof.column_depth)
ax.set_xlabel('log column depth')

ax.plot(lgy,prof.T/1e8,'k-',label=r'$T/10^8$K')
ax.plot(lgy,prof.mu,'b-',label=r'$\mu$')
ax.plot(lgy,1+(prof.luminosity-prof.Lrad)*Lsun/1e39,'r-',lw=0.5,alpha=0.5,label=r'$1+(L-L_{\rm rad})/10^{39}$')

colors = ['k','g','m','c','y']
cur_type = prof.conv_mixing_type[0]
ktop = 0

for k in range(1,len(lgy)):
    if prof.conv_mixing_type[k] != cur_type:
        if cur_type > 0:
            ax.axvspan(lgy[k-1], lgy[ktop], color=colors[cur_type], alpha=0.5)
        ktop = k
        cur_type = prof.conv_mixing_type[k]


ax.legend(loc=3)
ax.set_ylim([0,2.5])
ax.set_xlim([6.8,7.7])
plt.show()

