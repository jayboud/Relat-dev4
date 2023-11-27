import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


M = 1


def equations(p,a):
    r,h,k = p
    V_eff = -M/r + (h**2 - a**2*(k**2 - 1))/(2*r**2) - M*(h-a*k)**2/r**3
    return [V_eff - (k**2-1)/2,
            M - (h**2 - a**2*(k**2 - 1))/r + 3*M*(h-a*k)**2/r**2,
            h**2 - a**2*(k**2 -1) - 6*M*(h-a*k)**2/r]


# valeur h+
r_sp,h_sp,k_sp = [6],[3],[-0.05]  # guesses
a = np.linspace(0,1,300)
for a_i in a:
    rs,hs,ks = fsolve(equations, (r_sp[-1],h_sp[-1],k_sp[-1]),args=a_i)  # guess avec derniere valeur trouvee
    r_sp.append(rs)
    h_sp.append(hs)
    k_sp.append(ks)
r_sp.pop(0)  # enleve le guess initial
h_sp.pop(0)
k_sp.pop(0)
# valeur h-
r_sm,h_sm,k_sm = [6],[-3],[-0.05]  # guesses
a = np.linspace(0,1,300)
for a_i in a:
    rs,hs,ks = fsolve(equations, (r_sm[-1],h_sm[-1],k_sm[-1]),args=a_i)  # guess avec derniere valeur trouvee
    r_sm.append(rs)
    h_sm.append(hs)
    k_sm.append(ks)
r_sm.pop(0)  # enleve le guess initial
h_sm.pop(0)
k_sm.pop(0)
# graphique
fig,axs = plt.subplots(3)
for ax in axs:
    ax.tick_params(length=4, width=1)
axs[0].plot(a, np.array(r_sp),label="co-rotation")
axs[0].plot(a, np.array(r_sm),label="contre-rotation")
axs[0].set_ylabel(r"$r$")
axs[0].legend()
axs[0].tick_params(
    axis='x',
    which='both',
    bottom=False,
    top=False,
    labelbottom=False)
axs[1].plot(a, np.array(h_sp))
axs[1].plot(a, np.array(h_sm))
axs[1].set_ylabel(r"$h$")
axs[1].tick_params(
    axis='x',
    which='both',
    bottom=False,
    top=False,
    labelbottom=False)
axs[2].plot(a, np.array(k_sp)-1)
axs[2].plot(a, np.array(k_sm)-1)
axs[2].set_ylabel(r"$k-1$")
axs[2].set_xlabel(r"$a$")
plt.savefig("num2")
