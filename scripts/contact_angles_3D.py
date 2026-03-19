import numpy as np
import matplotlib.pyplot as plt

exact = np.array([160, 145, 130, 115, 100, 90, 80, 65, 50, 35, 20])
print(exact)

bb = np.array([
    [162.38523865697954 , 0.0008468718978636731],
    [145.80772439433525 , 0.0007154181427742473],
    [129.58576209463848 , 0.0005063159001423385],
    [114.23684051674628 , 0.000529950312874128],
    [99.36526721613015 , 0.0005232812114189957],
    [89.48806175860143 , 0.0005142026844008902],
    [79.5126336080325 , 0.0005038959592742817],
    [64.20691985220006 , 0.0005105433791813824],
    [48.14710966297012 , 0.0005111514509904475],
    [29.789579404605888 , 0.0007503264952192841],
])

nebb = np.array([
    [159.85846799349312 , 0.00025571215447963173],
    [144.6770199815064 , 0.0003614713203314116],
    [129.40656153199913 , 0.0005196429679911048],
    [114.49208477236327 , 0.0005242918154432095],
    [99.83183613215081 , 0.0005197799009767198],
    [90.0424361085652 , 0.0005135192962953246],
    [80.14285301575178 , 0.0005133140889038994],
    [65.04058241226939 , 0.0005091566885457212],
    [49.58716536238589 , 0.0005043912969318185],
    [33.137024688513776 , 0.0004956808300225869],
    [9.534904748404967 , 0.0004984940389841679]
])

plt.axhline(0, ls="--", color="black", zorder=0)
plt.plot(exact, np.abs(nebb[:,0]-exact), zorder=1)
plt.plot(exact[:-1], np.abs(bb[:,0]-exact[:-1]), zorder=1)
plt.scatter(exact, np.abs(nebb[:,0]-exact), label="NEBB", zorder=2)
plt.scatter(exact[:-1], np.abs(bb[:,0]-exact[:-1]), label="BB", zorder=2)
plt.legend()
plt.savefig("contact_angles_3D_tau=1.pdf")
plt.savefig("contact_angles_3D_tau=1.png")
plt.close()

plt.axhline(0, ls="--", color="black", zorder=0)
plt.plot(exact, nebb[:,1], zorder=1)
plt.plot(exact[:-1], bb[:,1], zorder=1)
plt.scatter(exact, nebb[:,1], label="NEBB", zorder=2)
plt.scatter(exact[:-1], bb[:,1], label="BB", zorder=2)
plt.legend()
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.savefig("spurious_currents_3D_tau=1.pdf")
plt.savefig("spurious_currents_3D_tau=1.png")