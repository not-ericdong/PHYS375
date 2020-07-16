import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def solving(y, x, n):
    return[y[1]*x**2, -y[1]*2*x - y[0]**n]

c0 = [0,1]
xs = np.linspace(0,10,200)
n=0
soln0 = odeint(solving, c0, xs, args=(n,))
n=1
soln1 = odeint(solving, c0, xs, args=(n,))
n=2
soln2 = odeint(solving, c0, xs, args=(n,))
n=3
soln3 = odeint(solving, c0, xs, args=(n,))
n=4>>>
soln4 = odeint(solving, c0, xs, args=(n,))
n=5
soln5 = odeint(solving, c0, xs, args=(n,))
ys0 = soln0[:,0]
ys1 = soln1[:,0]
ys2 = soln2[:,0]
ys3 = soln3[:,0]
ys4 = soln4[:,0]
ys5 = soln5[:,0]

plt.plot(xs,ys0)
plt.plot(xs,ys1)
plt.plot(xs,ys2)
plt.plot(xs,ys3)
plt.plot(xs,ys4)
plt.plot(xs,ys5)

plt.show()
