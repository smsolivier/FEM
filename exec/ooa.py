#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

h, e = np.loadtxt('err', unpack=True)

fit = np.polyfit(np.log(h), np.log(e), 1)
print(fit[0])

f = lambda x: np.exp(fit[1]) * x**(fit[0])

plt.loglog(h, e, '-o')
plt.loglog(h, f(h))
plt.xlabel(r'$\sqrt{\mathrm{Average \ Element \ Area}}$')
plt.ylabel('Max Absolute Error')
plt.show()