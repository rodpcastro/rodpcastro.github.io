```python
from expint import expi
import mpmath
import matplotlib.pyplot as plt
import numpy as np

mpmath.dps = 16 # mpmath precision

xv = np.logspace(-10, 2, 1000, dtype=np.float64)
ei_expint = np.empty(xv.shape, dtype=np.float64)
ei_mpmath = np.empty(xv.shape, dtype=np.float64)

for i, x in enumerate(xv):
    ei_expint[i] = expi(x)
    ei_mpmath[i] = mpmath.ei(x)

ei_diff = np.abs(ei_expint - ei_mpmath)/np.abs(ei_mpmath)

fig, ax = plt.subplots(figsize=(7,4))
ax.loglog(xv, ei_diff, 'k')
ax.set_xticks(np.logspace(-10, 2, 5))
ax.set_yticks(np.logspace(-16, -13, 4))
ax.set_xlabel('x')
ax.set_ylabel(r'$ \frac{|Ei_{expint} - Ei_{mpmath}|}{|Ei_{mpmath}|} $')
ax.set_title('Relative error')
plt.savefig('expi_error.svg')
```
