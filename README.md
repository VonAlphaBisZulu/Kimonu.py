# Kimonu.py
Identify concordant and **KI**netic **MO**dules through sparse exact **NU**llspace computations

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/VonAlphaBisZulu/Kimonu.py/main?labpath=examples%2Fkimonu_examples.ipynb)

## Try it

[`examples/kimonu_examples.ipynb`](examples/kimonu_examples.ipynb) runs every call on the networks of Figures 1
and 2 of the manuscript: concordance, kinetic modules, ACR and ACRR by both routes, their exact values at given
rate constants, the rate constants for given values, and the values as formulas in the rate constants. The
Binder badge opens it in the browser, with nothing to install.

```
pip install -r binder/requirements.txt
pip install .
```

Kimonu needs `RationalMatrix.from_fractions` from straindesign, which is on its GitHub main but not yet in a
release (1.19.1), so `binder/requirements.txt` installs straindesign from there.
