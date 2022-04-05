|License| |Travis| |Github Actions| |Code Coverage| |download|

.. raw:: html

   <p align="center">

.. raw:: html

   </p>

.. raw:: html

   <H2 align="center">

Mathematical Modeling for Optimization and Machine Learning

.. raw:: html

   </H2>

.. raw:: html

   <p align="center">

Created by Hassan Hijazi - hlh@lanl.gov.

.. raw:: html

   </p>

.. raw:: html

   <H2 align="center">

www.gravityopt.com

.. raw:: html

   </H2>

License
-------

Gravity is licensed under the BSD 3-Clause License. Please see the
`LICENSE <https://github.com/coin-or/Gravity/blob/master/LICENSE>`__
file for details.

` <https://paypal.me/hlhijazi>`__

Citing
------

The original paper was presentend at the Machine Learning Open Source
Software Workshop at NeurIPS 2018, a longer version of the paper can be
downloaded
`here <https://791a4f37-01ef-43ce-b940-f17c763418b1.filesusr.com/ugd/c6cff5_e4889c3e27b54023a70a8c0496ff90a0.pdf>`__.

Bibtex ref:
``@article{Gravity,   title={Gravity: A Mathematical Modeling Language for Optimization and Machine Learning},   author={Hassan Hijazi and Guanglei Wang and Carleton Coffrin},   journal={Machine Learning Open Source Software Workshop at NeurIPS 2018},   year={2018},   note = {Available at \url{www.gravityopt.com}.},   publisher={The Thirty-second Annual Conference on Neural Information Processing Systems (NeurIPS)} }``

Contributors
------------

See the list of contributors
`here <https://github.com/coin-or/Gravity/graphs/contributors>`__

Getting Started
---------------

First, you will need to install an IDE, I recommend to choose among the
following:

` <https://www.visualstudio.com/downloads/>`__ \|\|
` <https://www.jetbrains.com/clion/>`__ \|\|
` <https://developer.apple.com/xcode/downloads/>`__ \|\|
` <https://www.eclipse.org/downloads/packages/release/2018-09/r/eclipse-ide-cc-developers>`__

Then, follow the instructions presented in `INSTALL.md <INSTALL.md>`__.

After building, the Gravity library can be found under ``Gravity/lib``,
and the executables (from
```Gravity/examples`` <https://github.com/coin-or/Gravity/tree/master/examples>`__)
can be found under ``Gravity/bin/Release``

The model below was implemented in Xcode:

.. figure:: media/Kapture_Stable_Set.gif
   :alt: cover-example

   cover-example

Some Numerical Results:
-----------------------

Performance Profile on ACOPF
----------------------------

The first figure below is a performance profile illustrating percentage
of instances solved as a function of time. The figure compares Gravity,
`JuMP <http://www.juliaopt.org/JuMP.jl/latest/index.html>`__ and AMPL’s
NL interface (used by `AMPL <http://ampl.com/>`__ and
`Pyomo <http://www.pyomo.org/>`__) on all standard instances found in
the `PGLIB <https://github.com/power-grid-lib/pglib-opf>`__ benchmark
library.

.. figure:: https://static.wixstatic.com/media/c6cff5_9b2b29e8a33840c59902fc95ffabf3ed~mv2.png/v1/crop/x_0,y_0,w_1064,h_600/fill/w_869,h_490,al_c,usm_0.66_1.00_0.01/c6cff5_9b2b29e8a33840c59902fc95ffabf3ed~mv2.png
   :alt: Performance Profile on ACOPF

   Performance Profile on ACOPF

The figure below compares model build time between Gravity and
`JuMP <http://www.juliaopt.org/JuMP.jl/latest/index.html>`__ on the
`PGLIB <https://github.com/power-grid-lib/pglib-opf>`__ benchmarks.

.. figure:: https://static.wixstatic.com/media/c6cff5_27ee822625f24072b01110748c6f3923~mv2.jpg
   :alt: Model Build Time on ACOPF

   Model Build Time on ACOPF

+--------------------------------------------+
| Performance Profile on Inverse Ising Model |
+--------------------------------------------+

.. figure:: https://static.wixstatic.com/media/c6cff5_e38e7a012b104dc0ba19fec1e32c10ad~mv2.png/v1/crop/x_0,y_0,w_1058,h_600/fill/w_863,h_489,al_c,usm_0.66_1.00_0.01/c6cff5_e38e7a012b104dc0ba19fec1e32c10ad~mv2.png
   :alt: Performance Profile on Inverse Ising

   Performance Profile on Inverse Ising

Click `here <www.gravityopt.com>`__ for more details.

.. |License| image:: https://img.shields.io/badge/License-BSD--3-brightgreen.svg
   :target: https://opensource.org/licenses/BSD-3-Clause
.. |Travis| image:: https://travis-ci.com/coin-or/Gravity.svg?branch=master
   :target: https://travis-ci.com/coin-or/Gravity
.. |Github Actions| image:: https://github.com/coin-or/Gravity/actions/workflows/cmake.yml/badge.svg
   :target: https://github.com/coin-or/Gravity/actions/workflows/cmake.yml
.. |Code Coverage| image:: https://codecov.io/gh/coin-or/gravity/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/coin-or/Gravity
.. |download| image:: https://img.shields.io/badge/download%20%20-latest-blue.svg
   :target: https://github.com/coin-or/Gravity/releases
