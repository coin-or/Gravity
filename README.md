[![License](https://img.shields.io/badge/License-BSD--3-brightgreen.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Build Status](https://travis-ci.com/coin-or/Gravity.svg?branch=master)](https://travis-ci.com/coin-or/Gravity)
![CMake](https://github.com/coin-or/Gravity/workflows/CMake/badge.svg)
[![Code Coverage](https://codecov.io/gh/coin-or/gravity/branch/master/graph/badge.svg)](https://codecov.io/gh/coin-or/Gravity)
[![download](https://img.shields.io/badge/download%20%20-latest-blue.svg)](https://github.com/coin-or/Gravity/releases)

<a href="https://goo.gl/f7QLcS"><img alt="Chat on Slack" height="40" width="172" src="https://platform.slack-edge.com/img/sign_in_with_slack.png" srcset="https://platform.slack-edge.com/img/sign_in_with_slack.png 1x, https://platform.slack-edge.com/img/sign_in_with_slack@2x.png 2x" /></a>

<p align="center">
<img src="https://static.wixstatic.com/media/c6cff5_dd7659693c6247dc8eb8605d3dca95e8~mv2_d_3300_2550_s_4_2.png/v1/crop/x_1058,y_575,w_1183,h_1225/fill/w_288,h_298,al_c,usm_0.66_1.00_0.01/c6cff5_dd7659693c6247dc8eb8605d3dca95e8~mv2_d_3300_2550_s_4_2.png" width="250">
</p>
<H2 align="center"> Mathematical Modeling for Optimization and Machine Learning </H2>

<p align="center"> Created by Hassan Hijazi. </p>

<H2 align="center"> www.gravityopt.com </H2>



## License

Gravity is licensed under the BSD 3-Clause License. Please see the [LICENSE](https://github.com/coin-or/Gravity/blob/master/LICENSE) file for details.

[<img 
src="https://static.wixstatic.com/media/c6cff5_083fff4f0fa94b4b98b6790b18e7af8b~mv2.png/v1/fill/w_210,h_137,al_c,usm_0.66_1.00_0.01/c6cff5_083fff4f0fa94b4b98b6790b18e7af8b~mv2.png" width="100">](https://paypal.me/hlhijazi)

## Citing
The original paper was presentend at the Machine Learning Open Source Software Workshop at NeurIPS 2018, a longer version of the paper can be downloaded [here](https://791a4f37-01ef-43ce-b940-f17c763418b1.filesusr.com/ugd/c6cff5_e4889c3e27b54023a70a8c0496ff90a0.pdf).

Bibtex ref:
`@article{Gravity,
  title={Gravity: A Mathematical Modeling Language for Optimization and Machine Learning},
  author={Hassan Hijazi and Guanglei Wang and Carleton Coffrin},
  journal={Machine Learning Open Source Software Workshop at NeurIPS 2018},
  year={2018},
  note = {Available at \url{www.gravityopt.com}.},
  publisher={The Thirty-second Annual Conference on Neural Information Processing Systems (NeurIPS)}
}`

## Contributors
Hassan Hijazi, Los Alamos National Laboratory, The Australian National University | hlh@lanl.gov

Guanglei Wang, Ksenia Bestuzheva, Carleton Coffrin, Smitha Gopinath, Mertcan Yetkin

*****************************
See [INSTALL.md](INSTALL.md) for instructions on compiling Gravity

After running make, the Gravity executables can be found under Gravity/bin/
*****************************

Getting Started
-----------
First, you will need to install an IDE, I recommend to choose among the following:

[<img src="media/visual_studio.jpg" width="70">](https://www.visualstudio.com/downloads/) ||
[<img src="media/clion.jpg" width="50">](https://www.jetbrains.com/clion/) || 
[<img src="media/Xcode.png" width="50">](https://developer.apple.com/xcode/downloads/) || 
[<img src="media/eclipse-800x188.png" width="120">](https://www.eclipse.org/downloads/packages/release/2018-09/r/eclipse-ide-cc-developers)

Then, follow the instructions presented in [INSTALL.md](INSTALL.md).

The model below was implemented in Xcode:

![cover-example](media/Kapture_Stable_Set.gif)


Some Numerical Results:
-----------
Performance Profile on ACOPF
-----------

The first figure below is a performance profile illustrating percentage of instances solved as a function of time.
The figure compares Gravity, [JuMP](http://www.juliaopt.org/JuMP.jl/latest/index.html) and AMPL's NL interface (used by [AMPL](http://ampl.com/) and [Pyomo](http://www.pyomo.org/)) on all standard instances found in the [PGLIB](https://github.com/power-grid-lib/pglib-opf) benchmark library.

![Performance Profile on ACOPF](https://static.wixstatic.com/media/c6cff5_9b2b29e8a33840c59902fc95ffabf3ed~mv2.png/v1/crop/x_0,y_0,w_1064,h_600/fill/w_869,h_490,al_c,usm_0.66_1.00_0.01/c6cff5_9b2b29e8a33840c59902fc95ffabf3ed~mv2.png)

The figure below compares model build time between Gravity and [JuMP](http://www.juliaopt.org/JuMP.jl/latest/index.html) on the [PGLIB](https://github.com/power-grid-lib/pglib-opf) benchmarks.

![Model Build Time on ACOPF](https://static.wixstatic.com/media/c6cff5_27ee822625f24072b01110748c6f3923~mv2.jpg)

-----------
Performance Profile on Inverse Ising Model
-----------


![Performance Profile on Inverse Ising](https://static.wixstatic.com/media/c6cff5_e38e7a012b104dc0ba19fec1e32c10ad~mv2.png/v1/crop/x_0,y_0,w_1058,h_600/fill/w_863,h_489,al_c,usm_0.66_1.00_0.01/c6cff5_e38e7a012b104dc0ba19fec1e32c10ad~mv2.png)


Click [here](www.gravityopt.com) for more details.
