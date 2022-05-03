[![License](https://img.shields.io/badge/License-BSD--3-brightgreen.svg)](https://opensource.org/licenses/BSD-3-Clause)
<p align="left">
<img src="https://static.wixstatic.com/media/c6cff5_dd7659693c6247dc8eb8605d3dca95e8~mv2_d_3300_2550_s_4_2.png/v1/crop/x_1058,y_575,w_1183,h_1225/fill/w_288,h_298,al_c,usm_0.66_1.00_0.01/c6cff5_dd7659693c6247dc8eb8605d3dca95e8~mv2_d_3300_2550_s_4_2.png" width="150">
</p>
<!-- <H2 align="center"> Mathematical Modeling for Optimization and Machine Learning</H2> -->



Getting Started
-----------

To compile Gravity see [INSTALL.md](INSTALL.md).

After building, the Gravity library can be found under `Gravity/lib`, and the executables (from [`Gravity/examples`](https://github.com/coin-or/Gravity/tree/master/examples)) can be found under `Gravity/bin/Release`

If you're new to coding, a few editor recommendations are below:

[Visual Studio](https://www.visualstudio.com/downloads/) | [Clion](https://www.jetbrains.com/clion/) | [Xcode](https://developer.apple.com/xcode/downloads/) | [Eclipse](https://www.eclipse.org/downloads/packages/release/2018-09/r/eclipse-ide-cc-developers)
:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:

Have a look at [`Gravity_test.cpp`](https://github.com/coin-or/Gravity/blob/master/examples/Gravity_test.cpp) and other models [`here`](https://github.com/coin-or/Gravity/tree/master/examples) to learn from example.

### Some Numerical Results:
#### Performance Profile on ACOPF

The first figure below is a performance profile illustrating percentage of instances solved as a function of time.
The figure compares Gravity, [JuMP](http://www.juliaopt.org/JuMP.jl/latest/index.html) and AMPL's NL interface (used by [AMPL](http://ampl.com/) and [Pyomo](http://www.pyomo.org/)) on all standard instances found in the [PGLIB](https://github.com/power-grid-lib/pglib-opf) benchmark library.

![Performance Profile on ACOPF](https://static.wixstatic.com/media/c6cff5_9b2b29e8a33840c59902fc95ffabf3ed~mv2.png/v1/crop/x_0,y_0,w_1064,h_600/fill/w_869,h_490,al_c,usm_0.66_1.00_0.01/c6cff5_9b2b29e8a33840c59902fc95ffabf3ed~mv2.png)

The figure below compares model build time between Gravity and [JuMP](http://www.juliaopt.org/JuMP.jl/latest/index.html) on the [PGLIB](https://github.com/power-grid-lib/pglib-opf) benchmarks.

![Model Build Time on ACOPF](https://static.wixstatic.com/media/c6cff5_27ee822625f24072b01110748c6f3923~mv2.jpg)

#### Performance Profile on Inverse Ising Model


![Performance Profile on Inverse Ising](https://static.wixstatic.com/media/c6cff5_e38e7a012b104dc0ba19fec1e32c10ad~mv2.png/v1/crop/x_0,y_0,w_1058,h_600/fill/w_863,h_489,al_c,usm_0.66_1.00_0.01/c6cff5_e38e7a012b104dc0ba19fec1e32c10ad~mv2.png)

## Citing
The original paper was presentend at the Machine Learning Open Source Software Workshop at NeurIPS 2018, a longer version of the paper can be downloaded [here](https://791a4f37-01ef-43ce-b940-f17c763418b1.filesusr.com/ugd/c6cff5_e4889c3e27b54023a70a8c0496ff90a0.pdf).

Bibtex ref:
```
@article{Gravity,
  title={Gravity: A Mathematical Modeling Language for Optimization and Machine Learning},
  author={Hassan Hijazi and Guanglei Wang and Carleton Coffrin},
  journal={Machine Learning Open Source Software Workshop at NeurIPS 2018},
  year={2018},
  note = {Available at \url{www.gravityopt.com}.},
  publisher={The Thirty-second Annual Conference on Neural Information Processing Systems (NeurIPS)}
}
```

## License

Gravity is licensed under the BSD 3-Clause License. Please see the [LICENSE](https://github.com/coin-or/Gravity/blob/master/LICENSE) file for details.

[<img 
src="https://static.wixstatic.com/media/c6cff5_083fff4f0fa94b4b98b6790b18e7af8b~mv2.png/v1/fill/w_210,h_137,al_c,usm_0.66_1.00_0.01/c6cff5_083fff4f0fa94b4b98b6790b18e7af8b~mv2.png" width="100">](https://paypal.me/hlhijazi)

## Contributors
See the list of contributors [here](https://github.com/coin-or/Gravity/graphs/contributors) 
