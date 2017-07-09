# Multi-Scale (Sample) Entropy analysis (MSE)

The basis of this code was taken from the PhysioNet software base
(see [sampen](http://www.physionet.org/physiotools/sampen/) and [mse](http://www.physionet.org/physiotools/mse/)).

The concept of MSE using Sample Entropy (SampEn) is described in a [tutorial](https://physionet.org/physiotools/mse/tutorial/).

This package contains a C-library and a Python wrapper around it.

### C-library compilation

No external dependencies (except for the math library), just:
```bash
gcc -shared -o libsampen.so -O -Wall -fPIC sampen.c -lm
```
This command produces file "libsampen.so" in the current directory.
NOTE: do not change the name of the file, the Python wrapper searches for
this name.

### Python wrapper

In python, assuming the current directory with the repository files
and `libsampen.so` is **MultiScaleEntropy**:

```python
import MultiScaleEntropy
```

Now `MultiScaleEntropy` module contains the two most important functions:
* _sampen_: calculates SampEn estimate of a series
* _mse_: performs MSE analysis with SampEn as an entropy measure
