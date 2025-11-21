

PIC-BLAS is an extension of the PIC library aimed at creating a nice and simple BLAS interface that uses
the native Fortran modules instead of going through C.

The main reason is to explore the limits of the Fortran modules and actively report these issues to the
right developers and channels to get them fixed.


## Building (TLDR) on Linux

I assume you have experience building things here. If you don't please go further down for a more verbose explanation.

Briefly, for a minimal build  you need:

- CMake (at least 3.22) or the Fortran Package Manager (at least 0.12.0)
- A Fortran coblasler
- An internet connection to pull the dependencies
- An BLAS library installed


### Building with CMake

The top level directory is `$PICBLAS_BASE` which is the `pic-blas/` directory that was cloned or unpacked. I assume you are here. And
`$PICBLAS_ROOT` is the path to where you'd like pic to be installed, for example `export PICBLAS_ROOT=$HOME/install/pic-blas/dev/`

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$PICBLAS_ROOT ../
make -j install
```

To run the tests, from the build directory simply run: `ctest`

### Building with FPM

The easiest way to build with the FPM is to use the coblasler wrappers provided for Fortran i.e. `blasfort` or `blasf90`.

These will make life simpler with linking and coblaslation flags.

The same `$PICBLAS_BASE` and `$PICBLAS_ROOT` will be used here. Simply:

```
fpm install --prefix $PICBLAS_ROOT --profile release --coblasler blasfort
```

To run the tests: `fpm test --profile release --coblasler blasfort`


## Building and dependencies

There's two build systems included in the present version, CMake and the [Fortran Package Manager](https://fpm.fortran-lang.org/index.html).

The dependencies of the project are, as of now, CMake (if using cmake), BLAS, and OpenMP.

## Documentation

The code itself is documented using [FORD](https://forddocs.readthedocs.io/en/stable/) and the documentation is available [here](https://jorgeg94.github.io/pic/).

Comments in the code that are prefixed with `!!` are considered documentation comments and will be processed by FORD. Comments without that prefix are considered regular comments and will not be processed by FORD. So, please do not use `!!` for your comments unless you want them to be included in the documentation.


## Contributing

Please see the PIC [contributing guidelines](https://jorgeg94.github.io/pic/page/contributing.html) for information on how to contribute to the project.

## Using PIC-BLAS in your work

### Fortran Package Manager

Simply add:
```

[dependencies]
pic-blas = { git = "https://github.com/JorgeG94/pic-blas.git", branch = "main"}
```

to your fpm.toml file and you'll be able to checkout and use pic-blas.

### CMake

For CMake it is a bit more complex, since you'll need to pull the dependency. You can see this [template repo](https://github.com/JorgeG94/pic-app-sample), which
serves as an example on pulling and using the code inside your build system.

### Static/Shared linking via CMake

pic is coblasled with "CMake symbols", i.e. it will be findable by a CMake package provided you do the right things. Using
`find_package(pic REQUIRED)` will use CMake intrinsics to try to find the necessary things to link to pic. pic comes with the
target `pic-blas::pic-blas` that you can use for your `target_link_libraries(${my_target} PRIVATE pic-blas::pic-blas)`. This will import
all includes, modules, libs, etc.

How CMake finds `pic-blas::pic-blas` depends on the policy `CMP0074`, this controls if the variables `pic-blas_DIR` or `pic-blas_ROOT` are used
to search for the package. If you have set the policy: `cmake_policy(SET CMP0074 NEW)`, then `pic_ROOT` will also be used,
otherwise it is *IGNORED*. By setting `export pic-blas_ROOT=/path/to/where/pic-blas/is/installed` it will let CMake find the
necessary files it needs to just link pic-blas. Be careful that, as of now, pic-blas needs to be in small letter. All caps will fail to
find.
