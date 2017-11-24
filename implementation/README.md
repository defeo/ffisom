# Software

This software contains implementations of the algorithms described in
the [paper](../paper).

Some algorithms are implemented in C, using the [Flint
library](http://flintlib.org/), others are implemented in
[SageMath](http://sagemath.org/).

We provide a unified SageMath interface to run all algorithms.

## Running the software

The easiest way to run this software is through the Jupyter interface
offered by
[Binder](https://mybinder.org/v2/gh/defeo/ffisom/master?filepath=notebooks%2Fexample.ipynb).
Just follow the link and start playing with it!

If you prefer running the software on your own hardware, keep reading.

### Requirements

SageMath, version 6.0 or greater, is required to run the software.

### Build and running instructions

Clone this repo

	git clone https://github.com/defeo/ffisom.git

CD to this directory

	cd ffisom/implementation
	
Update git submodules

	git submodule init
	git submodule update

Build the library within SageMath

	sage -sh -c make

Update the library paths (must be done each time a new terminal is open)

	export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
	export PYTHONPATH=$PWD:$PYTHONPATH

Lauch SageMath

	sage

or

	sage --notebook

Use `import ffisom` to import the library.  See examples in [the
Jupyter notebook](../notebooks/example.ipynb).

### Alternative installation using Docker

If you do not have SageMath installed, but you have Docker, you can
use the provided Docker container to install everything in one go.

Import the image

	docker pull defeo/ffisom

Run SageMath within the container

	docker run -it --rm jupyter notebook --notebook-dir=notebooks

After starting up, it will display a URL in the terminal. Copy it in
your browser and you're ready to go!
