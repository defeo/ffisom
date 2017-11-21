FROM sagemath/sagemath:8.0

# Inspired from https://mybinder.readthedocs.io/en/latest/dockerfile.html#preparing-your-dockerfile

ENV NB_USER=sage
ENV HOME /home/sage

# Upgrade to jupyter 5.* as required by mybinder
RUN sage -pip install "notebook>=5" "ipykernel>=4.6"
# Install ffisom and its dependencies
# Clone the repo
RUN git clone https://github.com/defeo/ffisom.git
RUN cd ffisom
# Make sure ellmul sources are pulled
RUN git submodule init
RUN git submodule update
# Build the library
RUN cd implementation
RUN sage -sh -c make
# Make it available system-wide
ENV LD_LIBRARY_PATH "${HOME}/ffisom/implementation:$LD_LIBRARY_PATH"
ENV PYTHONPATH      "${HOME}/ffisom/implementation:$PYTHONPATH"
RUN cd ../..

USER root
# This will eventually be lifted upstream to sagemath/sagemath
COPY jupyter jupyter-notebook /usr/bin/

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
RUN chown -R ${NB_USER}:${NB_USER} ${HOME}
USER ${NB_USER}
