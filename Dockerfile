FROM sagemath/sagemath:8.0

# Inspired from https://mybinder.readthedocs.io/en/latest/dockerfile.html#preparing-your-dockerfile

ENV NB_USER=sage
ENV HOME /home/sage

# Upgrade to jupyter 5.* as required by mybinder
RUN sage -pip install "notebook>=5" "ipykernel>=4.6"
# Install ffisom and its dependencies
RUN git clone https://github.com/defeo/ffisom.git && \
  cd ffisom && \
  git submodule init && \
  git submodule update && \
  cd implementation && \
  sage -sh -c make # Build the library

# Make it available system-wide
ENV LD_LIBRARY_PATH "${HOME}/ffisom/implementation:$LD_LIBRARY_PATH"
ENV PYTHONPATH      "${HOME}/ffisom/implementation:$PYTHONPATH"

USER root
# This will eventually be lifted upstream to sagemath/sagemath
#COPY jupyter jupyter-notebook /usr/bin/

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
RUN chown -R ${NB_USER}:${NB_USER} ${HOME}
USER ${NB_USER}
