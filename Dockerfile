FROM sagemath/sagemath:8.0-2

ENV NB_USER=sage
ENV HOME /home/sage

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
USER root
RUN chown -R ${NB_USER}:${NB_USER} ${HOME}
USER ${NB_USER}

# Install ffisom and its dependencies
RUN git submodule init && \
  git submodule update && \
  cd implementation && \
  sage -sh -c make # Build the library

# Make it available system-wide
ENV LD_LIBRARY_PATH "${HOME}/implementation:$LD_LIBRARY_PATH"
ENV PYTHONPATH      "${HOME}/implementation:$PYTHONPATH"

# Install dependencies for plots
RUN sage -pip install pandas==0.20

EXPOSE 8888
CMD ["jupyter", "notebook", "--notebook-dir=notebooks", "--ip", "'*'", "--port", "8888"]
