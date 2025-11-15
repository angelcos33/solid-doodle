FROM ubuntu:latest
LABEL Description="Build Environment"

ENV HOME=/root
ENV DISPLAY=host.docker.internal:0

SHELL [ "/bin/bash", "-c" ]

RUN apt-get update && apt-get -y --no-install-recommends install \
    build-essential \
    clang \
    cmake \
    cmake-curses-gui \
    gdb \
    wget \
    git-all \
    gnuplot \
    vim \
    cppcheck

## Cleanup cached apt data we don't need anymore
RUN apt-get autoremove -y && apt-get clean && \
rm -rf /var/lib/apt/lists/*