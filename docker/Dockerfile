FROM ubuntu:18.04

LABEL description="contains dependencies for Harlem installation and testing"

WORKDIR /root

# Install dependencies
# Dependencies for building:
#      autoconf make gcc gcc-c++ gcc-gfortran swig python python-devel zlib-devel
# Dependencies for testing:
#      python-pip
RUN \
    apt-get update && \
    apt-get install -y g++ python3 python3-dev python3-pip\
        make autoconf swig zlib1g-dev\
        vim mc && \
    pip3 install -U pytest pytest-datadir

# nsimakov/harlem_ready

CMD ["bash"]
