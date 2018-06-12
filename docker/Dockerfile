FROM centos:7

LABEL description="contains dependencies for Harlem installation and testing"

# Install dependencies
# Dependencies for building:
#      autoconf make gcc gcc-c++ gcc-gfortran swig python python-devel zlib-devel
# Dependencies for testing:
#      python-pip
RUN \
    yum -y update && \
    yum -y install --setopt=tsflags=nodocs epel-release && \
    yum -y install --setopt=tsflags=nodocs \
        vim mc \
        autoconf make gcc gcc-c++ gcc-gfortran swig python python-devel zlib-devel \
        python-pip && \
    pip install -U pytest pytest-datadir

# nsimakov/harlem_ready

CMD ["bash"]
