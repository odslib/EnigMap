FROM debian:bookworm
ENV TZ=Europe/Lisbon
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update -y && apt upgrade -y
RUN apt install -y build-essential lld make git wget unzip cmake ninja-build
RUN apt install -y clang clang-13 libc++-13-dev libc++abi-13-dev 
RUN apt install -y libgtest-dev libboost-dev libssl-dev
RUN apt install -y lcov gcovr cppcheck
RUN apt install -y python3 python3-pip
RUN apt install -y doxygen

RUN wget 'https://github.com/boyter/scc/releases/download/v3.0.0/scc-3.0.0-x86_64-unknown-linux.zip'
RUN unzip 'scc-3.0.0-x86_64-unknown-linux.zip' -d /
RUN chmod +x /scc

RUN pip install cppcheck-codequality
RUN pip install matplotlib numpy


RUN apt install -y ocaml ocamlbuild automake autoconf libtool wget python-is-python3 libssl-dev git cmake perl
RUN mkdir /dockerfiles
RUN cd /dockerfiles
RUN git clone 'https://github.com/intel/linux-sgx' /dockerfiles/linux-sgx
RUN cd /dockerfiles/linux-sgx && make preparation
RUN cd /dockerfiles/linux-sgx/external/toolset/ubuntu20.04; cp ./* /usr/local/bin
RUN cd /dockerfiles/linux-sgx && make -j sdk && make -j sdk_install_pkg
RUN cd /dockerfiles/linux-sgx/linux/installer/bin/ && echo yes | ./sgx_linux_x64_sdk_2.17.*.bin
RUN echo "source /dockerfiles/linux-sgx/linux/installer/bin/sgxsdk/environment" > /startsgxenv.sh
RUN chmod +x /startsgxenv.sh

RUN apt install -y gdb
RUN apt install -y clang-format clang-format-13

WORKDIR /builder/
CMD ["bash"]
