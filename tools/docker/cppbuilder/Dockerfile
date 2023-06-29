FROM ubuntu:22.04
ENV TZ=Europe/Lisbon
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update -y && apt upgrade -y \
&& apt install -y build-essential lld make git wget unzip cmake ninja-build gdb \
        ocaml ocamlbuild automake autoconf libtool wget python-is-python3 libssl-dev git cmake perl \
        python3 python3-pip lcov gcovr cppcheck \
&& apt install -y clang clang-13 libc++-13-dev libc++abi-13-dev clang-format clang-format-13 \
        libx86-dev libclang-common-13-dev libclang-common-14-dev \
&& apt install -y libgtest-dev libboost-dev libssl-dev \
&& apt install -y doxygen \
&& apt install -y python3-matplotlib python3-numpy \
&& apt install -y libssl-dev libcurl4-openssl-dev protobuf-compiler libprotobuf-dev \
        debhelper cmake reprepro unzip pkgconf libboost-dev libboost-system-dev libboost-thread-dev \
        protobuf-c-compiler libprotobuf-c-dev lsb-release libsystemd0 \
&& rm -rf /var/lib/apt/lists/*

RUN pip install cppcheck-codequality


RUN wget 'https://github.com/boyter/scc/releases/download/v3.0.0/scc-3.0.0-x86_64-unknown-linux.zip'
RUN unzip 'scc-3.0.0-x86_64-unknown-linux.zip' -d /
RUN chmod +x /scc

RUN mkdir /dockerfiles
RUN cd /dockerfiles
RUN git clone 'https://github.com/intel/linux-sgx' /dockerfiles/linux-sgx && cd /dockerfiles/linux-sgx && git checkout 7cd5179
RUN cd /dockerfiles/linux-sgx && make preparation
RUN cd /dockerfiles/linux-sgx/external/toolset/ubuntu20.04; cp ./* /usr/local/bin
RUN cd /dockerfiles/linux-sgx && make sdk && make sdk_install_pkg
RUN cd /dockerfiles/linux-sgx/linux/installer/bin/ && echo yes | ./sgx_linux_x64_sdk_2.18.*.bin
RUN echo "source /dockerfiles/linux-sgx/linux/installer/bin/sgxsdk/environment" > /startsgxenv.sh
RUN chmod +x /startsgxenv.sh

SHELL ["/bin/bash", "-c"] 
RUN source /startsgxenv.sh && cd /dockerfiles/linux-sgx && make psw && make psw_install_pkg
RUN source /startsgxenv.sh && cd /dockerfiles/linux-sgx && make deb_psw_pkg
RUN cd /dockerfiles/linux-sgx/linux/installer/bin/ && echo yes | ./sgx_linux_x64_psw_2.18.*.bin || echo ignoring error
RUN ln -s /usr/lib/libsgx_urts.so /usr/lib/libsgx_urts.so.2

WORKDIR /builder/
CMD ["bash"]
