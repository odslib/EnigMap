uname -a
lsb_release -a
echo 'deb [arch=amd64] https://download.01.org/intel-sgx/sgx_repo/ubuntu jammy main' | sudo tee /etc/apt/sources.list.d/intel-sgx.list
wget -qO - https://download.01.org/intel-sgx/sgx_repo/ubuntu/intel-sgx-deb.key | sudo apt-key add -

sudo apt update -y && sudo apt upgrade -y
sudo apt install -y build-essential lld make git wget unzip cmake ninja-build gdb \
        ocaml ocamlbuild automake autoconf libtool wget python-is-python3 libssl-dev git cmake perl \
        python3 python3-pip lcov gcovr cppcheck \
clang clang-13 libc++-13-dev libc++abi-13-dev clang-format clang-format-13 \
libx86-dev libclang-common-13-dev libclang-common-14-dev libgtest-dev libboost-dev libssl-dev  doxygen  python3-matplotlib python3-numpy  libssl-dev libcurl4-openssl-dev protobuf-compiler libprotobuf-dev \
debhelper cmake reprepro unzip pkgconf libboost-dev libboost-system-dev libboost-thread-dev \
protobuf-c-compiler libprotobuf-c-dev lsb-release libsystemd0 

sudo apt-get install -y libsgx-urts libsgx-launch libsgx-enclave-common
sudo apt-get install -y libsgx-epid libsgx-quote-ex libsgx-dcap-ql
mkdir opt
cd opt
wget https://download.01.org/intel-sgx/latest/linux-latest/distro/ubuntu22.04-server/sgx_linux_x64_sdk_2.19.100.3.bin
chmod +x ./*
echo yes | ./sgx_linux_x64_sdk_2.19.100.3.bin
