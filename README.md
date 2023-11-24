
for osqp:
needs CMake (sudo apt install CMake)
needs gcc

sudo apt-get upgrade libstdc++6

if osqp does not compile
LD_PRELOAD="/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/x86_64-linux-gnu/libcurl.so.4" matlab &
