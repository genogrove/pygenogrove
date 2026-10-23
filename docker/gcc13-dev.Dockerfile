# Dev toolchain for pygenogrove: GCC 13 (genogrove's tested floor), htslib,
# and a Python with headers so pybind11 can configure.
#
# docker build -t pygenogrove-dev:gcc13 -f docker/gcc13-dev.Dockerfile docker
# CLion: Settings > Toolchains > Docker, image pygenogrove-dev:gcc13; then
# cmake -S . -B build -G Ninja && cmake --build build
# PYTHONPATH=build pytest tests/
FROM gcc:13.3

RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake ninja-build gdb pkg-config \
    python3 python3-dev python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --break-system-packages pytest

# htslib into /usr/local (zlib/bz2/lzma/curl dev libs come with the gcc image).
# Version and checksum match .github/scripts/install-htslib-linux.sh.
RUN wget -q https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 \
    && echo "84b510e735f4963641f26fd88c8abdee81ff4cb62168310ae716636aac0f1823  htslib-1.21.tar.bz2" | sha256sum -c - \
    && tar -xf htslib-1.21.tar.bz2 && cd htslib-1.21 \
    && ./configure && make -j"$(nproc)" && make install && ldconfig \
    && cd .. && rm -rf htslib-1.21 htslib-1.21.tar.bz2
