#!/usr/bin/env bash
#
# Build htslib from source inside the cibuildwheel manylinux container.
#
# genogrove (pulled in via add_subdirectory) does `pkg_check_modules(HTSLIB
# REQUIRED htslib)`, so htslib must be present at configure/build time even
# though pygenogrove only binds the interval/grove surface. The manylinux_2_28
# image (AlmaLinux 8) has no htslib package, so we build it from source into
# /usr/local; auditwheel then bundles libhts into the repaired wheel.
#
# Runs once per container via cibuildwheel's `before-all`.
set -euo pipefail

HTSLIB_VERSION="${HTSLIB_VERSION:-1.21}"

# Build deps for htslib (+ the compression/network backends it links).
dnf install -y \
    autoconf automake make gcc perl \
    zlib-devel bzip2-devel xz-devel libcurl-devel openssl-devel

curl -fsSL -o /tmp/htslib.tar.bz2 \
    "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"
mkdir -p /tmp/htslib
tar -xjf /tmp/htslib.tar.bz2 -C /tmp/htslib --strip-components=1

cd /tmp/htslib
./configure
make -j"$(nproc)"
make install            # -> /usr/local/lib/libhts.so, /usr/local/lib/pkgconfig/htslib.pc
ldconfig