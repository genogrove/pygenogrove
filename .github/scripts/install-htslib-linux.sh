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
# SHA256 of the official htslib-${HTSLIB_VERSION}.tar.bz2 release asset. Pinned so
# the download is verified before extraction (supply-chain guard in the publish
# path). Keep in sync with HTSLIB_VERSION — re-verify against the samtools release
# when bumping. An override of HTSLIB_VERSION without HTSLIB_SHA256 fails closed.
HTSLIB_SHA256="${HTSLIB_SHA256:-84b510e735f4963641f26fd88c8abdee81ff4cb62168310ae716636aac0f1823}"

# Build deps for htslib (+ the compression/network backends it links).
dnf install -y \
    autoconf automake make gcc perl \
    zlib-devel bzip2-devel xz-devel libcurl-devel openssl-devel

curl -fsSL -o /tmp/htslib.tar.bz2 \
    "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"
# Verify the tarball against the pinned checksum; abort (no extraction) on mismatch.
echo "${HTSLIB_SHA256}  /tmp/htslib.tar.bz2" | sha256sum -c -
mkdir -p /tmp/htslib
tar -xjf /tmp/htslib.tar.bz2 -C /tmp/htslib --strip-components=1

cd /tmp/htslib
./configure
make -j"$(nproc)"
make install            # -> /usr/local/lib/libhts.so, /usr/local/lib/pkgconfig/htslib.pc
ldconfig