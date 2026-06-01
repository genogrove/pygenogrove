#!/usr/bin/env python3
"""
Example usage of pygenogrove Python bindings.

This script demonstrates basic operations with the Grove and Interval classes.
"""

import sys

try:
    import pygenogrove as pg
except ImportError:
    print("Error: pygenogrove module not found.")
    print("Please build the module first:")
    print("  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release")
    print("  cmake --build build")
    print("\nThen add the build directory to your PYTHONPATH:")
    print("  export PYTHONPATH=$(pwd)/build:$PYTHONPATH")
    sys.exit(1)


def main():
    """Run example operations."""
    print(f"pygenogrove version: {pg.__version__}")
    print()

    # Create a grove with order 100
    print("Creating grove with order 100...")
    grove = pg.Grove(100)
    print(f"  {grove}")
    print()

    # Create some intervals
    print("Creating intervals...")
    intervals_chr1 = [
        pg.Interval(100, 200),
        pg.Interval(150, 250),
        pg.Interval(300, 400),
        pg.Interval(350, 450),
    ]

    intervals_chr2 = [
        pg.Interval(1000, 1500),
        pg.Interval(2000, 2500),
    ]

    # Insert intervals
    print("Inserting intervals into chr1...")
    for interval in intervals_chr1:
        key = grove.insert("chr1", interval)
        print(f"  Inserted: {interval} -> Key: {key.value}")

    print("\nInserting intervals into chr2...")
    for interval in intervals_chr2:
        key = grove.insert("chr2", interval)
        print(f"  Inserted: {interval} -> Key: {key.value}")

    print(f"\nTotal intervals in grove: {grove.size()}")
    print()

    # Query for overlapping intervals
    print("Querying for overlaps...")

    # Query chr1
    query1 = pg.Interval(175, 325)
    print(f"\nQuery chr1 with {query1}:")
    results1 = grove.intersect(query1, "chr1")
    print(f"  Found {len(results1)} overlapping intervals:")
    for key in results1:
        print(f"    {key.value}")

    # Query chr2
    query2 = pg.Interval(1200, 1800)
    print(f"\nQuery chr2 with {query2}:")
    results2 = grove.intersect(query2, "chr2")
    print(f"  Found {len(results2)} overlapping intervals:")
    for key in results2:
        print(f"    {key.value}")

    # Query all chromosomes
    query3 = pg.Interval(200, 250)
    print(f"\nQuery all chromosomes with {query3}:")
    results3 = grove.intersect(query3)
    print(f"  Found {len(results3)} overlapping intervals across all chromosomes:")
    for key in results3:
        print(f"    {key.value}")

    # Test overlap detection
    print("\nTesting interval overlap detection...")
    interval_a = pg.Interval(100, 200)
    interval_b = pg.Interval(150, 250)
    interval_c = pg.Interval(300, 400)

    print(f"  {interval_a} overlaps {interval_b}? {pg.Interval.overlaps(interval_a, interval_b)}")
    print(f"  {interval_a} overlaps {interval_c}? {pg.Interval.overlaps(interval_a, interval_c)}")

    # Graph overlay: connect intervals with directed edges
    print("\nBuilding a graph overlay (e.g. exon -> exon -> enhancer)...")
    gene = pg.Grove(3)
    exon1 = gene.insert("chr1", pg.Interval(1000, 1200))
    exon2 = gene.insert("chr1", pg.Interval(1400, 1600))
    enhancer = gene.add_external_key(pg.Interval(5000, 5500))  # not in the index

    gene.add_edge(exon1, exon2)
    gene.add_edge(exon2, enhancer)
    print(f"  edges: {gene.edge_count()}, out_degree(exon1): {gene.out_degree(exon1)}")
    print(f"  exon2 -> {[str(n.value) for n in gene.get_neighbors(exon2)]}")

    # Serialize the connected grove and read it back
    out_path = "example_gene.gg"
    print(f"\nSerializing grove to {out_path} and reloading...")
    gene.serialize(out_path)
    reloaded = pg.Grove.deserialize(out_path)
    src = list(reloaded.intersect(pg.Interval(1000, 1200), "chr1"))[0]
    print(f"  reloaded size: {reloaded.size()}, edges: {reloaded.edge_count()}")
    print(f"  reloaded exon1 -> {[str(n.value) for n in reloaded.get_neighbors(src)]}")

    print("\nExample completed successfully!")


if __name__ == "__main__":
    main()
