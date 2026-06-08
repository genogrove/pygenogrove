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

    # Data-carrying grove: each interval carries an associated BED record
    print("\nBuilding a BedGrove (intervals with associated BED data)...")
    bed = pg.BedGrove(100)
    genes = [
        ("chr1", pg.Interval(1000, 1999), "BRCA1", 900, "+"),
        ("chr1", pg.Interval(3000, 3999), "TP53", 850, "-"),
    ]
    for chrom, iv, name, score, strand in genes:
        entry = pg.BedEntry(chrom, iv.start, iv.end + 1)  # BED end is half-open
        entry.name = name
        entry.score = score
        entry.strand = strand
        bed.insert(chrom, iv, entry)

    print("  Querying chr1 with Interval(1500, 1600)...")
    for hit in bed.intersect(pg.Interval(1500, 1600), "chr1"):
        d = hit.data
        print(f"    {hit.value} -> {d.name} (score={d.score}, strand={d.strand})")

    bed_path = "example_genes.gg"
    print(f"\n  Serializing BedGrove to {bed_path} and reloading...")
    bed.serialize(bed_path)
    bed_reloaded = pg.BedGrove.deserialize(bed_path)
    hit = list(bed_reloaded.intersect(pg.Interval(1500, 1600), "chr1"))[0]
    print(f"  reloaded payload survived round-trip: {hit.data.name} "
          f"(score={hit.data.score})")

    # Data-carrying grove with GFF/GTF records (column-9 attributes)
    print("\nBuilding a GffGrove (intervals with associated GFF/GTF data)...")
    gff = pg.GffGrove(100)
    exon = pg.GffEntry("chr1", 1000, 2000, "exon")  # GFF coords are 1-based inclusive
    exon.source = "ensembl"
    exon.strand = "+"
    exon.attributes = {"gene_id": "ENSG1", "transcript_id": "ENST1", "exon_number": "1"}
    gff.insert("chr1", exon)  # 2-arg insert derives Interval(999, 1999) from the GFF coords

    print("  Querying chr1 with Interval(1500, 1600)...")
    for hit in gff.intersect(pg.Interval(1500, 1600), "chr1"):
        d = hit.data
        print(f"    {d.type} gene_id={d.get_gene_id()} "
              f"exon#={d.get_exon_number()} attrs={dict(d.attributes)}")

    # Reading a BED file with BedReader and loading it into a grove
    print("\nReading a BED file into a BedGrove with BedReader...")
    bed_file = "example_peaks.bed"
    with open(bed_file, "w") as fh:
        fh.write("chr1\t1000\t2000\tpeakA\t100\t+\n")
        fh.write("chr1\t3000\t4000\tpeakB\t200\t-\n")
        fh.write("chr2\t500\t900\tpeakC\t50\t.\n")

    peaks = pg.BedGrove(256)
    for e in pg.BedReader(bed_file):
        peaks.insert(e.chrom, e)              # 2-arg insert derives the key from BED coords
    print(f"  loaded {peaks.size()} peaks from {bed_file}")

    overlaps = list(peaks.intersect(pg.Interval(1500, 1500), "chr1"))
    print(f"  peak overlapping chr1:1500 -> {overlaps[0].data.name}")

    # Bulk insert: load many sorted records into a data-carrying grove at once.
    # Passing bare entries lets insert_bulk derive each Interval key from the
    # entry's own coordinates.
    print("\nBulk-loading sorted records into a BedGrove (insert_bulk)...")
    big = pg.BedGrove(256)
    entries = [pg.BedEntry("chr1", i * 100, i * 100 + 50) for i in range(1000)]
    keys = big.insert_bulk("chr1", entries, presorted=True)  # 10-20x faster than a loop
    print(f"  bulk-inserted {len(keys)} records; grove size = {big.size()}")

    # Nearest-neighbour (flanking) query: closest features on either side of a gap
    print("\nFinding the nearest non-overlapping neighbours (flanking)...")
    nn = peaks.flanking(pg.Interval(2100, 2200), "chr1")   # a gap with no peak
    if nn.predecessor is not None:
        pred = nn.predecessor
        gap = 2100 - pred.value.end - 1
        print(f"  nearest upstream: {pred.data.name} ({pred.value}), gap={gap}")
    if nn.successor is not None:
        succ = nn.successor
        gap = succ.value.start - 2200 - 1
        print(f"  nearest downstream: {succ.data.name} ({succ.value}), gap={gap}")

    print("\nExample completed successfully!")


if __name__ == "__main__":
    main()
