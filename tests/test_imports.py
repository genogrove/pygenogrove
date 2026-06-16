"""
Build sanity check: the compiled extension imports and exposes its public types.

This has no genogrove C++ counterpart — it guards the Python packaging/build
itself. It is a hard failure (not a skip) on purpose: if the compiled extension
is not installable/importable, the rest of the suite would otherwise pass by
silently skipping every test, hiding a broken build.
"""


def test_imports_core():
    """The standard genomic-coordinate Grove (JSON payload) and its types."""
    import pygenogrove
    assert hasattr(pygenogrove, 'GenomicCoordinate')
    assert hasattr(pygenogrove, 'Grove')
    assert hasattr(pygenogrove, 'Key')
    assert hasattr(pygenogrove, 'QueryResult')
    assert hasattr(pygenogrove, 'FlankingResult')


def test_interval_is_removed():
    """The old Interval-keyed surface is gone (genomic_coordinate is standard)."""
    import pygenogrove
    assert not hasattr(pygenogrove, 'Interval')
    assert not hasattr(pygenogrove, 'GenomicCoordinateGrove')


def test_imports_bed_data_grove():
    """The typed BED data grove (genomic_coordinate keyed) and its value types."""
    import pygenogrove
    assert hasattr(pygenogrove, 'BedGrove')
    assert hasattr(pygenogrove, 'BedKey')
    assert hasattr(pygenogrove, 'BedQueryResult')
    assert hasattr(pygenogrove, 'BedFlankingResult')
    assert hasattr(pygenogrove, 'BedEntry')
    assert hasattr(pygenogrove, 'BlockInfo')
    assert hasattr(pygenogrove, 'ThickInfo')
    assert hasattr(pygenogrove, 'RgbColor')


def test_imports_gff_data_grove():
    """The typed GFF data grove (genomic_coordinate keyed) and its value types."""
    import pygenogrove
    assert hasattr(pygenogrove, 'GffGrove')
    assert hasattr(pygenogrove, 'GffKey')
    assert hasattr(pygenogrove, 'GffQueryResult')
    assert hasattr(pygenogrove, 'GffFlankingResult')
    assert hasattr(pygenogrove, 'GffEntry')
    assert hasattr(pygenogrove, 'GffFormat')


def test_imports_file_readers():
    """The BED/GFF/BAM/FASTA file readers are exposed."""
    import pygenogrove
    assert hasattr(pygenogrove, 'BedReader')
    assert hasattr(pygenogrove, 'GffReader')
    assert hasattr(pygenogrove, 'BamReader')
    assert hasattr(pygenogrove, 'FastaReader')
    assert hasattr(pygenogrove, 'FastaEntry')
    assert hasattr(pygenogrove, 'FastaIndex')


def test_imports_sam_types():
    """The SAM/BAM alignment value types are exposed."""
    import pygenogrove
    assert hasattr(pygenogrove, 'SamEntry')
    assert hasattr(pygenogrove, 'AlignmentFlags')
    assert hasattr(pygenogrove, 'SamFlags')


def test_imports_point_key_types():
    """The point key types (numeric, kmer) and their groves are exposed."""
    import pygenogrove
    for name in ('Numeric', 'NumericGrove', 'NumericKey', 'NumericQueryResult',
                 'NumericFlankingResult',
                 'Kmer', 'KmerGrove', 'KmerKey', 'KmerQueryResult',
                 'KmerFlankingResult'):
        assert hasattr(pygenogrove, name), name


def test_imports_registry():
    """The universal interning registry singleton is exposed."""
    import pygenogrove
    assert hasattr(pygenogrove, 'Registry')


def test_imports_filetype_detector():
    """The file-type detector and its enums are exposed."""
    import pygenogrove
    assert hasattr(pygenogrove, 'FiletypeDetector')
    assert hasattr(pygenogrove, 'Filetype')
    assert hasattr(pygenogrove, 'CompressionType')