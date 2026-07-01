"""
Tests for the VCF/BCF reader (VcfReader) and its value types (VcfEntry,
SampleGenotype). Ports the behaviourally-observable cases from genogrove's
vcffile-test over the bound surface: record fields, SNP/indel predicates, INFO
and per-sample genotype parsing, header metadata, skip_filtered, and the
grove-loading helpers (to_coordinate / to_dict).
"""

import shutil
import subprocess

import pytest


def _pg():
    return pytest.importorskip("pygenogrove")


def _bgzip_tabix_vcf(tmp_path):
    """bgzip + tabix-index the shared _VCF fixture. Skips if the CLI tools are
    not on PATH. Returns the .vcf.gz path."""
    if not (shutil.which("bgzip") and shutil.which("tabix")):
        pytest.skip("bgzip/tabix not available")
    vcf = tmp_path / "calls.vcf"
    vcf.write_text(_VCF)
    subprocess.run(["bgzip", "-f", str(vcf)], check=True)
    gz = tmp_path / "calls.vcf.gz"
    subprocess.run(["tabix", "-p", "vcf", str(gz)], check=True)
    return str(gz)


_VCF = "\n".join([
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    '##FILTER=<ID=q10,Description="Quality below 10">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
    "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
               "FORMAT", "S1", "S2"]),
    "\t".join(["chr1", "100", "rs1", "A", "G", "50", "PASS", "DP=30",
               "GT:DP", "0/1:20", "1|1:25"]),
    "\t".join(["chr1", "200", ".", "AT", "A", ".", "q10", "DP=10",
               "GT:DP", "0/0:5", "./.:."]),
    "",
])


def _write_vcf(tmp_path):
    p = tmp_path / "calls.vcf"
    p.write_text(_VCF)
    return str(p)


def test_reads_records_and_core_fields(tmp_path):
    pg = _pg()
    records = list(pg.VcfReader(_write_vcf(tmp_path)))
    assert len(records) == 2

    a, b = records
    # POS is 1-based in the file; start is the 0-based htslib position.
    assert (a.chrom, a.start, a.ref, list(a.alt)) == ("chr1", 99, "A", ["G"])
    assert a.id == "rs1"
    assert a.end == 100  # start + len(REF)

    assert (b.start, b.ref, list(b.alt)) == (199, "AT", ["A"])
    assert b.end == 201
    assert b.id == ""  # "." -> empty


def test_qual_and_filter(tmp_path):
    pg = _pg()
    a, b = list(pg.VcfReader(_write_vcf(tmp_path)))
    assert a.qual_missing is False
    assert a.qual == pytest.approx(50.0)
    assert a.passed_filter() is True

    assert b.qual_missing is True  # QUAL "."
    assert b.passed_filter() is False
    assert list(b.filter) == ["q10"]


def test_snp_indel_predicates(tmp_path):
    pg = _pg()
    a, b = list(pg.VcfReader(_write_vcf(tmp_path)))
    assert a.is_snp() is True and a.is_indel() is False   # A -> G
    assert b.is_indel() is True and b.is_snp() is False    # AT -> A
    assert pg.VcfEntry.is_symbolic_allele("<DEL>") is True
    assert pg.VcfEntry.is_symbolic_allele("G") is False


def test_info_parsed(tmp_path):
    pg = _pg()
    a = next(iter(pg.VcfReader(_write_vcf(tmp_path))))
    # Integer INFO field (Number=1) -> list[int].
    assert a.info["DP"] == [30]


def test_samples_and_genotypes(tmp_path):
    pg = _pg()
    reader = pg.VcfReader(_write_vcf(tmp_path))
    assert reader.get_sample_names() == ["S1", "S2"]

    a, b = list(reader)
    assert list(a.format) == ["GT", "DP"]        # FORMAT column order
    assert [s.gt_string() for s in a.samples] == ["0/1", "1|1"]
    assert a.samples[1].phased is True          # '|' separator
    assert a.samples[0].fields["DP"] == [20]

    assert b.samples[0].is_hom_ref() is True     # 0/0
    assert b.samples[0].gt_string() == "0/0"
    assert b.samples[1].gt_string() == "./."     # missing genotype


def test_sites_only_skips_samples(tmp_path):
    pg = _pg()
    a = next(iter(pg.VcfReader(_write_vcf(tmp_path), parse_samples=False)))
    assert a.samples == []


def test_skip_filtered(tmp_path):
    pg = _pg()
    records = list(pg.VcfReader(_write_vcf(tmp_path), skip_filtered=True))
    # Only the PASS record survives.
    assert len(records) == 1
    assert records[0].start == 99


def test_header_and_contigs(tmp_path):
    pg = _pg()
    reader = pg.VcfReader(_write_vcf(tmp_path))
    assert "##fileformat=VCFv4.2" in reader.get_header()
    assert reader.get_contigs() == ["chr1"]


def test_to_coordinate_and_to_dict(tmp_path):
    pg = _pg()
    a = next(iter(pg.VcfReader(_write_vcf(tmp_path))))

    coord = a.to_coordinate()
    assert (coord.strand, coord.start, coord.end) == (".", 99, 99)  # closed, single base

    d = a.to_dict()
    assert d == {
        "id": "rs1", "ref": "A", "alt": ["G"], "qual": pytest.approx(50.0),
        "filter": ["PASS"], "is_snp": True, "is_indel": False,
    }


def test_loads_into_universal_grove(tmp_path):
    pg = _pg()
    g = pg.Grove()
    for v in pg.VcfReader(_write_vcf(tmp_path)):
        g.insert(v.chrom, v.to_coordinate(), v.to_dict())
    assert g.size() == 2
    hit = list(g.intersect(pg.GenomicCoordinate(".", 99, 99), "chr1"))[0]
    assert hit.data["id"] == "rs1"


def test_reader_accessors_after_clean_iteration(tmp_path):
    pg = _pg()
    reader = pg.VcfReader(_write_vcf(tmp_path))
    assert list(reader)  # consume both records
    # Clean EOF leaves no error; the line counter reflects records consumed.
    assert reader.get_error_message() == ""
    assert reader.get_current_line() == 2


def test_missing_file_raises(tmp_path):
    pg = _pg()
    with pytest.raises((RuntimeError, IOError, OSError)):
        pg.VcfReader(str(tmp_path / "nope.vcf"))


def test_region_filters_records(tmp_path):
    """A region string restricts iteration via htslib's tbx index (1-based
    inclusive). The fixture has variants at POS 100 and 200 on chr1."""
    pg = _pg()
    gz = _bgzip_tabix_vcf(tmp_path)
    records = list(pg.VcfReader(gz, region="chr1:150-250"))
    # start is the 0-based position; only POS 200 (start 199) overlaps.
    assert [r.start for r in records] == [199]


def test_empty_region_reads_all(tmp_path):
    """The default empty region reads the whole indexed file."""
    pg = _pg()
    gz = _bgzip_tabix_vcf(tmp_path)
    assert len(list(pg.VcfReader(gz, region=""))) == 2