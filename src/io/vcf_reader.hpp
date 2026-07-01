/*
 * Bindings for the VCF/BCF variant reader: SampleGenotype / VcfEntry (value
 * types) and VcfReader (single-pass iterator). Mirrors genogrove/io/vcf_reader.hpp
 * and follows the BamReader binding's shape.
 *
 * vcf_entry is not serializable (it holds vectors / variant-valued maps / nested
 * sample genotypes), so there is no typed `VcfGrove`. The intended flow is to
 * load variants into the universal Grove as JSON: build the key with
 * VcfEntry.to_coordinate() and a payload with VcfEntry.to_dict() (or your own
 * dict from the full record).
 *
 * INFO and per-sample FORMAT values are htslib-typed variants
 * (bool / list[int] / list[float] / str), converted to native Python objects by
 * pybind11/stl.h's variant + map casters. Read-only — the reader yields fresh
 * copies, so there is nothing to mutate back.
 */
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // vector / unordered_map / variant casters

#include <memory>
#include <stdexcept>
#include <string>

#include <genogrove/data_type/genomic_coordinate.hpp>
#include <genogrove/io/vcf_reader.hpp>

namespace py = pybind11;
namespace gio = genogrove::io;
namespace gdt = genogrove::data_type;

inline void bind_vcf_reader(py::module_& m) {
    // ---- Per-sample genotype + FORMAT data ----
    py::class_<gio::sample_genotype>(m, "SampleGenotype", R"pbdoc(
        Genotype and FORMAT data for one sample at one variant record.
    )pbdoc")
        .def(py::init<>())
        .def_readonly("gt_alleles", &gio::sample_genotype::gt_alleles,
                      "Decoded GT allele indices (0 = REF, 1.. = ALT index, "
                      "-1 = missing '.'); empty when no GT field.")
        .def_readonly("phased", &gio::sample_genotype::phased,
                      "Whether the genotype is phased (the '|' separator).")
        .def_readonly("has_gt", &gio::sample_genotype::has_gt,
                      "True if a GT field was present for this sample.")
        .def_readonly("fields", &gio::sample_genotype::fields,
                      "Other FORMAT fields keyed by tag (e.g. 'DP', 'GQ'), "
                      "excluding GT. Values are list[int] / list[float] / str.")
        .def("gt_string", &gio::sample_genotype::gt_string,
             "GT as a VCF-style string ('0/1', '0|1', './.', '' if no GT).")
        .def("is_hom_ref", &gio::sample_genotype::is_hom_ref,
             "True if every called allele is the reference (all 0); a missing "
             "genotype is not hom-ref.")
        .def("__repr__", [](const gio::sample_genotype& s) {
            return "SampleGenotype(gt='" + s.gt_string() + "')";
        });

    // ---- A single VCF/BCF record ----
    py::class_<gio::vcf_entry>(m, "VcfEntry", R"pbdoc(
        A single VCF/BCF variant record (0-based half-open `start`/`end`,
        decoupled from any grove key type).

        Load variants into the universal Grove as JSON, deriving the key from the
        record::

            g = pygenogrove.Grove()
            for v in pygenogrove.VcfReader("calls.vcf"):
                g.insert(v.chrom, v.to_coordinate(), v.to_dict())
    )pbdoc")
        .def(py::init<>())
        .def_readonly("chrom", &gio::vcf_entry::chrom, "CHROM (contig name).")
        .def_readonly("start", &gio::vcf_entry::start,
                      "0-based start (htslib position; VCF POS - 1).")
        .def_readonly("end", &gio::vcf_entry::end,
                      "0-based exclusive end (start + len(REF)).")
        .def_readonly("id", &gio::vcf_entry::id, "ID, empty when '.'.")
        .def_readonly("ref", &gio::vcf_entry::ref, "REF allele.")
        .def_readonly("alt", &gio::vcf_entry::alt,
                      "ALT alleles (empty for monomorphic records).")
        .def_readonly("qual", &gio::vcf_entry::qual,
                      "QUAL score (only valid when qual_missing is False).")
        .def_readonly("qual_missing", &gio::vcf_entry::qual_missing,
                      "True when QUAL is '.' (missing).")
        .def_readonly("filter", &gio::vcf_entry::filter,
                      "FILTER entries; ['PASS'] when passed, empty when '.'.")
        .def_readonly("info", &gio::vcf_entry::info,
                      "INFO fields (when parse_info): name -> "
                      "bool / list[int] / list[float] / str.")
        .def_readonly("format", &gio::vcf_entry::format,
                      "FORMAT keys in column order (incl. 'GT'); empty unless "
                      "parse_samples.")
        .def_readonly("samples", &gio::vcf_entry::samples,
                      "Per-sample SampleGenotype list, parallel to "
                      "VcfReader.get_sample_names().")
        .def("passed_filter", &gio::vcf_entry::passed_filter,
             "True if FILTER is PASS or unset ('.').")
        .def("is_snp", &gio::vcf_entry::is_snp,
             "True for a simple biallelic single-base SNP.")
        .def("is_indel", &gio::vcf_entry::is_indel,
             "True if any sequence ALT differs in length from REF.")
        .def_static("is_symbolic_allele", &gio::vcf_entry::is_symbolic_allele,
                    py::arg("allele"),
                    "True for symbolic / non-sequence ALT alleles "
                    "(<DEL>, *, ., breakends).")
        .def("to_coordinate",
             [](const gio::vcf_entry& e) {
                 if (e.end <= e.start) {
                     throw std::invalid_argument(
                         "VcfEntry spans no reference bases; cannot derive a "
                         "coordinate");
                 }
                 // VCF has no strand -> '.'; half-open [start, end) -> closed
                 // [start, end - 1].
                 return gdt::genomic_coordinate('.', e.start, e.end - 1);
             },
             R"pbdoc(
                 Derive the GenomicCoordinate key for this variant (unstranded
                 '.'; 0-based half-open [start, end) -> closed [start, end-1]).
                 Raises ValueError if the record spans no reference bases.
             )pbdoc")
        .def("to_dict",
             [](const gio::vcf_entry& e) {
                 py::dict d;
                 d["id"] = e.id;
                 d["ref"] = e.ref;
                 d["alt"] = e.alt;
                 d["qual"] = e.qual_missing ? py::object(py::none())
                                            : py::cast(e.qual);
                 d["filter"] = e.filter;
                 d["is_snp"] = e.is_snp();
                 d["is_indel"] = e.is_indel();
                 return d;
             },
             "A JSON-serializable dict of the core variant fields (id, ref, alt, "
             "qual, filter, is_snp, is_indel) — convenient as a Grove payload.")
        .def("__repr__", [](const gio::vcf_entry& e) {
            std::string alt = e.alt.empty() ? "." : e.alt.front();
            if (e.alt.size() > 1) alt += ",...";
            return "VcfEntry(chrom='" + e.chrom + "', start=" +
                   std::to_string(e.start) + ", ref='" + e.ref + "', alt='" +
                   alt + "')";
        });

    // ---- The reader ----
    py::class_<gio::vcf_reader>(m, "VcfReader", R"pbdoc(
        Single-pass iterator over a VCF/BCF file (plain, bgzip-ed, or binary BCF;
        htslib auto-detects). Yields VcfEntry. Not thread-safe — drive one reader
        per thread.

            for v in pygenogrove.VcfReader("calls.vcf", skip_filtered=True):
                ...
    )pbdoc")
        .def(py::init([](const std::string& path, bool parse_info,
                         bool parse_samples, bool skip_filtered,
                         const std::string& region) {
                 gio::vcf_reader_options opts;
                 opts.parse_info = parse_info;
                 opts.parse_samples = parse_samples;
                 opts.skip_filtered = skip_filtered;
                 opts.region = region;
                 return std::make_unique<gio::vcf_reader>(path, opts);
             }),
             py::arg("path"), py::arg("parse_info") = true,
             py::arg("parse_samples") = true, py::arg("skip_filtered") = false,
             py::arg("region") = "",
             "Open a VCF/BCF file. parse_info / parse_samples toggle INFO and "
             "per-sample parsing; skip_filtered drops non-PASS records. region "
             "is an htslib region string (\"chr:start-end\", 1-based inclusive); "
             "when set, only overlapping records are yielded and a CSI/TBI-indexed "
             "bgzip VCF or a BCF is required. Empty (default) reads the whole file.")
        .def("__iter__",
             [](gio::vcf_reader& r) -> gio::vcf_reader& { return r; })
        .def("__next__",
             [](gio::vcf_reader& r) {
                 gio::vcf_entry entry;
                 if (!r.read_next(entry)) {
                     throw py::stop_iteration();
                 }
                 return entry;
             },
             // htslib decode / parse touches no Python objects; the GIL is
             // reacquired before the returned entry is converted.
             py::call_guard<py::gil_scoped_release>())
        .def("get_header", &gio::vcf_reader::get_header,
             "Full VCF header text (## meta lines + the #CHROM column line).")
        .def("get_sample_names", &gio::vcf_reader::get_sample_names,
             "Sample names in column order (empty for sites-only VCFs).")
        .def("get_contigs", &gio::vcf_reader::get_contigs,
             "Contig names declared in the header.")
        .def("get_error_message", &gio::vcf_reader::get_error_message,
             "Error message from the most recent read; empty on clean EOF.")
        .def("get_current_line", &gio::vcf_reader::get_current_line,
             "1-based index of the most recently consumed record (counts records "
             "dropped by skip_filtered too); 0 before the first read.");
}