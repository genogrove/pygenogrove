/*
 * json_value — the payload type for the universal `Grove`
 * (grove<genomic_coordinate, json_value>).
 *
 * It holds a JSON text string and is converted to/from arbitrary Python objects
 * (dict / list / scalar / None) via the stdlib `json` module in the pybind11
 * type_caster below. Genogrove serializes it as a length-prefixed string
 * (byte-identical to serializer<std::string>), so the resulting `.gg` is a valid
 * genogrove file readable by a C++ grove<KeyT, std::string> — the payload is just
 * JSON text rather than a typed binary record.
 */
#pragma once

#include <pybind11/pybind11.h>

#include <istream>
#include <ostream>
#include <string>

#include <genogrove/data_type/serialization_traits.hpp>

namespace pygg {

struct json_value {
    // Always valid JSON. "null" decodes to Python None, so a default / no-data
    // payload round-trips to None.
    std::string json = "null";

    void serialize(std::ostream& os) const {
        genogrove::data_type::serialization_traits<std::string>::serialize(os, json);
    }
    static json_value deserialize(std::istream& is) {
        return json_value{
            genogrove::data_type::serialization_traits<std::string>::deserialize(is)};
    }
};

}  // namespace pygg

namespace pybind11 {
namespace detail {

// Converts arbitrary JSON-serializable Python objects <-> json_value.
template <>
struct type_caster<pygg::json_value> {
    PYBIND11_TYPE_CASTER(pygg::json_value, const_name("object"));

    // Python object -> json_value (json.dumps). Propagates a TypeError if the
    // object is not JSON-serializable.
    bool load(handle src, bool) {
        object json_mod = module_::import("json");
        value.json = json_mod.attr("dumps")(src).cast<std::string>();
        return true;
    }

    // json_value -> Python object (json.loads).
    static handle cast(const pygg::json_value& v, return_value_policy /*policy*/,
                       handle /*parent*/) {
        object json_mod = module_::import("json");
        return json_mod.attr("loads")(str(v.json)).release();
    }
};

}  // namespace detail
}  // namespace pybind11