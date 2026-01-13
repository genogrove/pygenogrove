# pygenogrove

Python bindings for the [genogrove](https://github.com/genogrove/genogrove) C++ library - a specialized B+ tree data structure optimized for genomic interval storage and querying.

## Features

- **Efficient Interval Storage**: B+ tree implementation optimized for genomic intervals
- **Multi-Index Support**: Separate trees for different chromosomes or indices
- **Fast Queries**: Overlap-based interval queries with efficient tree traversal
- **Sorted Insertion**: Optimized insertion path for pre-sorted data (10-20x speedup)
- **Memory Efficient**: Stable memory addresses with deque-based storage

## Installation

### Building from Source

Requirements:
- C++20 compatible compiler
- CMake 3.15+
- Python 3.8+

```bash
# Clone with submodules
git clone --recursive https://github.com/yourusername/pygenogrove.git
cd pygenogrove

# Build using CMake
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

# The built module will be in build/pygenogrove.so (or .pyd on Windows)
```

### Using pip (when packaged)

```bash
pip install pygenogrove
```

## Quick Start

```python
import pygenogrove as pg

# Create a grove with order 100 (max 99 keys per node)
grove = pg.Grove(100)

# Create intervals
interval1 = pg.Interval(100, 200)  # [100, 200)
interval2 = pg.Interval(150, 250)  # [150, 250)
interval3 = pg.Interval(300, 400)  # [300, 400)

# Insert intervals into different chromosomes
grove.insert("chr1", interval1)
grove.insert("chr1", interval2)
grove.insert("chr2", interval3)

print(f"Total intervals: {grove.size()}")  # Output: Total intervals: 3

# Query for overlapping intervals
query = pg.Interval(175, 225)
results = grove.intersect(query, "chr1")

print(f"Found {len(results)} overlapping intervals")
for key in results:
    interval = key.value
    print(f"  {interval.start}-{interval.end}")
```

## Usage Examples

### Basic Operations

```python
import pygenogrove as pg

# Create a grove (default order is 3)
grove = pg.Grove()

# Create and insert intervals
interval = pg.Interval(1000, 2000)
key = grove.insert("chr1", interval)

# Access interval properties
print(f"Start: {interval.start}")  # Output: Start: 1000
print(f"End: {interval.end}")      # Output: End: 2000

# Modify intervals
interval.start = 1100
interval.end = 2100
```

### Sorted Insertion (Optimized)

For pre-sorted data, use `insert_sorted()` for better performance:

```python
import pygenogrove as pg

grove = pg.Grove(100)

# Insert sorted intervals (each must be > previous)
grove.insert_sorted("chr1", pg.Interval(100, 200))
grove.insert_sorted("chr1", pg.Interval(300, 400))
grove.insert_sorted("chr1", pg.Interval(500, 600))
```

**Note**: `insert_sorted()` assumes each interval is greater than all existing intervals in that index. Using it incorrectly may corrupt the tree structure.

### Querying Intervals

```python
import pygenogrove as pg

grove = pg.Grove(100)

# Insert some intervals
grove.insert("chr1", pg.Interval(100, 200))
grove.insert("chr1", pg.Interval(300, 400))
grove.insert("chr2", pg.Interval(100, 200))

# Query specific chromosome
query = pg.Interval(150, 350)
results = grove.intersect(query, "chr1")

print(f"Found {len(results)} overlaps in chr1")
for key in results:
    print(f"  Interval: {key.value}")

# Query across all chromosomes
all_results = grove.intersect(query)
print(f"Found {len(all_results)} overlaps across all chromosomes")
```

### Overlap Detection

```python
import pygenogrove as pg

# Static method for checking overlap
interval1 = pg.Interval(100, 200)
interval2 = pg.Interval(150, 250)
interval3 = pg.Interval(300, 400)

print(pg.Interval.overlap(interval1, interval2))  # True (they overlap)
print(pg.Interval.overlap(interval1, interval3))  # False (no overlap)
```

## API Reference

### Interval

```python
Interval(start: int, end: int)
```

A genomic interval with start and end coordinates (0-based, half-open `[start, end)`).

**Attributes**:
- `start`: Start position (inclusive)
- `end`: End position (exclusive)

**Methods**:
- `Interval.overlap(a, b)`: Static method to check if two intervals overlap

### Grove

```python
Grove(order: int = 3)
```

A B+ tree container for genomic intervals with multi-index support.

**Parameters**:
- `order`: Maximum branching factor (max keys per node = order - 1)

**Methods**:
- `size()`: Get total number of intervals across all indices
- `get_order()`: Get the order (branching factor) of the tree
- `insert(index: str, interval: Interval) -> Key`: Insert an interval at the specified index
- `insert_sorted(index: str, interval: Interval) -> Key`: Insert pre-sorted interval (optimized)
- `intersect(query: Interval) -> QueryResult`: Find overlapping intervals across all indices
- `intersect(query: Interval, index: str) -> QueryResult`: Find overlapping intervals in specific index

### Key

Wrapper object for intervals stored in the grove. Returned by insert operations.

**Attributes**:
- `value`: The interval value of this key

### QueryResult

Result object containing matching intervals from a query.

**Attributes**:
- `query`: The query interval used for the search
- `keys`: List of matching keys

**Methods**:
- `__len__()`: Number of results
- `__iter__()`: Iterate over matching keys

## Current Status

This is an early development version. Currently exposed features:

- Basic grove and interval operations
- Insert and query functionality
- Multi-index support (per chromosome)
- Sorted insertion optimization

**Not yet exposed**:
- Graph overlay operations
- Genomic coordinates with strand information
- File I/O (BED, GFF/GTF readers)
- Bulk insertion operations
- Serialization/deserialization
- Data associations (key-value pairs)

## Performance Tips

1. **Choose appropriate order**: Higher order (e.g., 100-500) reduces tree height for large datasets
2. **Use sorted insertion**: If your data is pre-sorted, use `insert_sorted()` for 10-20x speedup
3. **Separate by chromosome**: Use index parameter to maintain separate trees per chromosome
4. **Query specific indices**: Query specific chromosomes instead of all indices when possible

## Development

### Building and Testing

```bash
# Build in debug mode
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build

# Run tests (when added)
cd build
ctest
```

### Project Structure

```
pygenogrove/
├── external/
│   └── genogrove/          # C++ library (git submodule)
├── src/
│   └── bindings.cpp        # pybind11 bindings
├── tests/                  # Python tests (to be added)
├── CMakeLists.txt          # Build configuration
├── README.md              # This file
└── CLAUDE.md              # Development guide
```

## Contributing

Contributions are welcome! Areas for contribution:
- Additional bindings for genogrove features
- Python tests and examples
- Documentation improvements
- Performance benchmarks

## License

This project inherits the license from the genogrove C++ library. See the LICENSE file for details.

## Related Projects

- [genogrove](https://github.com/seschwar/genogrove): The underlying C++ library

## Citation

If you use pygenogrove in your research, please cite the original genogrove library
