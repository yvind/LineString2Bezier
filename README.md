

# linestring2bezier

[![crates.io version](https://img.shields.io/crates/v/linestring2bezier.svg)](https://crates.io/crates/linestring2bezier)
[![docs.rs docs](https://docs.rs/linestring2bezier/badge.svg)](https://docs.rs/linestring2bezier)

A Rust crate for converting a `geo_types::LineString<T: CoordFloat>` into a piecewise cubic Bezier spline, and vice versa.
The conversion is controlled by a maximum allowed error, ensuring the Bezier approximation stays within a specified distance from the original line string.

## Features

- **LineString to Bezier**: Approximate a polyline as a sequence of cubic Bezier segments (`BezierString`).
- **Bezier to LineString**: Approximate a Bezier spline as a polyline within a given error.
- **Error control**: Specify the maximum deviation allowed between the original and the approximation.
- **Handles straight lines**: Segments with no curvature are represented as straight lines.

## Types

- `BezierSegment<T: CoordFloat>`: Either a cubic Bezier or a straight line.
	- For Bezier: `start`, `handle1`, `handle2`, `end` (`geo_types::Coord<T>`)
	- For Line: `start`, `end` (`geo_types::Coord<T>`)
- `BezierString<T: CoordFloat>`: A vector of `BezierSegment<T>`.

## Example

```rust
use geo_types::{Coord, LineString};
use linestring2bezier::{BezierString};

let line = LineString(vec![
		Coord { x: 0.0, y: 0.0 },
		Coord { x: 1.0, y: 2.0 },
		Coord { x: 2.0, y: 0.0 },
]);

let bezier = BezierString::from_line_string(line, 0.1).unwrap();
let approx_line = bezier.to_line_string(0.1).unwrap();
```

## Error Handling

- Returns errors if the error margin is too small or if the input is empty.

## Documentation

- [API Documentation (docs.rs)](https://docs.rs/linestring2bezier)

## License

Licensed under the MIT license ([LICENSE.md](LICENSE.md)).
