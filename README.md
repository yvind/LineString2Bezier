For converting geo_types::LineString into a piecewise bezier, given as a `Vec<BezierSegment>`.

A BezierSegment is just a tuple: `(geo_types::Coord, Option<(geo_types::Coord, geo_types::Coord)>, geo_types::Coord)`

The option is `None` if the segment is a straight line, otherwise the segment is a Cubic Bezier