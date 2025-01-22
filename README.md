For converting geo_types::LineString into a piecewise bezier, given as a Vec\<BezierSegment\>.

A BezierSegment is just a tuple: (geo_types::Coord, Option\<geo_types::Coord\>, Option\<geo_types::Coord\>, geo_types::Coord)

The options are eiter Some or None together
If they are None the segment is a line, if they are Some the segment is a Cubic Bezier