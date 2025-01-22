For converting geo LineString into a piecewise bezier, given as a Vec\<BezierSegment\>.

A BezierSegment is just a tuple: (geo::Coord, Option\<geo::Coord\>, Option\<geo::Coord\>, geo::Coord)

The options are eiter Some or None together
If they are None the segment is a line, if they are Some the segment is a Cubic Bezier