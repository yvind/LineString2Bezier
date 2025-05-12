For converting `geo_types::LineString<T: CoordFloat>` into a piecewise cubic bezier, given as a `Vec<BezierSegment<T>>`.  

`BezierSegment<T: CoordFloat>` has fields:  
 start: `geo_types::Coord<T>`  
 handles: `Option<(geo_types::Coord<T>, geo_types::Coord<T>)>`  
 end: `geo_types::Coord<T>`  

The handles field is `None` if the segment is a straight line, otherwise the segment is a Cubic Bezier.