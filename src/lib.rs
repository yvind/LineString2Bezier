//! Approximate a [LineString] with a cubic bezier spline [BezierString]
//! or approximate a [BezierString] with a [LineString]
//!
//! Pass the line_string and a maximum allowed error to [BezierString::from_line_string]
//! The error is the maximum distance between the line_string and the bezier spline
//! evaluated at the line_string vertices and at the middle point between the vertices

#![deny(
    elided_lifetimes_in_paths,
    explicit_outlives_requirements,
    keyword_idents,
    macro_use_extern_crate,
    meta_variable_misuse,
    missing_abi,
    missing_debug_implementations,
    missing_docs,
    non_ascii_idents,
    noop_method_call,
    rust_2021_incompatible_closure_captures,
    rust_2021_incompatible_or_patterns,
    rust_2021_prefixes_incompatible_syntax,
    rust_2024_prelude_collisions,
    single_use_lifetimes,
    trivial_casts,
    trivial_numeric_casts,
    unreachable_pub,
    unsafe_code,
    unsafe_op_in_unsafe_fn,
    unused_crate_dependencies,
    unused_extern_crates,
    unused_import_braces,
    unused_lifetimes,
    unused_qualifications,
    unused_results,
    clippy::unwrap_used,
    warnings
)]

use geo_types::{Coord, CoordFloat, Line, LineString};

type Result<T> = std::result::Result<T, Error>;

/// Error enum for the crate
#[derive(Debug, Clone, thiserror::Error)]
pub enum Error {
    /// The given error margin is too small
    #[error("The given error margin is too small")]
    TooSmallErrorGiven,
    /// Tried to convert an empty [BezierString] to [LineString]
    #[error("An empty BezierString cannot be made into a LineString")]
    EmptyBezierString,
    /// Tried to convert an empty [LineString] to [BezierString]
    #[error("An empty LineString cannot be made into a BezierString")]
    EmptyLineString,
}

/// Simple struct to represent a cubic bezier curve
#[derive(Clone, Debug)]
pub struct BezierCurve<T: CoordFloat = f64> {
    /// The start coord for the bezier curve
    pub start: Coord<T>,
    /// The first handle
    pub handle1: Coord<T>,
    /// The second handle
    pub handle2: Coord<T>,
    /// The end coord for the bezier curve
    pub end: Coord<T>,
}

impl<T: CoordFloat> BezierCurve<T> {
    /// Returns a new [BezierCurve]
    pub fn new(start: Coord<T>, handle1: Coord<T>, handle2: Coord<T>, end: Coord<T>) -> Self {
        BezierCurve {
            start,
            handle1,
            handle2,
            end,
        }
    }

    /// Returns a new [BezierCurve] with all fields set to [Coord::zero]
    pub fn zero() -> Self {
        BezierCurve {
            start: Coord::zero(),
            handle1: Coord::zero(),
            handle2: Coord::zero(),
            end: Coord::zero(),
        }
    }

    /// Get the bezier as an array of the coords
    pub fn to_array(self) -> [Coord<T>; 4] {
        [self.start, self.handle1, self.handle2, self.end]
    }

    /// Approximate this cubic bezier with a LineString with a max error of `error``
    pub fn to_line_string(self, error: T) -> Result<LineString<T>> {
        if error <= T::epsilon() * (T::one() + T::one()) {
            return Err(Error::TooSmallErrorGiven);
        }

        let mut num_sample_points = 3;

        let mut coords = Vec::with_capacity(2 * num_sample_points);
        let mut ts = Vec::with_capacity(2 * num_sample_points);

        let arr_bc = self.clone().to_array();
        let slice = arr_bc.as_slice();

        loop {
            ts.clear();
            coords.clear();

            #[allow(clippy::unwrap_used)]
            let sample_interval = T::one() / T::from(num_sample_points - 1).unwrap();

            coords.push(self.start);
            ts.push(T::zero());
            for i in 1..num_sample_points - 1 {
                #[allow(clippy::unwrap_used)]
                let t = T::from(i).unwrap() * sample_interval;

                ts.push(t);

                let coord = self.start * bezier_basis_0(t)
                    + self.handle1 * bezier_basis_1(t)
                    + self.handle2 * bezier_basis_2(t)
                    + self.end * bezier_basis_3(t);

                coords.push(coord);
            }
            ts.push(T::one());
            coords.push(self.end);

            let (mut current_error, _) =
                compute_max_error(&coords, 0, coords.len() - 1, slice, &ts);
            if current_error > error {
                // let's try a newton step
                reparameterize(&mut ts, &coords, 0, coords.len() - 1, slice);
                (current_error, _) = compute_max_error(&coords, 0, coords.len() - 1, slice, &ts);
            }
            if current_error <= error {
                break;
            }
            num_sample_points += 1;
        }

        Ok(LineString(coords))
    }
}

/// A BezierSegment is either a straight [Line] or a [BezierCurve]
#[derive(Debug, Clone)]
pub enum BezierSegment<T: CoordFloat = f64> {
    /// The segment is a cubic bezier
    Bezier(BezierCurve<T>),
    /// The segment is a line
    Line(Line<T>),
}

impl<T: CoordFloat> BezierSegment<T> {
    /// Is this segment a cubic bezier?
    pub fn is_bezier_curve(&self) -> bool {
        matches!(self, BezierSegment::Bezier(_))
    }

    /// Construct a Bezier segment given the start, end and optionally two handles
    pub fn new(start: Coord<T>, handles: Option<(Coord<T>, Coord<T>)>, end: Coord<T>) -> Self {
        if let Some((handle1, handle2)) = handles {
            BezierSegment::Bezier(BezierCurve {
                start,
                handle1,
                handle2,
                end,
            })
        } else {
            BezierSegment::Line(Line::new(start, end))
        }
    }

    /// Get an approximation of the euclidean length of the bezier segment
    /// For the cubic bezier this is taken as the mean of the
    /// distance between the end points and the distance
    /// from start to end through the handles
    /// For the line it is simple the euclidean line length
    pub fn length_approximation(&self) -> T {
        match self {
            BezierSegment::Bezier(bc) => {
                ((bc.end - bc.handle2).magnitude()
                    + (bc.handle2 - bc.handle1).magnitude()
                    + (bc.handle1 - bc.start).magnitude()
                    + (bc.end - bc.start).magnitude())
                    / (T::one() + T::one())
            }
            BezierSegment::Line(line) => line.delta().magnitude(),
        }
    }

    /// Get a LineString with at most error deviation from the Bezier segment
    /// The deviation is measured at the mid point between neighboring nodes
    pub fn to_line_string(self, error: T) -> Result<LineString<T>> {
        match self {
            BezierSegment::Line(line) => Ok(LineString::from(line)),
            BezierSegment::Bezier(bc) => bc.to_line_string(error),
        }
    }
}

impl<T> From<[Coord<T>; 4]> for BezierSegment<T>
where
    T: CoordFloat,
{
    fn from(value: [Coord<T>; 4]) -> Self {
        BezierSegment::Bezier(BezierCurve {
            start: value[0],
            handle1: value[1],
            handle2: value[2],
            end: value[3],
        })
    }
}

impl<T> From<[Coord<T>; 2]> for BezierSegment<T>
where
    T: CoordFloat,
{
    fn from(value: [Coord<T>; 2]) -> Self {
        BezierSegment::Line(Line {
            start: value[0],
            end: value[1],
        })
    }
}

// Bezier basis functions
#[inline]
fn bezier_basis_0<T: CoordFloat>(t: T) -> T {
    (T::one() - t).powi(3)
}

#[inline]
fn bezier_basis_1<T: CoordFloat>(t: T) -> T {
    (T::one() + T::one() + T::one()) * t * (T::one() - t).powi(2)
}

#[inline]
fn bezier_basis_2<T: CoordFloat>(t: T) -> T {
    (T::one() + T::one() + T::one()) * t.powi(2) * (T::one() - t)
}

#[inline]
fn bezier_basis_3<T: CoordFloat>(t: T) -> T {
    t.powi(3)
}

/// A BezierString is simply a vector of [BezierSegment]s
#[derive(Debug, Clone)]
pub struct BezierString<T: CoordFloat = f64>(pub Vec<BezierSegment<T>>);

impl<T: CoordFloat> BezierString<T> {
    /// Returns a [BezierString] with the given segments
    pub fn new(segments: Vec<BezierSegment<T>>) -> Self {
        BezierString(segments)
    }

    /// Returns an empty [BezierString]
    pub fn empty() -> Self {
        BezierString(Vec::new())
    }

    /// Return an iterator yielding the segments
    pub fn segments(&self) -> std::slice::Iter<'_, BezierSegment<T>> {
        self.0.iter()
    }

    /// Return an iterator yielding the segments as mutable
    pub fn segments_mut(&mut self) -> std::slice::IterMut<'_, BezierSegment<T>> {
        self.0.iter_mut()
    }

    /// The number of vertices and handles in this [BezierString]
    /// The endpoints of each [BezierSegment] are not counted double
    pub fn num_points(&self) -> usize {
        let mut num_points = 0;

        for segment in self.0.iter() {
            if segment.is_bezier_curve() {
                num_points += 3;
            } else {
                num_points += 1;
            }
        }
        num_points + 1
    }

    /// Get the number of segments that make up this [BezierString]
    pub fn num_segments(&self) -> usize {
        self.0.len()
    }

    /// Return the coordinates of a [BezierString] as a Vec of [BezierSegment]s
    pub fn into_inner(self) -> Vec<BezierSegment<T>> {
        self.0
    }

    /// Get a length approximation of the [BezierString]
    pub fn length_approximation(&self) -> T {
        self.0
            .iter()
            .fold(T::zero(), |acc, e| acc + e.length_approximation())
    }

    /// Convert a [BezierString] to a [LineString] with a maximum of `error` deviation between the two
    pub fn to_line_string(self, error: T) -> Result<LineString<T>> {
        let mut line = LineString::new(Vec::with_capacity(self.num_points()));

        let mut segments = self.0.into_iter();
        {
            let first_segment = segments.next().ok_or(Error::EmptyBezierString)?;
            let iter = first_segment.to_line_string(error)?.0.into_iter();
            line.0.extend(iter);
        }
        for segment in segments {
            let iter = segment.to_line_string(error)?.0.into_iter();
            line.0.extend(iter.skip(1));
        }

        Ok(line)
    }

    /// Convert a [LineString] to a [BezierString] with a maximum of `error` deviation between the two at the line_string vertices
    pub fn from_line_string(line_string: LineString<T>, error: T) -> Result<BezierString<T>> {
        if error <= T::epsilon() * (T::one() + T::one()) {
            return Err(Error::TooSmallErrorGiven);
        }
        let n_pts = line_string.0.len();
        if n_pts < 2 {
            return Err(Error::EmptyLineString);
        }
        if n_pts == 2 || (n_pts == 3 && line_string.is_closed()) {
            let mut segments = Vec::with_capacity(n_pts - 1);

            let mut prev_coord = line_string.0[0];
            for p in line_string.0.into_iter().skip(1) {
                segments.push(BezierSegment::new(prev_coord, None, p));
                prev_coord = p;
            }

            return Ok(BezierString(segments));
        }

        let mut tangent_right =
            compute_right_tangent(line_string.0.as_slice(), 0).unwrap_or(Coord::zero());
        let mut tangent_left =
            compute_left_tangent(line_string.0.as_slice(), n_pts - 1).unwrap_or(Coord::zero());

        if line_string.is_closed() {
            let diff = tangent_right - tangent_left;
            if let Some(norm_diff) = diff.try_normalize() {
                tangent_right = norm_diff;
                tangent_left = -tangent_right;
            }
        }

        let mut bezier_segments = Vec::with_capacity(5);

        // recursively tries to fit bezier segments to the line_string
        fit_cubic(
            line_string.0.as_slice(),
            0,
            n_pts - 1,
            tangent_right,
            tangent_left,
            error,
            &mut bezier_segments,
        );
        Ok(BezierString(bezier_segments))
    }
}

// recursive function to fit bezier curve segments to line_string
fn fit_cubic<T: CoordFloat>(
    polyline: &[Coord<T>],
    first: usize,
    last: usize,
    tangent_start: Coord<T>,
    tangent_end: Coord<T>,
    error: T,
    bezier_string: &mut Vec<BezierSegment<T>>,
) {
    // Handle two-point case, recursion base case
    if last - first == 1 {
        bezier_string.push(BezierSegment::new(polyline[first], None, polyline[last]));
        return;
    }

    // Parameterize points and attempt to fit curve,
    let mut ts = chord_length_parameterize(polyline, first, last);
    let mut bez_curve = generate_bezier(
        polyline,
        first,
        last,
        ts.as_slice(),
        &tangent_start,
        &tangent_end,
    );

    // Find max deviation
    let (max_error, _) = compute_max_error(polyline, first, last, &bez_curve, &ts);

    if max_error < error {
        bezier_string.push(bez_curve.into());
        return;
    }

    // Let's try one step of Newton-Rhapson
    reparameterize(&mut ts, polyline, first, last, &bez_curve);
    bez_curve = generate_bezier(polyline, first, last, &ts, &tangent_start, &tangent_end);
    let (max_error, split_point) = compute_max_error(polyline, first, last, &bez_curve, &ts);

    if max_error < error {
        bezier_string.push(bez_curve.into());
        return;
    }

    // Fitting failed, split at the point of max error and fit each part recursively
    let tangent_center = match compute_center_tangent(polyline, split_point) {
        Some(t) => t,
        None => {
            let rt = compute_right_tangent(polyline, split_point);
            let lt = compute_left_tangent(polyline, split_point);

            // computing the center tangent failed => the point before and after the split point is the same
            // or some infinite coords are encountered
            // lets assign a tangent perpendicular to rt and lt if they are non-zero
            // otherwise we have three identical points in a row and might as well assign zero to the tangent
            if let Some(rt) = rt
                && let Some(lt) = lt
            {
                // rt and lt must be equal
                (rt + lt)
                    .try_normalize()
                    .map(|c| Coord { x: -c.y, y: c.x })
                    .unwrap_or(Coord::zero())
            } else {
                Coord::zero()
            }
        }
    };

    let tangent_center_neg = -tangent_center;
    fit_cubic(
        polyline,
        first,
        split_point,
        tangent_start,
        tangent_center_neg,
        error,
        bezier_string,
    );
    fit_cubic(
        polyline,
        split_point,
        last,
        tangent_center,
        tangent_end,
        error,
        bezier_string,
    );
}

// fit a bezier-segment to the polyline between first and last
// using least-squares method to find the Bezier handles for region
// with t-values for evaluation and end point tangents given
fn generate_bezier<T: CoordFloat>(
    polyline: &[Coord<T>],
    first: usize,
    last: usize,
    ts: &[T],
    start_tangent: &Coord<T>,
    end_tangent: &Coord<T>,
) -> [Coord<T>; 4] {
    let n_pts = last - first + 1;
    let mut a: Vec<[Coord<T>; 2]> = Vec::with_capacity(n_pts);

    // Compute the A's
    for &t in ts {
        a.push([
            *start_tangent * bezier_basis_1(t),
            *end_tangent * bezier_basis_2(t),
        ]);
    }

    // Create the C and X matrices
    // C is symmetric 2x2, sum of indices in 2x2 gives index in flat array
    // X is a 2x1 vector
    let mut c = [T::zero(); 3];
    let mut x = [T::zero(); 2];

    for (i, &t) in ts.iter().enumerate() {
        c[0] = c[0] + a[i][0].dot_product(&a[i][0]);
        c[1] = c[1] + a[i][0].dot_product(&a[i][1]);
        c[2] = c[2] + a[i][1].dot_product(&a[i][1]);

        let tmp = polyline[first + i]
            - (polyline[first] * (bezier_basis_0(t) + bezier_basis_1(t))
                + polyline[last] * (bezier_basis_2(t) + bezier_basis_3(t)));

        x[0] = x[0] + a[i][0].dot_product(&tmp);
        x[1] = x[1] + a[i][1].dot_product(&tmp);
    }

    // Compute the determinants
    let det_c0_c1 = c[0] * c[2] - c[1] * c[1];
    let det_c0_x = c[0] * x[1] - c[1] * x[0];
    let det_x_c1 = x[0] * c[2] - x[1] * c[1];

    // Derive alpha values
    let (mut alpha_l, mut alpha_r) = if det_c0_c1.abs() < T::epsilon() {
        (T::zero(), T::zero())
    } else {
        (det_x_c1 / det_c0_c1, det_c0_x / det_c0_c1)
    };

    // If alpha negative, use the Wu/Barsky heuristic
    let seg_length = (polyline[last] - polyline[first]).magnitude();
    #[allow(clippy::unwrap_used)]
    let epsilon = T::from(1.0e-6).unwrap() * seg_length;
    if alpha_l < epsilon || alpha_r < epsilon {
        let dist = seg_length / (T::one() + T::one() + T::one());
        alpha_l = dist;
        alpha_r = dist;
    }

    [
        polyline[first],
        polyline[first] + *start_tangent * alpha_l,
        polyline[last] + *end_tangent * alpha_r,
        polyline[last],
    ]
}

// Vertex tangent functions
#[inline]
fn compute_right_tangent<T: CoordFloat>(polyline: &[Coord<T>], end: usize) -> Option<Coord<T>> {
    (polyline[end + 1] - polyline[end]).try_normalize()
}

#[inline]
fn compute_left_tangent<T: CoordFloat>(polyline: &[Coord<T>], end: usize) -> Option<Coord<T>> {
    (polyline[end - 1] - polyline[end]).try_normalize()
}

#[inline]
fn compute_center_tangent<T: CoordFloat>(polyline: &[Coord<T>], center: usize) -> Option<Coord<T>> {
    (polyline[center + 1] - polyline[center - 1]).try_normalize()
}

// normalized length along line_string from start of segment to every vertex in segment
// the t values to be used for computing error of fitted bezier
fn chord_length_parameterize<T: CoordFloat>(
    polyline: &[Coord<T>],
    first: usize,
    last: usize,
) -> Vec<T> {
    let mut ts = Vec::with_capacity(last - first + 1);

    ts.push(T::zero());
    for i in first..last {
        ts.push(ts[i - first] + (polyline[i + 1] - polyline[i]).magnitude());
    }

    let t_last = ts[last - first];
    for t in ts.iter_mut() {
        *t = *t / t_last;
    }
    ts
}

// given a set of points and their parameterization on the bez curve
// use Newton-Rhapson to try and refine the parameterization
fn reparameterize<T: CoordFloat>(
    ts: &mut [T],
    polyline: &[Coord<T>],
    first: usize,
    last: usize,
    bez_curve: &[Coord<T>],
) {
    for i in first..=last {
        ts[i - first] = newton_raphson(bez_curve, polyline[i], ts[i - first]);
    }
}

// bez_curve Q(u) at time t is supposed to be p, refine t
fn newton_raphson<T: CoordFloat>(bez_curve: &[Coord<T>], p: Coord<T>, t: T) -> T {
    let two = T::one() + T::one();
    let three = two + T::one();

    // Q(t)
    let bez_t = evaluate_bezier(3, bez_curve, t);

    // Cubic bez prime is quadratic
    let mut bez_prime = [Coord::<T>::zero(); 3];
    // Cubic bez double prime is linear
    let mut bez_double_prime = [Coord::<T>::zero(); 2];

    // Generate control vertices for Q'
    for i in 0..3 {
        bez_prime[i].x = (bez_curve[i + 1].x - bez_curve[i].x) * three;
        bez_prime[i].y = (bez_curve[i + 1].y - bez_curve[i].y) * three;
    }

    // Generate control vertices for Q''
    for i in 0..2 {
        bez_double_prime[i].x = (bez_prime[i + 1].x - bez_prime[i].x) * two;
        bez_double_prime[i].y = (bez_prime[i + 1].y - bez_prime[i].y) * two;
    }

    // Compute Q'(t) and Q''(t)
    let qp_t = evaluate_bezier(2, &bez_prime, t);
    let qpp_t = evaluate_bezier(1, &bez_double_prime, t);

    // f(t) is Q linearized around the root Q(t) - p == 0

    // Compute f'(t)
    let denominator = (qp_t.x) * (qp_t.x)
        + (qp_t.y) * (qp_t.y)
        + (bez_t.x - p.x) * (qpp_t.x)
        + (bez_t.y - p.y) * (qpp_t.y);
    if denominator.abs() <= two * T::epsilon() {
        return t;
    }
    // Compute f(t)
    let numerator = (bez_t.x - p.x) * (qp_t.x) + (bez_t.y - p.y) * (qp_t.y);

    // t = t - f(t)/f'(t)
    t - (numerator / denominator)
}

// max distance between polyline vertices and the fitted curve and index of that vertex
fn compute_max_error<T: CoordFloat>(
    polyline: &[Coord<T>],
    first: usize,
    last: usize,
    bez_curve: &[Coord<T>],
    ts: &[T],
) -> (T, usize) {
    //let two = T::one() + T::one();

    let mut split_point = first + 1;
    let mut max_dist = T::zero();
    for i in first..=last {
        let p = evaluate_bezier(3, bez_curve, ts[i - first]);
        let dist = (p - polyline[i]).magnitude_squared();

        if dist >= max_dist {
            max_dist = dist;
            split_point = i;
        }

        // let p = evaluate_bezier(3, bez_curve, (ts[i - first] + ts[i - first + 1]) / two);
        // let dist = (p - (polyline[i] + polyline[i + 1]) / two).magnitude_squared();
        // if dist >= max_dist {
        //     max_dist = dist;
        //     split_point = i.max(first + 1);
        // }
    }
    (max_dist.sqrt(), split_point)
}

// evaluate a bezier of the given degree at time t
fn evaluate_bezier<T: CoordFloat>(degree: usize, bezier_segment: &[Coord<T>], t: T) -> Coord<T> {
    // Create a temporary vector to store the control points
    let mut v_temp = bezier_segment[..=degree].to_vec();

    // De Casteljau algorithm, just lerp-ing between lower degree beziers
    for i in 1..=degree {
        for j in 0..=(degree - i) {
            v_temp[j] = v_temp[j] * (T::one() - t) + v_temp[j + 1] * t;
        }
    }
    v_temp[0]
}

// impl the used geo functions myself to avoid
// having the entire geo crate as a dependency
trait VectorTraits<T = f64, Rhs = Self>
where
    Self: Sized,
    T: CoordFloat,
{
    fn magnitude_squared(&self) -> T;
    fn magnitude(&self) -> T;
    fn dot_product(&self, other: &Rhs) -> T;
    fn try_normalize(&self) -> Option<Self>;
    fn is_finite(&self) -> bool;
}

impl<T: CoordFloat> VectorTraits<T> for Coord<T> {
    #[inline]
    fn magnitude_squared(&self) -> T {
        self.x * self.x + self.y * self.y
    }

    #[inline]
    fn magnitude(&self) -> T {
        self.magnitude_squared().sqrt()
    }

    #[inline]
    fn dot_product(&self, other: &Self) -> T {
        self.x * other.x + self.y * other.y
    }

    fn try_normalize(&self) -> Option<Self> {
        let magnitude = self.magnitude();
        let result = *self / magnitude;

        if magnitude.is_finite() && result.is_finite() {
            Some(result)
        } else {
            None
        }
    }

    #[inline]
    fn is_finite(&self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }
}

#[cfg(test)]
mod tests {
    use crate::{BezierCurve, BezierSegment, BezierString};
    use geo::{Distance, Euclidean};
    use geo_types::{Coord, LineString, coord, line_string};

    fn actual_max_dist_between_bezier_and_line_string<T: geo_types::CoordFloat>(
        curve: &BezierCurve<T>,
        line_string: &LineString<T>,
        inner_sample_points: usize,
    ) -> T {
        let sample_points = inner_sample_points + 1;
        let ts = (0..=sample_points)
            .map(|e| T::from(e).expect("always ok") / T::from(sample_points).expect("always ok"))
            .collect::<Vec<_>>();
        let curve_array = curve.clone().to_array();
        let mut max_dist = T::zero();
        for t in ts {
            let bc = crate::evaluate_bezier(3, curve_array.as_slice(), t);
            let dist = Euclidean.distance(line_string, &geo_types::Point::from(bc));
            if dist > max_dist {
                max_dist = dist;
            }
        }
        max_dist
    }

    fn line_string_f64() -> LineString {
        line_string![
            coord! {x: -61.026, y: 23.042},
            coord! {x: -59.296, y: -8.572},
            coord! {x: -26.109, y: -21.155},
            coord! {x: 5.977, y: -6.370},
            coord! {x: -31.772, y: 26.031},
            coord! {x: 7.235, y: 49.309},
            coord! {x: 36.647, y: 16.908},
            coord! {x: 21.548, y: 12.347},
            coord! {x: 26.581, y: -40.186},
            coord! {x: -16.043, y: -38.142},
            coord! {x: -61.656, y: -28.547},
            coord! {x: -72.823, y: -46.320},
            coord! {x: -85.248, y: 6.685},
            coord! {x: -77.384, y: 21.469},
            coord! {x: -68.733, y: 27.289},
            coord! {x: -61.026, y: 23.042},
        ]
    }

    fn bezier_string_f64() -> BezierString {
        BezierString(vec![
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: -61.026, y: 23.042 },
                handle1: coord! { x: -57.879, y: 11.085 },
                handle2: coord! { x: -66.961, y: 1.130 },
                end: coord! { x: -59.296, y: -8.572 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: -59.296, y: -8.572 },
                handle1: coord! { x: -50.704, y: -19.448 },
                handle2: coord! { x: -39.962, y: -21.622 },
                end: coord! { x: -26.109, y: -21.155 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: -26.109, y: -21.155 },
                handle1: coord! { x: -12.320, y: -20.690 },
                handle2: coord! { x: 7.621, y: -20.068 },
                end: coord! { x: 5.977, y: -6.370 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: 5.977, y: -6.370 },
                handle1: coord! { x: 3.662, y: 12.919 },
                handle2: coord! { x: -32.211, y: 6.608 },
                end: coord! { x: -31.772, y: 26.031 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: -31.772, y: 26.031 },
                handle1: coord! { x: -31.371, y: 43.766 },
                handle2: coord! { x: -10.349, y: 51.654 },
                end: coord! { x: 7.235, y: 49.309 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: 7.235, y: 49.309 },
                handle1: coord! { x: 24.174, y: 47.050 },
                handle2: coord! { x: 30.476, y: 32.844 },
                end: coord! { x: 36.647, y: 16.908 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: 36.647, y: 16.908 },
                handle1: coord! { x: 38.871, y: 11.164 },
                handle2: coord! { x: 22.617, y: 18.413 },
                end: coord! { x: 21.548, y: 12.347 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: 21.548, y: 12.347 },
                handle1: coord! { x: 17.970, y: -7.949 },
                handle2: coord! { x: 38.889, y: -23.655 },
                end: coord! { x: 26.581, y: -40.186 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: 26.581, y: -40.186 },
                handle1: coord! { x: 16.629, y: -53.553 },
                handle2: coord! { x: 0.479, y: -40.321 },
                end: coord! { x: -16.043, y: -38.142 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: -16.043, y: -38.142 },
                handle1: coord! { x: -34.090, y: -35.762 },
                handle2: coord! { x: -43.639, y: -25.952 },
                end: coord! { x: -61.656, y: -28.547 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: -61.656, y: -28.547 },
                handle1: coord! { x: -69.769, y: -29.716 },
                handle2: coord! { x: -68.262, y: -53.131 },
                end: coord! { x: -72.823, y: -46.320 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: -72.823, y: -46.320 },
                handle1: coord! { x: -84.652, y: -28.654 },
                handle2: coord! { x: -83.821, y: -14.528 },
                end: coord! { x: -85.248, y: 6.685 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: -85.248, y: 6.685 },
                handle1: coord! { x: -85.687, y: 13.210 },
                handle2: coord! { x: -81.474, y: 16.366 },
                end: coord! { x: -77.384, y: 21.469 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: -77.384, y: 21.469 },
                handle1: coord! { x: -74.837, y: 24.646 },
                handle2: coord! { x: -72.786, y: 26.899 },
                end: coord! { x: -68.733, y: 27.289 },
            }),
            BezierSegment::Bezier(BezierCurve {
                start: coord! { x: -68.733, y: 27.289 },
                handle1: coord! { x: -65.312, y: 27.618 },
                handle2: coord! { x: -61.901, y: 26.365 },
                end: coord! { x: -61.026, y: 23.042 },
            }),
        ])
    }

    fn debug_print_line_string<T: geo_types::CoordFloat>(line_string: &LineString<T>) {
        for c in line_string.0.iter() {
            print!("{};", debug_print_point(c));
        }
        println!(" ");
        println!("{:?}", line_string.0.len());
    }

    fn debug_print_point<T: geo_types::CoordFloat>(c: &Coord<T>) -> String {
        format!(
            "{:.0?} {:.0?}",
            (c.x * T::from(1000).expect("Always good")),
            (c.y * T::from(1000).expect("Always good")),
        )
    }

    fn debug_print_bezier_string<T: geo_types::CoordFloat>(bezier_string: &BezierString<T>) {
        for (i, segment) in bezier_string.0.iter().enumerate() {
            match segment {
                BezierSegment::Bezier(bc) => {
                    if i == 0 {
                        print!("{}", debug_print_point(&bc.start));
                    }
                    print!(
                        " 1;{};{};{}",
                        debug_print_point(&bc.handle1),
                        debug_print_point(&bc.handle2),
                        debug_print_point(&bc.end)
                    )
                }
                BezierSegment::Line(line) => {
                    if i == 0 {
                        print!("{}", debug_print_point(&line.start));
                    }
                    print!(";{}", debug_print_point(&line.end));
                }
            }
        }
        println!(";");
        println!("{:?}", bezier_string.num_points())
    }

    #[test]
    fn single_bezier_curve_to_line_string_f64() {
        let segment = bezier_string_f64().0[0].clone();
        let curve = match segment {
            BezierSegment::Bezier(bezier_curve) => bezier_curve,
            BezierSegment::Line(_) => panic!("We want to test the bezier not line"),
        };
        let line_string = curve
            .clone()
            .to_line_string(0.1)
            .expect("Approximation failed");

        let max_dist = actual_max_dist_between_bezier_and_line_string(&curve, &line_string, 1);
        println!("{:?}", max_dist);

        assert!(max_dist < 0.1);
    }

    #[test]
    fn bezier_string_to_line_string_f64() {
        let curve = bezier_string_f64();
        let line_string = curve
            .clone()
            .to_line_string(0.1)
            .expect("Approximation failed");

        let mut max_error = 0.;
        for segment in curve.0.iter() {
            match segment {
                BezierSegment::Line(_) => continue,
                BezierSegment::Bezier(bc) => {
                    // find the start and end index in line_string for bc
                    // find the max error
                    // reparameterize
                    // check against 0.1
                    let max_dist =
                        actual_max_dist_between_bezier_and_line_string(bc, &line_string, 1);
                    println!("{:?}", max_dist);
                    if max_dist > max_error {
                        max_error = max_dist
                    }
                }
            }
        }
        debug_print_bezier_string(&curve);
        debug_print_line_string(&line_string);
        assert!(max_error < 0.1);
    }

    #[test]
    fn single_line_to_bezier_segment_f64() {
        let segment = LineString(
            line_string_f64()
                .0
                .into_iter()
                .enumerate()
                .filter(|(i, _)| *i < 2)
                .map(|(_, c)| c)
                .collect(),
        );
        let curve =
            BezierString::from_line_string(segment, 0.1).expect("Could not convert to bezier");
        assert!(curve.num_points() == 2)
    }

    #[test]
    fn line_string_to_bezier_string_f64() {
        let line_string = line_string_f64();
        let curve = BezierString::from_line_string(line_string.clone(), 0.1)
            .expect("Could not convert to bezier");

        let mut max_error = 0.;
        for segment in curve.0.iter() {
            match segment {
                BezierSegment::Line(_) => continue,
                BezierSegment::Bezier(bc) => {
                    let max_dist =
                        actual_max_dist_between_bezier_and_line_string(bc, &line_string, 0);
                    println!("{:?}", max_dist);
                    if max_dist > max_error {
                        max_error = max_dist
                    }
                }
            }
        }
        debug_print_bezier_string(&curve);
        debug_print_line_string(&line_string);
        assert!(max_error < 0.1);
    }
}
