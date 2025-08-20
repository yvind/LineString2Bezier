//! Approximate a geo_types::LineString<T> with a cubic bezier spline
//!
//! Pass the line_string and a maximum allowed error to BezierString::from_line_string
//! The error is the maximum distance between the line_string and the bezier spline
//! evaluated at the line_string vertices and at the middle point between the vertices

use geo_types::{Coord, CoordFloat, LineString};

/// Simple struct to represent a bezier segment
/// If the `Option` is `None` then the segment is a simple straight line segment
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BezierSegment<T: CoordFloat = f64> {
    pub start: Coord<T>,
    pub handles: Option<(Coord<T>, Coord<T>)>,
    pub end: Coord<T>,
}

impl<T: CoordFloat> BezierSegment<T> {
    /// Is this segment a cubic bezier?
    pub fn is_bezier_segment(&self) -> bool {
        self.handles.is_some()
    }

    pub fn new(
        start: Coord<T>,
        handles: Option<(Coord<T>, Coord<T>)>,
        end: Coord<T>,
    ) -> BezierSegment<T> {
        BezierSegment {
            start,
            handles,
            end,
        }
    }

    /// Get an approximation of the length of the cubic bezier segment
    /// Taken as the mean of the distance between the end points
    /// and the distance from start to end through the handles
    pub fn length_approximation(&self) -> T {
        if let Some((handle_1, handle_2)) = self.handles {
            ((self.end - handle_2).magnitude()
                + (handle_2 - handle_1).magnitude()
                + (handle_1 - self.start).magnitude()
                + (self.end - self.start).magnitude())
                / (T::one() + T::one())
        } else {
            (self.end - self.start).magnitude()
        }
    }

    /// Get a LineString with at most error deviation from the Bezier segment
    /// The deviation is measured at the mid point between neighboring nodes
    pub fn to_line_string(self, error: T) -> LineString<T> {
        if let Some((handle_1, handle_2)) = self.handles {
            let two = T::one() + T::one();

            let mut num_sample_points = 2;

            let mut coords = Vec::with_capacity(num_sample_points + 1);
            let mut ts = Vec::with_capacity(num_sample_points);

            let bez_curve_slice = &self.as_slice().unwrap();

            let mut current_error = two * error;
            while current_error > error {
                ts.clear();
                coords.clear();

                let sample_interval = T::one() / T::from(num_sample_points).unwrap();

                coords.push(self.start);
                ts.push(T::zero());
                for i in 1..num_sample_points {
                    let t = T::from(i).unwrap() * sample_interval;

                    ts.push(t);

                    let coord = self.start * bezier_basis_0(t)
                        + handle_1 * bezier_basis_1(t)
                        + handle_2 * bezier_basis_2(t)
                        + self.end * bezier_basis_3(t);

                    coords.push(coord);
                }
                ts.push(T::one());
                coords.push(self.end);

                (current_error, _) =
                    compute_max_error(&coords, 0, coords.len() - 1, bez_curve_slice, &ts);

                num_sample_points += 2;
            }

            LineString(coords)
        } else {
            LineString(vec![self.start, self.end])
        }
    }

    fn as_slice(&self) -> Option<[Coord<T>; 4]> {
        if let Some((handle_1, handle_2)) = self.handles {
            Some([self.start, handle_1, handle_2, self.end])
        } else {
            None
        }
    }
}

impl<T> From<[Coord<T>; 4]> for BezierSegment<T>
where
    T: CoordFloat,
{
    fn from(value: [Coord<T>; 4]) -> Self {
        BezierSegment {
            start: value[0],
            handles: Some((value[1], value[2])),
            end: value[3],
        }
    }
}

// Bezier basis functions
#[inline]
fn bezier_basis_0<T: CoordFloat>(t: T) -> T {
    (T::one() - t).powi(3)
}

#[inline]
fn bezier_basis_1<T: CoordFloat>(t: T) -> T {
    T::from(3.0).unwrap() * t * (T::one() - t).powi(2)
}

#[inline]
fn bezier_basis_2<T: CoordFloat>(t: T) -> T {
    T::from(3.0).unwrap() * t.powi(2) * (T::one() - t)
}

#[inline]
fn bezier_basis_3<T: CoordFloat>(t: T) -> T {
    t.powi(3)
}

/// A BezierString is simply a vector of [BezierSegment]s
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BezierString<T: CoordFloat = f64>(pub Vec<BezierSegment<T>>);

impl<T: CoordFloat> BezierString<T> {
    /// The number of vertices and handles in this [BezierString]
    /// The endpoints of each [BezierSegment] are not counted double
    pub fn num_points(&self) -> usize {
        let mut num_points = 0;

        for segment in self.0.iter() {
            if segment.is_bezier_segment() {
                num_points += 3;
            } else {
                num_points += 1;
            }
        }
        num_points + 1
    }

    /// Convert a [BezierString] to a [geo_types::LineString] with a maximum of `error` deviation between the two
    pub fn to_line_string(self, error: T) -> LineString<T> {
        let mut line = LineString::new(Vec::new());

        let mut first_segment = true;
        for segment in self.0 {
            let iter = segment.to_line_string(error).0.into_iter();
            if first_segment {
                line.0.extend(iter);
                first_segment = false;
            } else {
                line.0.extend(iter.skip(1));
            }
        }

        line
    }

    /// Convert a [geo_types::LineString] to a [BezierString] with a maximum of `error` deviation between the two
    /// evaluated at the line_string vertices and at the middle point between the vertices
    pub fn from_line_string(line_string: LineString<T>, error: T) -> BezierString<T> {
        let n_pts = line_string.0.len();
        if n_pts <= 2 || (n_pts <= 3 && line_string.is_closed()) {
            let mut segments = vec![];

            let mut prev_coord = line_string.0[0];
            for p in line_string.0.iter().skip(1) {
                segments.push(BezierSegment::new(prev_coord, None, *p));
                prev_coord = *p;
            }

            return BezierString(segments);
        }

        let mut bezier_segments = vec![];

        let mut tangent_right = compute_right_tangent(line_string.0.as_slice(), 0);
        let mut tangent_left = compute_left_tangent(line_string.0.as_slice(), n_pts - 1);

        if line_string.is_closed() {
            let diff = tangent_right - tangent_left;
            if diff != Coord::zero() {
                tangent_right = diff.try_normalize().unwrap();
                tangent_left = -tangent_right;
            }
        }

        // recursively tries to fit bezier segments to the line_string
        fit_cubic(
            &line_string.0,
            0,
            n_pts - 1,
            tangent_right,
            tangent_left,
            error,
            &mut bezier_segments,
        );
        BezierString(bezier_segments)
    }
}
// helper functions

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
    let (mut max_error, mut split_point) =
        compute_max_error(polyline, first, last, &bez_curve, &ts);

    if max_error < error {
        bezier_string.push(bez_curve.into());
        return;
    }

    // If error not too large, try one step of Newton-Rhapson
    if max_error < T::from(2.).unwrap() * error {
        ts = reparameterize(polyline, first, last, &ts, &bez_curve);
        bez_curve = generate_bezier(polyline, first, last, &ts, &tangent_start, &tangent_end);
        (max_error, split_point) = compute_max_error(polyline, first, last, &bez_curve, &ts);

        if max_error < error {
            bezier_string.push(bez_curve.into());
            return;
        }
    }

    // Fitting failed, split at the point of max error and fit each part recursively
    let tangent_center = compute_center_tangent(polyline, split_point);
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
        c[0] = c[0] + a[i][0].dot_product(a[i][0]);
        c[1] = c[1] + a[i][0].dot_product(a[i][1]);
        c[2] = c[2] + a[i][1].dot_product(a[i][1]);

        let tmp = polyline[first + i]
            - (polyline[first] * (bezier_basis_0(t) + bezier_basis_1(t))
                + polyline[last] * (bezier_basis_2(t) + bezier_basis_3(t)));

        x[0] = x[0] + a[i][0].dot_product(tmp);
        x[1] = x[1] + a[i][1].dot_product(tmp);
    }

    // Compute the determinants
    let det_c0_c1 = c[0] * c[2] - c[1] * c[1];
    let det_c0_x = c[0] * x[1] - c[1] * x[0];
    let det_x_c1 = x[0] * c[2] - x[1] * c[1];

    // Derive alpha values
    let alpha_l = if det_c0_c1 == T::zero() {
        T::zero()
    } else {
        det_x_c1 / det_c0_c1
    };
    let alpha_r = if det_c0_c1 == T::zero() {
        T::zero()
    } else {
        det_c0_x / det_c0_c1
    };

    // If alpha negative, use the Wu/Barsky heuristic
    let seg_length = (polyline[last] - polyline[first]).magnitude();
    let epsilon = T::from(1.0e-6).unwrap() * seg_length;

    if alpha_l < epsilon || alpha_r < epsilon {
        let dist = seg_length / T::from(3.0).unwrap();
        return [
            polyline[first],
            polyline[first] + *start_tangent * dist,
            polyline[last] + *end_tangent * dist,
            polyline[last],
        ];
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
fn compute_right_tangent<T: CoordFloat>(polyline: &[Coord<T>], end: usize) -> Coord<T> {
    (polyline[end + 1] - polyline[end]).try_normalize().unwrap()
}

#[inline]
fn compute_left_tangent<T: CoordFloat>(polyline: &[Coord<T>], end: usize) -> Coord<T> {
    (polyline[end - 1] - polyline[end]).try_normalize().unwrap()
}

#[inline]
fn compute_center_tangent<T: CoordFloat>(polyline: &[Coord<T>], center: usize) -> Coord<T> {
    (polyline[center + 1] - polyline[center - 1])
        .try_normalize()
        .unwrap()
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
    for i in (first + 1)..=last {
        ts.push(ts[i - first - 1] + (polyline[i] - polyline[i - 1]).magnitude());
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
    polyline: &[Coord<T>],
    first: usize,
    last: usize,
    ts: &[T],
    bez_curve: &[Coord<T>],
) -> Vec<T> {
    let mut new_ts = vec![T::zero(); last - first + 1];

    for i in first..=last {
        new_ts[i - first] = newton_raphson(bez_curve, polyline[i], ts[i - first]);
    }
    new_ts
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
    if denominator == T::zero() {
        return t;
    }
    // Compute f(t)
    let numerator = (bez_t.x - p.x) * (qp_t.x) + (bez_t.y - p.y) * (qp_t.y);

    // t = t - f(t)/f'(t)
    t - (numerator / denominator)
}

// max of distance between polyline points and the fitted curve
// checks even the halfway point of every line segment
fn compute_max_error<T: CoordFloat>(
    polyline: &[Coord<T>],
    first: usize,
    last: usize,
    bez_curve: &[Coord<T>],
    ts: &[T],
) -> (T, usize) {
    let two = T::one() + T::one();

    let mut split_point = (last + first).div_ceil(2);
    let mut max_dist = T::zero();
    for i in first..last {
        let p = evaluate_bezier(3, bez_curve, ts[i - first]);
        let dist = (p - polyline[i]).magnitude_squared();

        if dist >= max_dist {
            max_dist = dist;
            split_point = i;
        }

        let p = evaluate_bezier(3, bez_curve, (ts[i - first] + ts[i - first + 1]) / two);
        let dist = (p - (polyline[i] + polyline[i + 1]) / two).magnitude_squared();
        if dist >= max_dist {
            max_dist = dist;
            split_point = i.max(first + 1);
        }
    }
    (max_dist, split_point)
}

// evaluate the bezier value at time t
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
    fn dot_product(&self, other: Rhs) -> T;
    fn try_normalize(&self) -> Option<Self>;
    fn is_finite(&self) -> bool;
}

impl<T: CoordFloat> VectorTraits<T> for Coord<T> {
    fn magnitude_squared(&self) -> T {
        self.x * self.x + self.y * self.y
    }

    fn magnitude(&self) -> T {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    fn dot_product(&self, other: Self) -> T {
        self.x * other.x + self.y * other.y
    }

    fn try_normalize(&self) -> Option<Self> {
        let magnitude = self.magnitude();
        let result = *self / magnitude;

        if result.is_finite() && magnitude.is_finite() {
            Some(result)
        } else {
            None
        }
    }

    fn is_finite(&self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }
}
