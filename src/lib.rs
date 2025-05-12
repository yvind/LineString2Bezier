//! Approximate a geo_types::LineString<T> with a cubic bezier spline
//!
//! Pass the linestring and a maximum allowed error to BezierString::from_linestring
//! The error is the maximum distance between the linestring and the bezier spline
//! evaulated at the linestring vertices and at the middle point between the vertices

use geo_types::{Coord, CoordFloat, LineString};

/// Simple struct to represent a bezier segement
/// If the `Option` is `None` then the segment is a simple straight line segment
#[derive(Clone, Debug)]
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

/// A BezierString is simply a vector of [BezierSegment]s
#[derive(Debug)]
pub struct BezierString<T: CoordFloat = f64>(pub Vec<BezierSegment<T>>);

impl<T: CoordFloat> BezierString<T> {
    /// The number of vertices and handles in this [BezierString]
    /// The endpoints of each [BezierSegment] are not counted double
    pub fn num_points(&self) -> usize {
        let mut num_points = 0;

        for segment in self.0.iter() {
            num_points += 1;
            if segment.is_bezier_segment() {
                num_points += 2
            }
        }
        num_points + 1
    }

    /// Convert a [geo_types::LineString] to a [BezierString] with a maximum of `error` deviation between the two
    /// evaulated at the linestring vertices and at the middle point between the vertices
    pub fn from_linestring(linestring: LineString<T>, error: T) -> BezierString<T> {
        let n_pts = linestring.0.len();
        if n_pts <= 2 || (n_pts <= 3 && linestring.is_closed()) {
            let mut segments = vec![];

            let mut prev_coord = linestring.0[0];
            for p in linestring.0.iter().skip(1) {
                segments.push(BezierSegment::new(prev_coord, None, *p));
                prev_coord = *p;
            }

            return BezierString(segments);
        }

        let mut bezier_segments = vec![];

        let mut tangent_right = Self::compute_right_tangent(linestring.0.as_slice(), 0);
        let mut tangent_left = Self::compute_left_tangent(linestring.0.as_slice(), n_pts - 1);

        if linestring.is_closed() {
            tangent_right = (tangent_right - tangent_left).try_normalize().unwrap();
            tangent_left = -tangent_right;
        }

        // recursivly tries to fit bezier segments to the linestring
        Self::fit_cubic(
            &linestring.0,
            0,
            n_pts - 1,
            tangent_right,
            tangent_left,
            error,
            &mut bezier_segments,
        );
        BezierString(bezier_segments)
    }

    // recursive function to fit bezier curve segments to linestring
    fn fit_cubic(
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
        let mut ts = Self::chord_length_parameterize(polyline, first, last);
        let mut bez_curve = Self::generate_bezier(
            polyline,
            first,
            last,
            ts.as_slice(),
            &tangent_start,
            &tangent_end,
        );

        // Find max deviation
        let (mut max_error, mut split_point) =
            Self::compute_max_error(polyline, first, last, &bez_curve, &ts);

        if max_error < error {
            bezier_string.push(bez_curve.into());
            return;
        }

        // If error not too large, try one step of newton-rhapson
        if max_error < T::from(2.).unwrap() * error {
            ts = Self::reparameterize(polyline, first, last, &ts, &bez_curve);
            bez_curve =
                Self::generate_bezier(polyline, first, last, &ts, &tangent_start, &tangent_end);
            (max_error, split_point) =
                Self::compute_max_error(polyline, first, last, &bez_curve, &ts);

            if max_error < error {
                bezier_string.push(bez_curve.into());
                return;
            }
        }

        // Fitting failed, split at the point of max error and fit each part recursively
        let tangent_center = Self::compute_center_tangent(polyline, split_point);
        let tangent_center_neg = -tangent_center;
        Self::fit_cubic(
            polyline,
            first,
            split_point,
            tangent_start,
            tangent_center_neg,
            error,
            bezier_string,
        );
        Self::fit_cubic(
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
    fn generate_bezier(
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
            a.push([*start_tangent * Self::b1(t), *end_tangent * Self::b2(t)]);
        }

        // Create the C and X matrices
        // C is symmetric 2x2, sum of indecies in 2x2 gives index in flat array
        // X is a 2x1 vector
        let mut c = [T::zero(); 3];
        let mut x = [T::zero(); 2];

        for (i, &t) in ts.iter().enumerate() {
            c[0] = c[0] + a[i][0].dot_product(a[i][0]);
            c[1] = c[1] + a[i][0].dot_product(a[i][1]);
            c[2] = c[2] + a[i][1].dot_product(a[i][1]);

            let tmp = polyline[first + i]
                - (polyline[first] * (Self::b0(t) + Self::b1(t))
                    + polyline[last] * (Self::b2(t) + Self::b3(t)));

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

    // Bezier basis functions
    #[inline]
    fn b0(t: T) -> T {
        (T::one() - t).powi(3)
    }

    #[inline]
    fn b1(t: T) -> T {
        T::from(3.0).unwrap() * t * (T::one() - t).powi(2)
    }

    #[inline]
    fn b2(t: T) -> T {
        T::from(3.0).unwrap() * t.powi(2) * (T::one() - t)
    }

    #[inline]
    fn b3(t: T) -> T {
        t.powi(3)
    }

    // Vertex tangent functions
    #[inline]
    fn compute_right_tangent(polyline: &[Coord<T>], end: usize) -> Coord<T> {
        (polyline[end + 1] - polyline[end]).try_normalize().unwrap()
    }

    #[inline]
    fn compute_left_tangent(polyline: &[Coord<T>], end: usize) -> Coord<T> {
        (polyline[end - 1] - polyline[end]).try_normalize().unwrap()
    }

    #[inline]
    fn compute_center_tangent(polyline: &[Coord<T>], center: usize) -> Coord<T> {
        (polyline[center + 1] - polyline[center - 1])
            .try_normalize()
            .unwrap()
    }

    // normalized length along linestring from start of segment to every vertex in segment
    // the t values to be used for computing error of fitted bezier
    fn chord_length_parameterize(polyline: &[Coord<T>], first: usize, last: usize) -> Vec<T> {
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
    // use newton rhapson to try and refine the parameterization
    fn reparameterize(
        polyline: &[Coord<T>],
        first: usize,
        last: usize,
        ts: &[T],
        bez_curve: &[Coord<T>],
    ) -> Vec<T> {
        let mut new_ts = vec![T::zero(); last - first + 1];

        for i in first..=last {
            new_ts[i - first] = Self::newton_raphson(bez_curve, polyline[i], ts[i - first]);
        }
        new_ts
    }

    // bez_curve Q(u) at time t is supposed to be p, refine t
    fn newton_raphson(bez_curve: &[Coord<T>], p: Coord<T>, t: T) -> T {
        // Q(t)
        let bez_t = Self::evaluate_bezier(3, bez_curve, t);

        // Cubic bez prime is quadratic
        let mut bez_prime = [Coord::<T>::zero(); 3];
        // Cubic bez double prime is linear
        let mut bez_double_prime = [Coord::<T>::zero(); 2];

        // Generate control vertices for Q'
        for i in 0..3 {
            bez_prime[i].x = (bez_curve[i + 1].x - bez_curve[i].x) * T::from(3.0).unwrap();
            bez_prime[i].y = (bez_curve[i + 1].y - bez_curve[i].y) * T::from(3.0).unwrap();
        }

        // Generate control vertices for Q''
        for i in 0..2 {
            bez_double_prime[i].x = (bez_prime[i + 1].x - bez_prime[i].x) * T::from(2.0).unwrap();
            bez_double_prime[i].y = (bez_prime[i + 1].y - bez_prime[i].y) * T::from(2.0).unwrap();
        }

        // Compute Q'(t) and Q''(t)
        let qp_t = Self::evaluate_bezier(2, &bez_prime, t);
        let qpp_t = Self::evaluate_bezier(1, &bez_double_prime, t);

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
    fn compute_max_error(
        polyline: &[Coord<T>],
        first: usize,
        last: usize,
        bez_curve: &[Coord<T>],
        ts: &[T],
    ) -> (T, usize) {
        let mut split_point = (last + first).div_ceil(2);
        let mut max_dist = T::zero();
        for i in first..last {
            let p = Self::evaluate_bezier(3, bez_curve, ts[i - first]);
            let dist = (p - polyline[i]).magnitude_squared();

            if dist >= max_dist {
                max_dist = dist;
                split_point = i;
            }

            let p = Self::evaluate_bezier(
                3,
                bez_curve,
                (ts[i - first] + ts[i - first + 1]) / T::from(2.).unwrap(),
            );
            let dist =
                (p - (polyline[i] + polyline[i + 1]) / T::from(2.).unwrap()).magnitude_squared();
            if dist >= max_dist {
                max_dist = dist;
                split_point = i.max(first + 1); //if i == first { i + 1 } else { i };
            }
        }
        (max_dist, split_point)
    }

    // evaluate the bezier value at time t
    fn evaluate_bezier(degree: usize, bezier_segment: &[Coord<T>], t: T) -> Coord<T> {
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
