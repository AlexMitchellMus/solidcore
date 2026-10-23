//! Small dense linear solves, on the stack: every Newton step and velocity
//! here is a system of at most 4 rows and 4 columns. Heap vectors for these
//! cost more than the arithmetic, at a dozen per marching step.

pub(crate) const MAX: usize = 4;

/// J x = r with m rows and n columns, m and n at most 4.
#[derive(Clone, Copy, Debug)]
pub(crate) struct System {
    pub m: usize,
    pub n: usize,
    pub j: [[f64; MAX]; MAX],
    pub r: [f64; MAX],
}

impl System {
    /// m rows of n columns, every entry zero.
    pub fn new(m: usize, n: usize) -> System {
        System {
            m,
            n,
            j: [[0.0; MAX]; MAX],
            r: [0.0; MAX],
        }
    }

    /// The 3 × 2 system [a b] x = r: a vector's coordinates along two directions.
    pub fn columns2(a: [f64; 3], b: [f64; 3], r: [f64; 3]) -> System {
        let mut s = System::new(3, 2);
        for k in 0..3 {
            s.j[k] = [a[k], b[k], 0.0, 0.0];
            s.r[k] = r[k];
        }
        s
    }
}

/// Solve the square n × n system A x = b by Gaussian elimination with partial
/// pivoting; None when singular or not finite.
fn solve_square(n: usize, mut a: [[f64; MAX]; MAX], mut b: [f64; MAX]) -> Option<[f64; MAX]> {
    for col in 0..n {
        let mut pivot = col;
        for r in col + 1..n {
            if a[r][col].abs() > a[pivot][col].abs() {
                pivot = r;
            }
        }
        if !(a[pivot][col].abs() > 1e-300) {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for r in col + 1..n {
            let f = a[r][col] / a[col][col];
            for k in col..n {
                a[r][k] -= f * a[col][k];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = [0.0; MAX];
    for r in (0..n).rev() {
        let mut sum = b[r];
        for k in r + 1..n {
            sum -= a[r][k] * x[k];
        }
        x[r] = sum / a[r][r];
    }
    x[..n].iter().all(|v| v.is_finite()).then_some(x)
}

/// (JᵀJ + damping I) x = Jᵀr.
pub(crate) fn solve_regularized(s: &System, damping: f64) -> Option<[f64; MAX]> {
    let (m, n) = (s.m, s.n);
    let mut a = [[0.0; MAX]; MAX];
    let mut b = [0.0; MAX];
    for p in 0..n {
        for q in 0..n {
            a[p][q] = (0..m).map(|k| s.j[k][p] * s.j[k][q]).sum();
        }
        a[p][p] += damping;
        b[p] = (0..m).map(|k| s.j[k][p] * s.r[k]).sum();
    }
    solve_square(n, a, b)
}

/// The damped least-squares step for J x = r, damped by |r|^2 (Levenberg and
/// Marquardt with Yamashita and Fukushima's choice of damping): the Newton step
/// where J has full rank, and a step onto the nearest solution where it does
/// not. An edge lying along the other surface makes the curve/surface Jacobian
/// rank 2, and a plain solve then throws its root away. The damping vanishes
/// with the residual, so a grazing double root still converges at Newton's
/// rate; a fixed damping of 1e-12 of the normal matrix stalled there 2e-7 to
/// 1e-6 short of the touch (ssi-corpus c5_hemisphere_equator, where a plane
/// piece's edge grazes the circle), and each stalled root traced the arc again.
pub(crate) fn solve_damped(s: &System) -> Option<[f64; MAX]> {
    solve_regularized(s, (0..s.m).map(|k| s.r[k] * s.r[k]).sum())
}

/// The Newton step for J x = r: solved directly where J is square (least
/// squares where it is tall), and damped only where that fails or would throw
/// the parameters out of range (`fits` false). The direct solve comes first
/// because the normal equations square J's conditioning: at a crossing of
/// 2e-6 rad, a loop of radius 1e-6 (ssi-corpus c7_tiny_loop_parab_1e-6), that
/// lost the corrector its convergence. The damped step is for a singular J: an
/// edge lying along the other surface, or the free parameter at a pole.
pub(crate) fn newton_step(s: &System, fits: impl Fn(&[f64]) -> bool) -> Option<[f64; MAX]> {
    let direct = if s.m == s.n {
        solve_square(s.n, s.j, s.r)
    } else {
        solve_regularized(s, 0.0)
    };
    match direct {
        Some(step) if fits(&step[..s.n]) => Some(step),
        _ => solve_damped(s),
    }
}
