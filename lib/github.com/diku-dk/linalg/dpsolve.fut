-- | Generic solver for fix-point (Bellman) equations. This module implements a
-- generic solver for finding solutions to fix-point (Bellman) equations. The
-- module may be used for finding solutions to dynamic programming problems. The
-- implementation is based on John Rust's poly-algorithm, which combines the
-- strategies of successive applications and the Newton-Kantorovich method.
--
-- The module is enriched with a version of the solver that uses Futhark's
-- support for automatic differentiation (forward mode AD) to compute the
-- Jacobian for the passed fix-point function.
--
-- For more details, see Rust's description in Section 4 of
-- https://editorialexpress.com/jrust/sdp/ndp.pdf and a demonstration of its use
-- in https://elsman.com/pdf/dpsolve-2025-09-27.pdf.

import "lu"

local
-- | Module type specifying generic solvers for Bellman equations. Notice that
-- this module type is declared `local`, which means that it cannot directly be
-- referenced by name from the outside.  By not exposing the module type, the
-- module may be extended with new members in minor versions.
module type dpsolve = {
  -- | The scalar type.
  type t

  -- | The type of solver parameters.
  type param = {
    sa_max          : i64,  -- Maximum and minimum numbers of
    sa_min          : i64,  --   successive approximation steps.
    sa_tol          : t,    -- Stopping tolerance for successive
                            --   approximation.
    max_fxpiter     : i64,  -- Maximum number of times to switch
                            --   between Newton-Kantorovich and
                            --   successive approximation iterations.
    pi_max          : i64,  -- Maximum number of Newton-Kantorovich steps.
    pi_tol          : t,    -- Final exit tolerance in fixed point
                            --   algorithm, measured in units of
                            --   numerical precision
    tol_ratio       : t     -- Relative tolerance before switching
                            --   to N-K algorithm when discount factor is
                            --   supplied as input in `poly`.
    }

  -- | The default parameter value
  val default : param

  -- | Find a fix-point for the function `f` using successive approximation with
  -- the initial guess `v`, parameter `p`, and beta-value `b`. The tolerance,
  -- the minimal, and the maximal number of iterations can be adjusted by
  -- altering the parameter `p`. The argument `b` is the beta discount factor,
  -- which allows for stopping early (when the relative tolerance is close to
  -- `b`). The function returns a quintuple containing an approximate fix-point,
  -- a boolean specifying whether the algorithm converged (according to the
  -- values in `p`), the number of iterations used, and finally, the tolerance
  -- and the relative tolerance of the last two fix-point approximations and the
  -- last two tolerances, respectively (maximum of each dimension).
  val sa [m] : (f:[m]t->[m]t) -> (v:[m]t) -> (p:param) -> (b:t)
               -> ([m]t, bool, i64, t, t)

  -- | Find a fix-point for the function `f` using Newton-Kantorovich iterations
  -- with the initial guess `v`, and parameter `p`. The tolerance, the minimal,
  -- and the maximal number of iterations can be adjusted by altering the
  -- parameter `p`. The function `f` should return a pair of a new next
  -- approximation and the Jacobian matrix for the function `f` relative to the
  -- argument given. The function returns a quintuple containing an approximate
  -- fix-point, a Jacobian matrix for the fix-point, a boolean specifying whether
  -- the algorithm converged (according to the values in `p`), the number of
  -- iterations used, and finally, the tolerance of the last two fix-point
  -- approximations (maximum of each dimension).
  val nk [m] : (f: [m]t -> ([m]t,[m][m]t)) -> (v:[m]t) -> (p:param)
               ->  ([m]t, [m][m]t, bool, i64, t)

  -- | Find a fix-point for the function `f` using a combination of successive
  -- approximation iterations and Newton-Kantorovich iterations. The initial
  -- guess is `v` and the parameter `p` is passed to the calls to `sa` and
  -- `nk`. The function `f` should return a pair of a new next approximation and
  -- the Jacobian matrix for the function `f` relative to the argument
  -- given. The function returns a 7-tuple containing an approximate fix-point, a
  -- Jacobian matrix for the fix-point, a boolean specifying whether the
  -- algorithm converged (according to the values in `p`), the number of
  -- iterations used for the total sa iterations, the total nk iterations, and
  -- the number of round-trips. The 7'th element of the result tuple is the
  -- tolerance of the last two fix-point approximations (maximum of each
  -- dimension).
  val poly [m] : (f: [m]t -> ([m]t,[m][m]t)) -> (v:[m]t) -> (p:param)
                 -> (b:t) -> ([m]t,[m][m]t,bool,i64,i64,i64,t)

  -- | Find a fix-point for the function `f` using a combination of successive
  -- approximation iterations and Newton-Kantorovich iterations. The initial
  -- guess is `v` and the parameter `p` is passed to the calls to `sa` and
  -- `nk`. The function uses forward-mode automatic differentiation to compute
  -- the Jacobian matrix relative to the argument given to `f`. The function
  -- returns a 6-tuple containing an approximate fix-point, a boolean specifying
  -- whether the algorithm converged (according to the values in `p`), the
  -- number of iterations used for the total sa iterations, the total nk
  -- iterations, and the number of round-trips. The 6'th element of the result
  -- tuple is the tolerance of the last two fix-point approximations (maximum of
  -- each dimension).
  val polyad [m] : (f:[m]t->[m]t) -> (v:[m]t) -> (p:param) -> (b:t)
                   -> ([m]t, bool, i64, i64, i64, t)
}

-- | Parameterised module for creating generic solvers.

module mk_dpsolve (r:real) : dpsolve with t = r.t = {

  type t = r.t
  local module lu = mk_lu r

  def blksz : i64 = 16  -- 1 or 16

  def ols [n] (m:[n][n]t) (b:[n]t) : [n]t =
    lu.ols blksz m b

  def eye (n:i64) (m:i64) : [n][m]t =
    tabulate_2d n m (\i j -> r.bool (i == j))

  def zero (n:i64) (m:i64) : [n][m]t =
    replicate n (replicate m (r.i64 0))

  type param = {
    sa_max          : i64,  -- Maximum number of contraction steps
    sa_min          : i64,  -- Minimum number of contraction steps
    sa_tol          : t,    -- Absolute tolerance before (in dpsolve.poly: tolerance before switching
                            --   to N-K algorithm)
    max_fxpiter     : i64,  -- Maximum number of times to switch between Newton-Kantorovich iterations
                            --   and contraction iterations.
    pi_max          : i64,  -- Maximum number of Newton-Kantorovich steps
    pi_tol          : t,    -- Final exit tolerance in fixed point algorithm, measured in units of
                            --   numerical precision
    tol_ratio       : t     -- Relative tolerance before switching to N-K algorithm
                            --   when discount factor is supplied as input in dpsolve.poly
  }

  def default : param =
    {sa_max      = 20,
     sa_min      = 2,
     sa_tol      = r.f64 1.0e-3,
     max_fxpiter = 35,
     pi_max      = 40,
     pi_tol      = r.f64 1.0e-13,
     tol_ratio   = r.f64 1.0e-03
    }

  def sa [m] (bellman : [m]t -> [m]t)
             (V0:[m]t)
             (ap:param)
             (bet:t) : ([m]t, bool, i64, t, t) =
      loop (V0,converged,i,tol,_rtol) = (V0, false, 0, r.i64 0, r.i64 0)
      while !converged && i < ap.sa_max do
        let V = bellman V0
	let tol' = reduce r.max (r.i64 0) (map2 (\a b -> r.(abs(a-b))) V V0)
	let rtol' = if i == 1 then r.i64 1
	 	    else r.(tol' / tol)
	let converged =
             -- Rule 1
             (i > ap.sa_min && (r.(abs(bet-rtol') < ap.tol_ratio)))
          || -- Rule 2
             --let adj = f64.(maximum V0 |> abs |> log10 |> ceil)
             --let ltol = ap.sa_tol * f64.(10 ** adj)
	     let ltol = ap.sa_tol
             in (i > ap.sa_min && r.(tol' < ltol))
	in (V, converged, i+1, tol',rtol')

  def nk [m] (bellman : [m]t -> ([m]t,[m][m]t))
             (V0:[m]t)
             (ap:param) : ([m]t, [m][m]t, bool, i64, t) =
    loop (V0,_dV0,converged,i,_tol) = (V0, zero m m, false, 0, r.i64 1)
      while !converged && i < ap.pi_max do
        let (V1, dV) = bellman V0
        --let V = map2 (-) V0 (la.matvecmul_row (la.inv dV) V1)
	let F = map2 (map2 (r.-)) (eye m m) dV
  	--let _hermitian =
	--  map2 (map2 (r.==)) F (transpose F) |> map (reduce (&&) true) |> reduce (&&) true
	let V = map2 (r.-) V0 (ols F (map2 (r.-) V0 V1))  -- NK-iteration
	-- do additional SA iteration for stability and accurate measure of error bound
	let (V0, _) = bellman V
	let tol' = r.maximum (map2 (\a b -> r.(abs(a-b))) V V0) -- tolerance

	-- adjusting the N-K tolerance to the magnitude of ev
	let adj = r.(maximum V0 |> abs |> log10 |> ceil)
	let ltol = r.(ap.pi_tol * (i64 10 ** adj)) -- Adjust final tolerance
        -- ltol=ap.pi_tol  -- tolerance

        let converged = r.(tol' < ltol) -- Convergence achieved
        in (V0, dV, converged, i+1, tol')

  -- dpsolve.poly(f,v0,ap,bet): Solve for fixed point using a combination of
  -- Successive Approximations (SA) and Newton-Kantorovich (NK) iterations.  The
  -- argument `ap` holds algorithm parameters. The argument `bet` is a discount
  -- factor. Enters rule for stopping SA and switching to NK iterations.  SA
  -- should stop prematurely when relative tolerance is close to bet.  The
  -- argument `v0` is the initial value guess. The function returns a fix-point
  -- `V`, the Frechet derivative of the Bellman operator, and various count
  -- information including the obtained tolerance.

  def poly [m] (bellman : [m]t -> ([m]t,[m][m]t))
               (V0:[m]t)
               (ap:param)
               (bet:t) : ([m]t,[m][m]t,bool,i64,i64,i64,t) =
    loop (V0,_dV,converged,i,j,k,_tol) = (V0, zero m m, false, 0, 0, 0, r.i64 1)
      while !converged && k < ap.max_fxpiter do
        -- poly-algorithm loop (switching between sa and nk)
        let (V1,_,i',_,_) = sa ((.0) <-< bellman) V0 ap bet
	let (V2, dV, c2, j', tol) = nk bellman V1 ap
        in (V2,dV,c2,i+i',j+j',k+1,tol)

  def idd n i = tabulate n (\j -> if i==j then r.f32 1 else r.f32 0)

  def wrapj [n][m] (f: [n]t->[m]t) (x:[n]t) : ([m]t,[m][n]t) =
    (f x, #[sequential_outer] tabulate n (jvp f x <-< idd n) |> transpose)

  def polyad f x ap bet =
    let (y,_,b,i,j,k,tol) = poly (wrapj f) x ap bet
    in (y,b,i,j,k,tol)

}
