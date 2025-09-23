-- | Operations for performing LU-decomposition with partial (row) pivoting and
-- operations for solving systems of linear equations using the LU-decomposition
-- functionality. The module `mk_lup` is parameterised over an ordered
-- field. Examples of ordered fields include `f64` and `f32`.

import "linalg"

local module type lup = {

  -- | Type of elements
  type t

  -- | Type of permutations
  type perm[m]

  -- | Perform LU-decomposition with partial (row) pivoting. A call `lup A`
  -- returns a pair `(LU,P)` of a matrix `LU` and a permutation `P`, such that
  -- `permute P A = matmul (lower LU) (upper LU)`. Notice that both the lower
  -- and upper triangular matrices are embedded in the resulting matrix `LU`.
  val lup [m] : (mat:*[m][m]t) -> ([m][m]t, perm[m])

  -- | Perform a permutation of a vector given a permutation.
  val permute 'a [m] : perm[m] -> *[m]a -> *[m]a

  -- | Perform an inverse permutation of a vector given a permutation. We have
  -- `permute_inv p (permute p v) = v` for any permutation `p` of length `m` and
  -- vector `v` of length `m`.
  val permute_inv 'a [m] : perm[m] -> *[m]a -> *[m]a

  -- | Forward solving. Reads only lower part of argument matrix (excluding
  -- diagonal)
  val forsolve [n] : [n][n]t -> [n]t -> [n]t

  -- | Backward solving. Reads only upper part of argument matrix (including
  -- diagonal)
  val backsolve [n] : [n][n]t -> [n]t -> [n]t

  -- | Solve a linear system using LU-decomposition with partial (row) pivoting.
  val ols [m] : *[m][m]t -> *[m]t -> *[m]t

}

-- | LU-decomposition module parameterised over an ordered field.
module mk_lup (T: ordered_field) : lup with t = T.t
= {

  type t = T.t

  def z = T.i64 0

  def dotprod [n] (a: [n]T.t) (b: [n]T.t): T.t =
    map2 (T.*) a b |> reduce (T.+) (T.i64 0)

  def swap [m] 'a (p:*[m]a) (i:i64) (j:i64) : *[m]a  =
    let tmp = copy (p[i])
    let p[i] = copy (p[j])
    in p with [j] = tmp

  type perm[m] = *[m]i64

  def permute 'a [m] (p:perm[m]) (v:*[m]a) : *[m]a =
    map (\i -> v[i]) p

  def permute_inv 'a [m] (p:perm[m]) (v:*[m]a) : *[m]a =
    scatter (copy v) p v

  def step [m] (mat:*[m][m]t) (p:*[m]i64) (j:i64) : (*[m][m]t, *[m]i64) =
    let jp = reduce (\(a:(i64,t)) (b:(i64,t)) : (i64,t) -> if a.1 T.> b.1 then a else b) (-1,z)
		    (map2 (\i v -> (i+j,T.abs v)) (iota (m-j)) (mat[j:,j])) |> (.0)
    let pv = copy(mat[jp][j])
    let mat = swap mat jp j
    let p = swap p jp j
    let mat[j+1:,j] = map (T./ pv) mat[j+1:,j]
    let mat[j+1:,j+1:] =
      tabulate_2d (m-j-1) (m-j-1)
		  (\r c ->
		     mat[r+j+1][c+j+1] T.- (mat[r+j+1][j] T.* mat[j][c+j+1])
		  )
    in (mat,p)

  -- Perform LU-decomposition with partial (row) pivoting. A call `lup a`
  -- returns a pair `(LU,p)` of a matrix `LU` and a permutation `p`, such that
  -- `permute p a = LU`.
  def lup [m] (mat:*[m][m]t) : ([m][m]t, [m]i64) =
    loop (mat, p) = (mat, iota m) for j < m do step mat p j

  -- reads only lower part of L and assumes implicit diagonal elements are 1
  def forsolve [n] (L:[n][n]t) (b:[n]t) : [n]t =
    let y : *[n]t = replicate n (T.i64 0)
    in loop y for i in 0..<n do
       let sum = dotprod L[i,:i] y[:i]
       let y[i] = copy(b[i] T.- sum)
       in y

  -- reads only upper par of U including diagonal
  def backsolve [n] (U:[n][n]t) (y:[n]t) : [n]t =
    let x : *[n]t = replicate n (T.i64 0)
    in loop x for j in 0..<n do
       let i = n - j - 1
       let sum = dotprod U[i,i+1:n] x[i+1:n]
       let x[i] = copy(y[i] T.- sum) T./ U[i,i]
       in x

  def ols [n] (a: *[n][n]t) (b:*[n]t) : *[n]t =
    let (LU,p) = lup a
    in backsolve LU (forsolve LU (permute p b))

}
