-- | [LU-decomposition](https://en.wikipedia.org/wiki/LU_decomposition).

-- | The block size is a tunable parameter.  16 and 32 are decent
-- sizes.  This decomposition returns the L and U parts embedded into
-- a single square matrix.
module mk_lu(T: float) : {
  val lu [m] : (block_size: i64) -> (mat: [m][m]T.t) -> [m][m]T.t
} = {

def dotprod [n] (a: [n]T.t) (b: [n]T.t): T.t =
  map2 (T.*) a b |> reduce (T.+) (T.i64 0)

def lud_diagonal [b] (a: [b][b]T.t): *[b][b]T.t =
  map (\mat ->
         let mat = copy mat in
         loop (mat: *[b][b]T.t) for i < b-1 do
         let col_elem j =
           if j > i
           then (mat[j,i] T.- (dotprod mat[j,:i] mat[:i,i]))
                T./ mat[i,i]
           else mat[j,i]
         let col = tabulate b col_elem
         let mat[:,i] = col
         let row_elem j =
           if j > i
           then (mat[i+1, j] T.- (dotprod mat[:i+1, j] mat[i+1, :i+1]))
           else mat[i+1, j]
         let row = tabulate b row_elem
         let mat[i+1] = row
         in mat)
      (unflatten (opaque 1) b a)
  |> head

def lud_perimeter_upper [m][b] (diag: [b][b]T.t, a0s: [m][b][b]T.t): *[m][b][b]T.t =
    let a1s = map transpose a0s in
    let a2s =
        map (map (\row0 ->
                    loop row = copy row0 for i < b do
                    let sum =
                      (loop sum= T.i64 0 for k < i do sum T.+ diag[i,k] T.* row[k])
                    let row[i] = copy (row[i] T.- sum)
                    in  row))
            a1s
    in map transpose a2s

def lud_perimeter_lower [b][m] (diag: [b][b]T.t, mat: [m][b][b]T.t): *[m][b][b]T.t =
  map (map (\row0 ->
              loop row = copy row0 for j < b do
              let sum = loop sum=T.i64 0 for k < j do
                          sum T.+ diag[k,j] T.* row[k]
              let row[j] = copy (row[j] T.- sum) T./ diag[j,j]
              in  row))
      mat

def lud_internal [m][b] (top_per: [m][b][b]T.t,
                         lft_per: [m][b][b]T.t,
                         mat_slice: [m][m][b][b]T.t): *[m][m][b][b]T.t =
  let top_slice = map transpose top_per in
  map2 (\mat_arr lft ->
        map2 (\mat_blk top ->
               map2 (\mat_row lft_row ->
                      map2 (\mat_el top_row ->
                             let prods = map2 (T.*) lft_row top_row
                             let sum = T.sum prods
                             in mat_el T.- sum)
                           mat_row top)
                   mat_blk lft)
            mat_arr top_slice)
      mat_slice lft_per

def pad_to [n] 'a (m: i64) (x: a) (arr: [n]a) : [m]a =
  arr ++ replicate (m - n) x :> [m]a

def lu [m] (block_size: i64) (mat: [m][m]T.t): [m][m]T.t =
  let b = block_size
  let num_blocks = (m+b-1) / b -- rounding up
  let n = b * num_blocks
  -- Maybe pad the input to be a multiple of the block size.
  let padding = n - m
  let mat = if padding != 0
            then map (pad_to n (T.i64 0)) mat ++
                 replicate padding (replicate n (T.i64 0))
            else mat :> [n][n]T.t
  -- transform matrix in [n/b,n/b,b,b] block versions for upper and
  -- lower parts the blocks of the lower part
  let matb =
    tabulate_3d num_blocks num_blocks b
                (\i_b j_b i -> tabulate b (\j -> mat[i_b*b+i, j_b*b + j]))

  let matb = loop(matb) for step < ((n / b) - 1) do
      -- 1. compute the current diagonal block
      let diag = lud_diagonal(matb[step,step]) in

      -- 2. compute the top  perimeter
      let row_slice = matb[step,step+1:num_blocks]
      let top_per_irreg = lud_perimeter_upper(diag, row_slice)

      -- 3. compute the left perimeter and update matrix
      let col_slice = matb[step+1:num_blocks,step]
      let lft_per_irreg = lud_perimeter_lower(diag, col_slice)

      -- 4. compute the internal blocks
      let inner_slice = matb[step+1:num_blocks,step+1:num_blocks]
      let internal = lud_internal(top_per_irreg, lft_per_irreg, inner_slice)

      -- 5. update matrix in place
      let matb[step,step] = diag
      let matb[step, step+1:num_blocks] = top_per_irreg
      let matb[step+1:num_blocks, step] = lft_per_irreg
      let matb[step+1:num_blocks, step+1:num_blocks] = internal
      in matb

  let last_step = (n / b) - 1 in
  let matb[last_step,last_step] =
          lud_diagonal( matb[last_step, last_step] )

  let ret_padded = map (\i_ind  ->
                        map  (\j_ind  ->
                              let (ii, jj) = (i_ind/b, j_ind/b)
                              let ( i,  j) = (i_ind - ii*b, j_ind - jj*b)
                              in  matb[ii,jj,i,j]
                             ) (iota n)
                       ) (iota n)
  in take m (map (take m) ret_padded)
}
