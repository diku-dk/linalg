-- Random matrices are almost certainly invertible, so we should be
-- good here.

import "../lib/github.com/diku-dk/linalg/linalg"

module linalg_f32 = mk_linalg f32

-- ==
-- entry: bench_batch_inv
-- random input { [1][1000][1000]f32 }
-- random input { [1000][128][128]f32 }
-- random input { [10000][8][8]f32 }
-- random input { [10000][16][16]f32 }
entry bench_batch_inv ms = #[unsafe] map linalg_f32.inv ms

-- ==
-- entry: bench_inv
-- random input { [1000][1000]f32 }
entry bench_inv m = #[unsafe] linalg_f32.inv m
