-- Random systems of linear equations almost certainly have solutions,
-- so we should be good here.

import "../lib/github.com/diku-dk/linalg/linalg"

module linalg_f32 = mk_linalg f32

-- ==
-- entry: bench_ols
-- random input { [100][100]f32 [100]f32 }
-- random input { [200][200]f32 [200]f32 }
-- random input { [300][300]f32 [300]f32 }
-- random input { [400][400]f32 [400]f32 }
entry bench_ols m b = #[unsafe] linalg_f32.ols m b
