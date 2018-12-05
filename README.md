# fastLiquidAssociation
Forked from gundt/fastLiquidAssociation

## Changes:
1. Allow `threads = 1` for `fastMLA`
2. Allow `nvec` as a list with `z` and `x` for fastMLA, where `z` is the original `nvec`, `x` is the indices of first group.   
   Then the pairs to examine would be one mate from `x` and the other from the rest (all - `x` - `z`). Additionally, `z` will  
   not be involved in the pairing any more.
3. Allow group "Z" as discrete values, in this case, `Stein lemma` will not apply, `E(g'(z))` will be calulated manully.
