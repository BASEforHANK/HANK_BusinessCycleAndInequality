# HANKEstim Toolbox

Toolbox for estimating medium-scale HANK models based on "Shocks, Frictions, and Inequality in US Business Cycles" by Christian Bayer, Benjamin Born and Ralph Luetticke ([link to paper](https://www.benjaminborn.de/publication/bbl_inequality_2020/)).

Note that the paper is currently under revision and the updated toolbox does not replicate the results in the linked version of the paper.

There is also currently an issue with MKL on macOS. Therefore, we use OpenBLAS for non-LINUX machines, but it's important that you at use Julia 1.7.1 or newer. Note, that MKL offers considerable speed-ups.

To get started, clone or download the repository. The full documentation can be found [here](http://www.benjaminborn.de/HANK_BusinessCycleAndInequality/build/) or locally in `docs/build`.
