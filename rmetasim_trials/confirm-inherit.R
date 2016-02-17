confland <- function () 
{
    rland <- NULL
    rland <- landscape.new.empty()
    rland <- landscape.new.intparam(rland, h = 2, s = 2)
    rland <- landscape.new.switchparam(rland, mp = 0)
    rland <- landscape.new.floatparam(rland)
    S <- matrix(c(1, 0,
                  0, 1), byrow = TRUE, nrow = 2)
    R <- matrix(c(0, 10,
                  0, 0), byrow = TRUE, nrow = 2)
    M <- matrix(c(0, 0,
                  0, 1), byrow = TRUE, nrow = 2)
    rland <- landscape.new.local.demo(rland, S, R, M)
    S <- matrix(rep(0, 16), nrow = 4)
    R <- matrix(rep(0, 16), nrow = 4)
    M <- matrix(rep(0, 16), nrow = 4)
    rland <- landscape.new.epoch(rland, S = S, R = R, M = M, 
        carry = c(1000, 1000))
    rland <- landscape.new.locus(rland, type = 0, ploidy = 1, 
        mutationrate = 0.00, transmission = 1, numalleles = 20)
    rland <- landscape.new.locus(rland, type = 2, ploidy = 2, 
        mutationrate = 0.007, transmission = 0, numalleles = 6, 
        allelesize = 75)
    rland <- landscape.new.individuals(rland, c(0, 2, 0, 2))
    rland
}
