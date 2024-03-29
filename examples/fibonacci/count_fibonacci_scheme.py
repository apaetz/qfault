'''
Created on 2012-08-11

@author: adam
'''
from fibonacci_scheme import FibonacciSchemeAP09, FibonacciSchemeSyndrome
from qfault.counting import count_parallel
from qfault.noise import NoiseModelXSympy, NoiseModelZSympy, NoiseModelXZSympy
from qfault.qec.error import Pauli, xType
from qfault.util.plotting import plotList
import logging


def PlotPaccept(fibonacci, epsilons, j_max):
    paccept = []
    for j in range(1,j_max+1):
        paccept.append([fibonacci.PrAccept(j, e) for e in epsilons])
        
#    print epsilons
#    print 'p(j|j-1):'
#    for j, pj in enumerate(paccept):
#        print j+1, pj

    plotList(epsilons, paccept, labelList=('1','2','3','4','5'), xLabel=r'$\epsilon$', yLabel=r'$p(j|j-1)$', legendLoc='lower left', filename='fibonacci-paccept',
             xscale='log')

def PlotEpsilonCSS(fibonacci, epsilons, j_max):
    e_css = []
    for j in range(1,j_max+1):
        e_css.append([fibonacci.EpsilonCSS(j, e) for e in epsilons])
        
    plotList(epsilons, e_css, labelList=('1','2','3','4','5'), xLabel=r'$\epsilon$', yLabel=r'$\epsilon_{css}$', legendLoc='lower left', filename='fibonacci-e_css',
             xscale='log', yscale='log')
    
def PlotPrBadSBT(epsilons, j_max):
    fibonacci = FibonacciSchemeSyndrome()
    pbad = []
    for j in range(1,j_max+1):
        print j, len(fibonacci._SubblockTeleportComponent(xType, j).locations(Pauli.Y))
        pbad.append([fibonacci._SubblockTeleportPrBad(xType, j)(e/15.) for e in epsilons])
        
    plotList(epsilons, pbad, labelList=('1','2','3','4','5'), xLabel=r'$\epsilon$', yLabel=r'Pr[bad] (sub-block)', legendLoc='lower left', filename='fibonacci-pbad-sbt',
             xscale='log', yscale='log')   
        
        
def PlotPrBadBP1(epsilons):
    fibonacci = FibonacciSchemeSyndrome()
    pbad = []
    pr_bad = fibonacci.BP1().prBad(fibonacci._noise_models[Pauli.Y], Pauli.Y)
    for j in range(1):
        pbad.append([pr_bad(e/15.) for e in epsilons])
        
    plotList(epsilons, pbad, labelList=('1','2','3','4','5'), xLabel=r'$\epsilon$', yLabel=r'Pr[bad] (1-BP)', legendLoc='lower left', filename='fibonacci-pbad-bp1',
             xscale='log', yscale='log')
        


if __name__ == '__main__':
#    util.cache.enableFetch(False)
#    util.cache.enableMemo(False)
    
    count_parallel.setPool(count_parallel.DummyPool())
    
#    logging.getLogger('counting.component').setLevel(logging.WARNING)
#    logging.getLogger('counting.component').setLevel(logging.DEBUG)
    
    logging.getLogger('scheme.fibonacci').setLevel(logging.DEBUG)
    
    noise_models = {Pauli.X: NoiseModelXSympy(),
                              Pauli.Z: NoiseModelZSympy(),
                              Pauli.Y: NoiseModelXZSympy(),
                             }
    kBPT = {Pauli.Y: 1}
#    scheme = FibonacciScheme(kBPT)
    
#    count = scheme.count()
#    print count

#    e = Epsilon(xType, 2, 2, 2, 0.001)
#    print e

#    counting_nm = {Pauli.Y: CountingNoiseModelXZ()}
#    bp1 = BP1Max({Pauli.Y: 1})
#    result = bp1.count(counting_nm, Pauli.Y)
#    print result.counts
#    raise
    
    epsilons = [1e-4, 2e-4, 4e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 1e-2]
#    epsilons = (1e-3,)
    j_max = 5
    
    fib_ap09 = FibonacciSchemeAP09(epsilon_scale=8/15., disable_sbt=[])
    fib_syndrome = FibonacciSchemeSyndrome()
    
#    print fib_syndrome.EpsilonCSS(2, 1e-4/15)

#    print fib_ap09.EpsilonCSS(2, 1e-4/15)

#    print BP1LevelJ(kBPT, BP1Max(kBPT), 1).count(noise_models, Pauli.Y).counts
#    print BP1LevelJ(kBPT, BP1Max(kBPT), 2).count(noise_models, Pauli.Y).counts
#    
#    raise
    
#    print fib_syndrome.PrAccept(1, 1e-3)
#    PlotPaccept(fib_syndrome, epsilons, j_max)
    PlotPrBadSBT(epsilons, j_max)
#    PlotPrBadBP1(epsilons)
#    PlotEpsilonCSS(fib_syndrome, epsilons, j_max)
    
#    p_syndrome = []
#    p_ap09 = []
#    e_syndrome = []
#    e_ap09 = []
#    for j in range(2,3):
#        p_syndrome.append([fib_syndrome.PrAccept(j,e/15) for e in epsilons])
#        e_syndrome.append([fib_syndrome.EpsilonCSS(j, e) for e in epsilons])
#        p_ap09.append([fib_ap09.PrAccept(j,8/15. * e) for e in epsilons])
#        e_ap09.append([fib_ap09.EpsilonCSS(j, e) for e in epsilons])
#        
#    print epsilons
##    print 'p(j|j-1):'
#    for j in range(len(p_syndrome)):
##        print j+1, p_syndrome[j], p_ap09[j]
#        print j+1, e_syndrome[j], e_ap09[j]


#    generator = SyndromeKeyGenerator(ed422.ED412Code(gaugeType=error.xType), None)
#    xkey = generator.get_key(PauliError(4, 1, 0))
#    print xkey
#        
##    print fib_syndrome.Epsilon1(xkey, 1, 1, 1e-3/15)
#
#    print fib_syndrome.BlockTeleport(xType, 1, 3, 1e-3/15)
#    print fib_ap09.BlockTeleport(xType, 1, 3, 1e-3 * 8/15.)
    
#    print fib_ap09.Epsilon(zType, 0, 1, 1, 1e-3 * 8/15.)
#    print fib_ap09.Epsilon(zType, 1, 1, 1, 1e-3 * 8/15.)
    
#        
#
#    plotList(epsilons, paccept, labelList=('1','2','3','4','5'), xLabel=r'$\epsilon$', yLabel=r'$p(j|j-1)$', legendLoc='lower left', filename='fibonacci-paccept',
#             xscale='log')

#    e0 = 0.001
#    print Epsilon(xType, 0, 0, 0, e0)
#    print Epsilon(xType, 0, 1, 0, e0)
#    print Epsilon(xType, 1, 1, 0, e0)
#    print Epsilon(xType, 1, 1, 1, e0)

#    print f_s(xType, 2, 2, e0)
