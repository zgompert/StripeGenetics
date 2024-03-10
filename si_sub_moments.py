import moments
import scipy
import numpy
import matplotlib
import matplotlib.pyplot as pyplot
import sys
import re
import pylab

## set population ids, 2n sample sizes, and vcf file
#n1 = 238*2;
#n1 = 57*2;
#n2 = 602*2;
p1 = "R17";
#p1 = "refugio";
p2 = "main";
#p2 = "santaynez";
vv = "filtered_vars_noch8.vcf";

## downsample to 50
n1 = 50
n2 = 50

ns = numpy.array([int(n1), int(n2)])

## read data from filtered vcf and create downsampled/projected 2d sfs 
dd = moments.Misc.make_data_dict_vcf(vv,"/uufs/chpc.utah.edu/common/home/gompert-group4/projects/timema_SV_balance/demog/idsPlus.txt")
#dd = moments.Misc.make_data_dict_vcf(vv,"/uufs/chpc.utah.edu/common/home/gompert-group4/projects/timema_SV_balance/demog/ids.txt")
fs = moments.Spectrum.from_data_dict(dd,[p1, p2],projections=ns, polarized=False)
fs_proj = fs ## not currently projecting down

of = "dd_mom_"+p1+"_"+p2+"_si.txt"
pngf = "dd_mom_"+p1+"_"+p2+"_si.png"

ofile = open(of,"w")

ns =  fs_proj.sample_sizes

## define strict isolation (SI) model, includes exponential growth
def SI(params, ns):
    s,nu1,nu2,T = params
    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate(nu_func, T, m=numpy.array([[0, 0], [0, 0]]))
    return fs

## define strict isolation (SI) model, split into 2 pops 
def SI(params, ns):
    nu1,nu2,T = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T, m=numpy.array([[0, 0], [0, 0]]))
    return fs

func = SI
# Parameters are: (s, nu1, nu2, and T)
#upper_bound = [.9, 40, 40, 6]
#lower_bound = [.1, 0.1, 0.1, 0]

# Parameters are: (nu1, nu2, and T)
upper_bound = [100, 100, 10]
lower_bound = [0.1, 0.1, 0]

# This is our initial guess for the parameters, which is somewhat arbitrary.
#p0 = [0.5, 0.8, 1.2, 0.4]
p0 = [0.8, 1.2, 0.4]
np = len(p0)
p0 = moments.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
# powell method for apporximate solution
popt = moments.Inference.optimize_log_powell(p0, fs_proj, func, lower_bound=lower_bound, upper_bound=upper_bound,verbose=5, maxiter=10)

# optimization
ofile = open(of,"w")
print(ofile)
ofile.write("{0}".format(fs_proj.Fst()))
ofile.write("\n")

## initial defaults
ll_opt = -99999999999
theta = 1.0
mpopt = p0
popt = p0
for x in range(10):
    theta = 1.0
    popt = moments.Misc.perturb_params(popt, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
    # powell method for apporximate solution
    popt = moments.Inference.optimize_log_powell(popt, fs_proj, func, lower_bound=lower_bound, upper_bound=upper_bound,verbose=5, maxiter=30)
    # BFGS for final fit
    popt = moments.Inference.optimize_log(popt, fs_proj, func, lower_bound=lower_bound, upper_bound=upper_bound,verbose=5, maxiter=30)
    model = func(popt, ns)
    ll_model = moments.Inference.ll_multinom(model,fs_proj)
    theta = moments.Inference.optimal_sfs_scaling(model,fs_proj)
    ## write current
    ofile.write("{0}".format(ll_model))
    ofile.write(" {0}".format(theta))
    for a in range(np):
        ofile.write(" {0}".format(popt[a]))
        
    ofile.write("\n")

    if(ll_model > ll_opt):
        ll_opt = ll_model
        mpopt = popt		
        theta_opt = moments.Inference.optimal_sfs_scaling(model,fs_proj)

## write max
ofile.write("{0}".format(ll_opt))
ofile.write(" {0}".format(theta_opt))
for a in range(np):
    ofile.write(" {0}".format(mpopt[a]))

ofile.write("\n")
ofile.close()


pylab.figure(figsize=(8,6))
moments.Plotting.plot_2d_comp_multinom(model,fs_proj,vmin=1,resid_range=50,show=False)
pylab.savefig(pngf, dpi=400)

