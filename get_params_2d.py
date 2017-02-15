import re, os, sys
import numpy as np

os.chdir("data/")

params = []
for d in os.walk(".").next()[1]:
    m = re.match(r"W_R(\d+?)_amp2d_(\d+?)", d)
    if not m:
        continue

    print "Processing directory %s ..." % d
    R = int(m.group(1))
    ind = int(m.group(2))
    ## get parameters
    data = np.load(d+"/qarray_2d.npz")
    h = data['h']

    if 'chi_array' in data: # case 1: using chi, tth
        chi_array = data['chi_array']
        tth_array = data['tth_array']
        chi_range = (chi_array[-1]-chi_array[0])/np.pi*180.
        tth_range = (tth_array[-1]-tth_array[0])/np.pi*180.
        
        sys.path[0] = d
        from W_parameters import a0, D
        h = h / (2.*np.pi/a0)

        params.append((R, ind, "%.1f"%(D*1./R), "%d%d%d"%(h[0], h[1], h[2]), \
                       "TC", "%.2f\t"%chi_range, "%.2f\t"%tth_range))
        for looptype in ["int", "vac"]:
            os.system("cp %s/intensity_2d_%s_R%d.pdf "%(d, looptype, R) + \
                "plots/R%d_h%d%d%d_chi%.1f_tth%.1f_intensity_2d_%s.pdf"%(\
                 R, h[0], h[1], h[2], chi_range, tth_range, looptype))

        del sys.modules['W_parameters']

    elif 'hk_array' in data: # case 2: slicing hkl space on (1 -1 0) plane
        hk_array = data['hk_array']
        l_array = data['l_array']
        hk_range = (hk_array[-1]-hk_array[0])/(len(hk_array)-1)*len(hk_array)
        l_range = (l_array[-1]-l_array[0])/(len(l_array)-1)*len(l_array)
        
        sys.path[0] = d
        from W_parameters import a0, D
        h = h / (2.*np.pi/a0)

        params.append((R, ind, "%.1f"%(D*1./R), "%d%d%d"%(h[0], h[1], h[2]), \
                       "HKL", "%.3f\t"%hk_range, "%.3f\t"%l_range))
        for looptype in ["int", "vac"]:
            os.system("cp %s/intensity_2d_%s_R%d.pdf "%(d, looptype, R) + \
                "plots/R%d_h%d%d%d_hk%.3f_l%.3f_intensity_2d_%s.pdf"%(\
                 R, h[0], h[1], h[2], hk_range, l_range, looptype))

        del sys.modules['W_parameters']


params = sorted(params, key=lambda tup: (tup[0], tup[1]))
## write to file
with open("parameters.txt", 'w') as f:
    f.write("R\tind\tD/R\th\ttype\tchi/hk rg\ttth/l rg\n")
    for tup in params:
        f.write('\t'.join(str(item) for item in tup)+'\n')

os.system("cat parameters.txt")
