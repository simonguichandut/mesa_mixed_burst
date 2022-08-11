from utils import *

def h1_spline(column_depth, iso_data, x1, x2, Npts):
    xx = column_depth
    yy = iso_data['h1']

    imax = np.argmin(abs(xx-eval(x2)))

    if x1=='auto':
        imin = list(xx).index(xx[yy<max(yy)-0.005][0])-1
    else:
        imin = np.argmin(abs(xx-eval(x2)))

    xx = xx[imin:imax]
    yy = yy[imin:imax]

    # Sub-sample in logspace
    xxsub = np.logspace(np.log10(xx[0]), np.log10(xx[-1]), Npts)
    # These points don't actually exist, take the closest ones
    isub = [np.argmin(abs(xx-xi)) for xi in xxsub]

    # Interpolate a spline
    # The type of spline (Pchip) is a monotonic, piecewise cubic hermite polynomials
    spl = interpolate.PchipInterpolator(xx[isub], yy[isub])

    # New points
    all_isos = list(iso_data.keys())
    new_iso_data = {key:np.copy(iso_data[key]) for key in all_isos}   # new dictionary preserves isotope order
    new_iso_data['h1'] = np.concatenate((iso_data['h1'][:imin], spl(xx), iso_data['h1'][imax:]))
    
    # Then leave all the metals as is, and helium is Y=1-X-Z
    
    Z = np.zeros(len(new_iso_data['h1']))
    for iso in all_isos:
        if iso in ('h1','he4'):
            pass
        else:
            new_iso_data[iso] = iso_data[iso]
            Z += iso_data[iso]

    new_iso_data['he4'] = 1 - new_iso_data['h1'] - Z

    return new_iso_data, imin, imax


def write_new_model(iso_data, old_modfile, new_modfile):
    with open(old_modfile, "r") as fold:
        with open(new_modfile, "w") as fnew:
        
            # Header data stays the same
            prev_line = ''
            stop_next = False
            for line in fold:
                fnew.write(line)
                if stop_next:
                    break
                if line == '\n' and 'time' in prev_line:
                    stop_next = True
                prev_line = line
    
    
            # Data: replace abundance data
            ih1 = model.bulk_names.index('h1')
            for k,line in enumerate(fold):
            
                if line == '\n': # end of data
                    fnew.write(line)
                    break
                
                new_line = []
                for i,val in enumerate(line.split()):
                
                    # File doesn't change for all columns not abundance
                    if i<ih1:
                        new_val_str = val
                        # if eval(val)>=0:
                        #     new_val_str = ' ' + new_val_str # extra whitespace so aligns with negative values
    
                    else:
                        new_val = new_iso_data[list(new_iso_data.keys())[i-ih1]][k]
                        new_val_str = ("%.17e"%new_val).replace('e',"E")
    
                    if eval(new_val_str)>=0:
                        new_val_str = ' ' + new_val_str # extra whitespace so aligns with negative values
    
                    new_line.append(new_val_str)
                
                whitespace = 4-len(str(k+1))
                fnew.write(' '*whitespace + '   '.join(new_line)+'\n')
    
            # previous model data
            for line in fold:
                fnew.write(line)


        
# Command line call
parser = argparse.ArgumentParser(description="Smooth out abundance profiles")
parser.add_argument('-i','--input', type=str, help='Input model file', default=None)
parser.add_argument('-o','--output', type=str, help='Output model file (default is <input_model>_smoothed.mod', default=None)
parser.add_argument('-s','--show', action='store_true', help='show plot dont save')

parser.add_argument('-m','--method', type=str, help='method for smoothing. Currently only h1_spline', default='h1_spline')
parser.add_argument('-x1', type=str, help='lower boundary (column depth) of the smoothing interval. auto finds this automatically as the first point where h1 dips below constant', default='auto')
parser.add_argument('-x2', type=str, help='upper boundary (column depth) of the smoothing interval', default=None)
parser.add_argument('-n','--Npoints', type=int, help='number of points to subsample in the interval. Less points means more smoothing', default=10)

if __name__ == "__main__":
    args = parser.parse_args()

    model = mr.MesaData(args.input)

    # Get column depth array
    column_depth = integrate.cumtrapz(model.dq/model.R**2) * model.xmstar/(4*np.pi)
    column_depth = np.insert(column_depth, 0, 1e-99)

    # Generate abundances dict
    ih1 = model.bulk_names.index('h1')
    iso_data = { key:model.bulk_data[key] for key in list(model.bulk_data.keys())[ih1:] }    

    # Smooth
    if args.method == 'h1_spline':     
        new_iso_data,imin,imax = h1_spline(column_depth, iso_data, x1=args.x1, x2=args.x2, Npts=args.Npoints)

    # Plot or save
    if args.show:
        fig,ax = plt.subplots(1,1)
        ax.loglog(column_depth, iso_data['h1'], color='b', ls='-', lw=0.8)
        ax.loglog(column_depth, iso_data['he4'], color='r', ls='-', lw=0.8)
        ax.loglog(column_depth, new_iso_data['h1'], color='b', ls='--', lw=0.8)
        ax.loglog(column_depth, new_iso_data['he4'], color='r', ls='--', lw=0.8)
        ax.set_xlim([1e2,1e10])
        ax.set_ylim([1e-3,1.1])
        ax.axvline(column_depth[imin], color='k', ls='--', lw=0.5)
        ax.axvline(column_depth[imax], color='k', ls='--', lw=0.5)
        plt.show()

    else:
        if args.output is None:
            new_modfile = model.file_name.replace(".mod", "_smoothed.mod")
        else:
            new_modfile = args.output

        write_new_model(new_iso_data, args.input, new_modfile)
        
        # Also write the code that was used to smooth this out (in case this present code changes)
        info_file = new_modfile.replace(".mod", "_info.txt")
        import inspect
        from datetime import datetime
        with open(info_file, 'w') as f:
            f.write("Generated on " + datetime.today().strftime('%Y-%m-%d %H:%M:%S') + "\n")
            f.write("The following arguments were used: \n")
            for arg in vars(args):
                f.write(arg + " : " + str(getattr(args,arg)) + "\n")

            f.write("\nThe code for the smoothing function at this time was: \n")
            f.write('------------------------------------------------------------------------------\n')
            f.write(inspect.getsource(h1_spline))
            f.write('------------------------------------------------------------------------------\n')
