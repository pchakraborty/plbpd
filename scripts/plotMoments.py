#!/usr/bin/env python
# plotMoments.py


"""
Input:   hdf5 file containing velocity and density data, Y (plot is for
         x-section y=Y) and scale (for quiver)
Output:  matplotlib plots of moments (density, velocity, temperature)
"""

##############################################################

def processCommandLineArgs():
    """
    returns a list containing the command line args
    1st element: input hdf5 file
    2nd element: Y (x-section y=Y)
    3rd element: scale (for quiver)
    """

    import optparse, sys
    
    usage = "for usage, try: %prog -h"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-f", dest="file", default="file.h5",
                      help="input hdf5 file [default:file.h5]", metavar="FILE")
    parser.add_option("-y", type="int", help="x-section is at y=Y", metavar="Y")
    parser.add_option("-s", type="float", default=1.,
                      help="scale (for quiver) [default=1]", metavar="scale")
    parser.add_option("-w", type="string", default="all",
                      help="what to plot: density, velocity or temperature [default:all]")
    parser.add_option("--manual_label", action="store_true", dest="clabel",
                      help="add labels to isotherms manually")
    parser.add_option("-o", dest="ofile", default="vel.png",
                      help="output vel profile file [default=vel.png]",
                      metavar="FILE")

    (opts, args) = parser.parse_args()
    if opts.clabel != True:
        opts.clabel = False
    
    # exit if incorrect num of options are specified
    if len(args)>0:
        parser.error("additional arg(s) specified")

    if (opts.y==None):
        parser.error("incorrect options specified")

    # store the args in a buffer (list)
    commLineArgs = list()
    commLineArgs.append(opts.file)
    commLineArgs.append(opts.y)
    commLineArgs.append(opts.s)
    commLineArgs.append(opts.w)
    commLineArgs.append(opts.clabel)
    commLineArgs.append(opts.ofile)

    # print command line options provided
    print ''
    print 'COMMAND LINE ARGS:'.rjust(30)
    print 'input hdf5 file:'.rjust(30), opts.file
    print 'input Y (x-section is at y=Y):'.rjust(30), opts.y
    print 'input scale (quiver):'.rjust(30), opts.s
    print 'what to plot:'.rjust(30), opts.w
    print 'output vel profile file:'.rjust(30), opts.ofile
    print ''

    # return the buffer to the calling function
    return commLineArgs

##############################################################

def main():
    import h5py as h5
    import numpy as np
    import pylab as pl

    # read the command line arguments
    tmp = processCommandLineArgs()
    h5file, Y, scale, what2plot, CLABEL, ofile = tmp[0:6]

    # load hdf5 file
    f = h5.File(h5file,'r')

    # density
    if what2plot=="all" or what2plot=="density":
        rho = f['rho']
        nrows = rho.shape[0]
        n2 = rho.shape[1]
        ncols = rho.shape[2]
        shape = (nrows,ncols)
        rhoplane = np.zeros(shape)
        assert(Y>=0); assert(Y<=n2)
        for i in xrange(nrows):
            for j in xrange(ncols):
                rhoplane[i,j]= rho[i,Y,j]
        pl.figure(1)
        pl.imshow(rhoplane,alpha=1.0,origin='lower'), pl.colorbar()
        pl.title("density")

    # temperature (theta)
    if what2plot=="all" or what2plot=="temperature":
        if 'theta' in list(f):
            theta = f['theta']
            nrows = theta.shape[0]
            n2 = theta.shape[1]
            ncols = theta.shape[2]
            shape = (nrows,ncols)
            thetaplane = np.zeros(shape)
            for i in xrange(nrows):
                for j in xrange(ncols):
                    thetaplane[i,j] = theta[i,Y,j]

            pl.figure(2)
            pl.imshow(thetaplane,alpha=1.0,origin='lower'), pl.colorbar()
            pl.title("temperature")

            V = []
            x = 0
            while x<1.0:
                x += 0.05
                V.append(x)
            pl.figure(3)
            contourplot = pl.contour(thetaplane,V)
            pl.clabel(contourplot,inline=1,fontsize=10,manual=CLABEL)
            pl.title("isotherms")

        else:
            print "input file does not contain temperature data"

    # velocity
    if what2plot=="all" or what2plot=="velocity":
        u = f['u']
        nrows = u.shape[0]
        n2 = u.shape[1]
        ncols = u.shape[2]
        assert(u.shape[3]==3)
        shape = (nrows,ncols)
        ux = np.zeros(shape)
        uz = np.zeros(shape)
        for i in xrange(nrows):
            for j in xrange(ncols):
                ux[i,j] = u[i,Y,j,0]
                uz[i,j] = u[i,Y,j,2]
        [x,z] = pl.mgrid[0:ncols, 0:nrows]
        pl.figure(4)
        pl.quiver(x.transpose(),z.transpose(),ux,uz,scale=1/scale)
        pl.axis('tight')
        if(ofile!="vel.png"):
            pl.savefig(ofile)

    pl.ion()
    if(ofile=="vel.png"):
        pl.show()
    
##############################################################

if __name__=="__main__":
    import warnings
    warnings.filterwarnings("ignore")
    main()

##############################################################
