#!/usr/bin/env python
# plot


"""
Input:   hdf5 restart file containing velocity and density data,
         Y (plot is for x-section y=Y)
Output:  matplotlib plots of nodetype at y=Y plane
"""

##############################################################

def processCommandLineArgs():
    """
    returns a list containing the command line args
    1st element: input hdf5 file
    2nd element: Y (x-section y=Y)
    """

    import optparse, sys
    
    usage = "for usage, try: %prog -h"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-f", dest="file", default="file.h5",
                      help="input hdf5 file [default:file.h5]", metavar="FILE")
    parser.add_option("-y", type="int", help="x-section is at y=Y", metavar="Y")

    (opts, args) = parser.parse_args()

    # exit if incorrect num of options are specified
    if len(args)>0:
        parser.error("additional arg(s) specified")

    if (opts.y==None):
        parser.error("incorrect options specified")

    # store the args in a buffer (list)
    commLineArgs = list()
    commLineArgs.append(opts.file)
    commLineArgs.append(opts.y)

    # print command line options provided
    print ''
    print 'COMMAND LINE ARGS:'.rjust(30)
    print 'input hdf5 file:'.rjust(30), opts.file
    print 'input Y (x-section is at y=Y):'.rjust(30), opts.y
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
    h5file, Y = tmp[0:2]

    # load hdf5 file
    f = h5.File(h5file,'r')

    # plot a colormap of type
    type = f['type']
    nrows = type.shape[0]
    n2 = type.shape[1]
    ncols = type.shape[2]
    shape = (nrows,ncols)
    typeplane = np.zeros(shape)
    assert(Y>=0); assert(Y<=n2)
    for i in range(nrows):
        for j in range(ncols):
            typeplane[i,j]= type[i,Y,j]
    pl.figure(1)
    pl.imshow(typeplane,alpha=1.0), pl.colorbar()

    pl.show()
    
##############################################################

if __name__=="__main__":
    import warnings
    warnings.filterwarnings("ignore")
    main()

##############################################################
