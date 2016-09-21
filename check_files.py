#!/usr/bin/env python
import ROOT
import rat
import os
import fnmatch
import optparse

if __name__ == '__main__':
    """Script runs through all files in the directory given by the -d flag. It tries to open
       the root files to check if they are zombies.
    """

    parser = optparse.OptionParser( usage = "python %prog [flags]")
    parser.add_option("-d", dest="dir_name", help="Name of directory for input files")

    (options, args) = parser.parse_args()

    if not options.dir_name:
        print 'Need directory name!'
        raise Exception
    
    #run through directory
    for root, dirs, filenames in os.walk(options.dir_name):
        for x in filenames:
            if fnmatch.fnmatch(x, "*.root"):
                infile = options.dir_name + x

                #open file
                f = ROOT.TFile(infile)
                if f.IsZombie():
                    print "File %s is a zombie!" % x
