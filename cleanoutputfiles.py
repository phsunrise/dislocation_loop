import os, glob

os.system("ls *.err")
val = raw_input("Clean these and associated files? ")
if val in ['y', 'Y', 'yes']:
    for filename_err in glob.glob("*.err"):
        filename = filename_err[:-4]
        os.system("rm %s.err"%filename)
        os.system("rm %s.out"%filename)
        os.system("rm %s.sbatch"%filename)
        print "Finished %s!" % filename
else:
    print "Aborted!"
