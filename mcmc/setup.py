#!/usr/bin/python
import sys, os, glob, shutil

try:
    prjdir = sys.argv[1]
    objdir = sys.argv[2]
except:
    print "Usage:", sys.argv[0], "<project name> <object directory>"
    sys.exit(1)

# treat the object directory
if os.path.islink(objdir):
    print objdir, "is a link."
    sys.exit(1)
elif os.path.isfile(objdir):
    print objdir, "is a file."
    sys.exit(1)
elif os.path.isdir(objdir):
    print objdir, "is an existing directory."
else:
    print objdir, "does not exist."
    print "Creating a new directory ..."
    os.makedirs(objdir)

rootdir = os.getcwd()

prjdir = os.path.join(rootdir, "src/project", prjdir)
auxdir = os.path.join(rootdir, "src/aux")
mcmcdir = os.path.join(rootdir, "src/mcmc")

print "Creating links to project files ..."
os.chdir(prjdir)
source = glob.glob("*.c")+glob.glob("*.cc")+glob.glob("*.f90")
header = glob.glob("*.h")+glob.glob("*.H")
os.chdir(rootdir)
for files in source+header:
    if os.path.exists(os.path.join(objdir, files)):
        os.remove(os.path.join(objdir, files))
    os.symlink(os.path.join(prjdir, files), os.path.join(objdir, files))

print "Creating links to auxiliary files ..."
os.chdir(auxdir)
source = glob.glob("*.c")+glob.glob("*.cc")
header = glob.glob("*.h")+glob.glob("*.H")
os.chdir(rootdir)
for files in source+header:
    if os.path.exists(os.path.join(objdir, files)):
        os.remove(os.path.join(objdir, files))
    os.symlink(os.path.join(auxdir, files), os.path.join(objdir, files))

print "Creating links to mcmc files ..."
os.chdir(mcmcdir)
source = glob.glob("*.c")+glob.glob("*.cc")
header = glob.glob("*.h")+glob.glob("*.H")
os.chdir(rootdir)
for files in source+header:
    if os.path.exists(os.path.join(objdir, files)):
        os.remove(os.path.join(objdir, files))
    os.symlink(os.path.join(mcmcdir, files), os.path.join(objdir,files))
if os.path.exists(os.path.join(objdir, "Makefile")):
    os.remove(os.path.join(objdir, "Makefile"))
os.symlink(os.path.join(mcmcdir, "Makefile"), os.path.join(objdir,"Makefile"))

print "End"
