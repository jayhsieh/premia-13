
import sys;
import re
import string
import os.path 

def get_var (DATA, name):
    inside = 0;
    val = list ();
    for l in DATA:
        if ( not inside and re.match (r'.*' +  name, l) ):
            entry = re.sub (r'.*' + name + r'[ =]*' , '', l)
            entry = re.sub(r'[ ]*\\$', '', entry)
            val.append (entry)
            inside = 1
        elif (inside):
            entry = re.sub(r'[ ]*\\$', '', l)
            val.append (entry)
        if ( inside and not re.match (r'.*\\$', l) ):
            inside = 0;
            return val
    return val

def get_all_vars (fic):
    f = open(fic, "r")
    DATA = f.readlines()
    f.close()

    DATA = map (lambda x: re.sub(r'\n$', '', x), DATA)
    subdir = get_var (DATA, 'SUBDIRS')
    sources = get_var (DATA, 'SOURCES')
    libname = get_var (DATA, 'noinst_LTLIBRARIES')

    # Join all entries and simplify spaces
    subdir = re.sub (r'[ ]+', ' ', string.join (subdir))
    sources = re.sub (r'[ ]+', ' ', string.join (sources))

    # Extract libname
    (libname, ext) = os.path.splitext (libname[0])
    
    return (subdir, libname, sources)

def treat (fic):
    print "Reading " + fic + "...",
    (subdir, libname, sources) = get_all_vars (fic)
    dir = os.path.dirname (fic)
    cmakefile = os.path.join (dir, 'CMakeLists.txt')
    print "Writing " + cmakefile + "...",
    out = open (cmakefile, 'w')
    # out = sys.stdout

    if ( subdir != '' ):
        for e in string.split(subdir):
            out.write ("add_subdirectory (" + e + ")\n")
    if ( sources != '' ):
        out.write ("add_library (" + libname + " OBJECT\n")
        for e in string.split(sources):
            out.write ( '\t' + e + '\n' ) 
        out.write (")\n")
    out.close ()
    print "Done."

if __name__ == '__main__': 
    treat (sys.argv[1])

