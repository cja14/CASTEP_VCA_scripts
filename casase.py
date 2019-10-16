import ase.io
import sys

# To redirect the stdout/stderr (took ages to get this to work!!!)
class NullDevice():
    def write(self, s):
        pass

def casread(casfile):
    oldtargetout = sys.stdout
    oldtargeterr = sys.stderr
    try:
        sys.stdout = NullDevice()
        sys.stderr = NullDevice()
        atoms = ase.io.read(casfile)
        sys.stdout = oldtargetout
        sys.stderr = oldtargeterr
        return atoms
    except:
        sys.stdout = oldtargetout
        sys.stderr = oldtargeterr
        raise sys.exc_info()
