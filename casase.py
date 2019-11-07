import ase.io
import sys


class NullDevice():
    """ Blank output stream to redirect the stdout/stderr to """
    def write(self, s):
        pass


def casread(casfile):
    """ ase.io.read() except with no ouput if CASTEP not linked properly """
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
        """ Bare exception not a problem since error gets raised anyway """
        sys.stdout = oldtargetout
        sys.stderr = oldtargeterr
        raise sys.exc_info()
