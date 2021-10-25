import pwem.objects.data as data

class TopolStruct(data.EMFile):
    """Represents a Topol.top file. """

    def __init__(self, filename=None, **kwargs):
        data.EMFile.__init__(self, filename, **kwargs)
        self._origin = None

class PDB_test(data.AtomStruct):
    """ TEST """
    def __init__(self, filename=None, pseudoatoms=False, **kwargs):
        data.AtomStruct.__init__(self, filename, pseudoatoms, **kwargs)

class GroFile(data.AtomStruct):
    """Represents a Gro file """
    def __init__(self, filename=None, pseudoatoms=False, **kwargs):
        data.AtomStruct.__init__(self, filename, pseudoatoms, **kwargs)

class TopolFile(TopolStruct):
    def __init__(self, filename=None, **kwargs):
        TopolStruct.__init__(self, filename, **kwargs)