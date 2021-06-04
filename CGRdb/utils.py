from datasketch import MinHashLSH
organic_set = {'B', 'C', 'N', 'O', 'Si', 'P', 'S', 'Se', 'F', 'Cl', 'Br', 'I'}


def validate_molecule(mol):
    if int(mol):
        return False
    for n, a in mol.atoms():
        if a.atomic_symbol not in organic_set:
            return False
    return True


class MinHashLSH(MinHashLSH):

    def _byteswap(self, hs):
        return int.from_bytes(hs.byteswap().data, 'big')

    def _hashed_byteswap(self, hs):
        return self.hashfunc(int.from_bytes(hs.byteswap().data, 'big'))


__all__ = ['validate_molecule', 'MinHashLSH']
