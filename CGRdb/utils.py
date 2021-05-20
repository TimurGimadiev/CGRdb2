organic_set = {'B', 'C', 'N', 'O', 'Si', 'P', 'S', 'Se', 'F', 'Cl', 'Br', 'I'}


def validate_molecule(mol):
    if int(mol):
        return False
    for n, a in mol.atoms():
        if a.atomic_symbol not in organic_set:
            return False
    return True


__all__ = ['validate_molecule']
