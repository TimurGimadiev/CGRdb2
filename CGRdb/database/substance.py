from pony.orm import Database, PrimaryKey, Required, Set, Optional as PonyOptional, Json
from .config import Entity
from CGRtools import MoleculeContainer
from .molecule import NonOrganic, Molecule
from typing import Tuple, Iterable, Union, Dict, List, Optional
from functools import cached_property, reduce
from operator import or_

orgainc_set = ["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"]


class Substance(Entity):
    id = PrimaryKey(int, auto=True)
    reactions = Set('ReactionSubstance')
    components = Set('SubstanceStructure')
    name = PonyOptional(str)

    def __init__(self, substance: Optional[Iterable[Tuple[MoleculeContainer, Optional[float]]]] = None, /):
        """
        :param substance:
        """
        super().__init__()
        if substance is not None:
            if isinstance(substance, Iterable):
                raise ValueError("Iterable type consisted of tuples should be provided")
            substance = list(substance)
            for mol, frac in substance:
                if not isinstance(mol, MoleculeContainer):
                    raise ValueError('First element of tuple should be CGRtools.MoleculeContainer')
                if frac is not None and not isinstance(frac, float):
                    raise TypeError('molar_fraction should be of float type')

            for mol, frac in substance:
                ...
                SubstanceStructure(structure, self, molar_fraction=frac, mapping=mapping)

    @cached_property
    def structure(self) -> MoleculeContainer:
        return reduce(or_, (x.structure for x in self.components))


class SubstanceStructure(Entity):
    id = PrimaryKey(int, auto=True)
    molar_fraction = PonyOptional(float)
    _mapping = PonyOptional(Json, column='mapping', nullable=True)
    substance = Required('Substance')
    molecule = PonyOptional('Molecule')
    non_organic = PonyOptional('NonOrganic')

    def __init__(self, structure: Union[Molecule, NonOrganic], substance: 'Substance', /, *,
                 molar_fraction: Optional[float] = None,
                 mapping: Union[Dict[int, int], List[Tuple[int, int]], None] = None):
        if molar_fraction is not None and not isinstance(molar_fraction, float):
            raise TypeError('molar_fraction should be of float type')

        if isinstance(mapping, dict):
            if not all(isinstance(x, int) for x in mapping.items() for x in x):
                raise TypeError('Mapping is dict that contains following structure: '
                                '{int (atom of substance) : int (atom of NonOrganic or Molecule)}')
            mapping = list(mapping.items())
        elif isinstance(mapping, (tuple, list)):
            if not all(isinstance(x, (tuple, list)) and len(x) == 2 and all(isinstance(x, int) for x in x)
                       for x in mapping):
                raise TypeError('For the tuple pairs of integers following constructions should be applied: '
                                '[(int (atom of substance), int (atom of NonOrganic or Molecule)),...]')
        elif mapping is not None:
            raise TypeError('Mapping of Substance to NonOrganic or Molecule expected')

        if isinstance(structure, Molecule):
            molecule = structure
            non_organic = None
        elif isinstance(structure, NonOrganic):
            molecule = None
            non_organic = structure
        else:
            raise TypeError('At least one either Molecule or NonOrganic required')
        super().__init__(molar_fraction=molar_fraction, molecule=molecule, _mapping=mapping,
                         non_organic=non_organic, substance=substance)

    @cached_property
    def structure(self):
        if self.molecule:
            return self.molecule.canonic_structure.structure.remap(self.mapping, copy=True)
        else:
            return self.non_organic.structure.remap(self.mapping, copy=True)

    @cached_property
    def mapping(self):
        return dict(self._mapping) if self._mapping else {}


__all__ = ['SubstanceStructure', 'Substance']
