from functools import cached_property
from typing import Optional, Union, List, Dict, Tuple

from CGRtools.containers import ReactionContainer
from pony.orm import PrimaryKey, Required, Optional as PonyOptional, Set, Json

from .config import Entity
from .substance import Substance


class Reaction(Entity):
    id = PrimaryKey(int, auto=True)
    substances = Set('ReactionSubstance')

    def __init__(self, reaction: Optional[ReactionContainer] = None, /):
        super().__init__()
        if reaction is not None:
            if any(x.connected_components_count > 1 for x in reaction.molecules()):
                raise ValueError('Each CGRtools.MoleculeContainer should contain only one sole graph(molecule)')
            for x in reaction.reactants:
                ReactionSubstance(Substance(((x, None),)), self, False)
            for x in reaction.products:
                ReactionSubstance(Substance(((x, None),)), self, True)


class ReactionSubstance(Entity):
    id = PrimaryKey(int, auto=True)
    is_product = Required(bool)
    _mapping = PonyOptional(Json, column='mapping', nullable=True)
    reaction = Required('Reaction')
    substance = Required('Substance')

    def __init__(self, substance: Substance, reaction: Reaction, is_product: bool, /, *,
                 mapping: Union[Dict[int, int], List[Tuple[int, int]], None] = None):
        if isinstance(mapping, dict):
            if not all(isinstance(x, int) for x in mapping.items() for x in x):
                raise TypeError('Mapping is dict that contains following structure: '
                                '{int (atom of substance) : int (atom of Substance)}')
            mapping = list(mapping.items())
        elif isinstance(mapping, (tuple, list)):
            if not all(isinstance(x, (tuple, list)) and len(x) == 2 and all(isinstance(x, int) for x in x)
                       for x in mapping):
                raise TypeError('For the tuple pairs of integers following constructions should be applied: '
                                '[(int (atom of substance), int (atom of Substance)),...]')
        elif mapping is not None:
            raise TypeError('Mapping of Reaction to Substance expected')
        super().__init__(reaction=reaction, substance=substance, is_product=is_product, _mapping=mapping)

    @cached_property
    def structure(self):
        if self.mapping:
            return self.substance.structure.remap(self.mapping, copy=True)
        return self.substance.structure

    @cached_property
    def mapping(self):
        return dict(self._mapping) if self._mapping else {}


__all__ = ['Reaction', 'ReactionSubstance']
