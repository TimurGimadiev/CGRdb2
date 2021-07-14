from functools import cached_property
from typing import Optional, Union, List, Dict, Tuple

from CGRtools.containers import ReactionContainer, CGRContainer
from pony.orm import PrimaryKey, Required, Optional as PonyOptional, Set, Json, IntArray
from pickle import loads, dumps
from .config import Entity
from . import substance
from StructureFingerprint import LinearFingerprint
from ..utils import validate_molecule
from . import molecule
from . import config
from .indexes import CursorHolder, RequestPack
from functools import partial


class CGR(Entity):
    id = PrimaryKey(int, auto=True)
    cgrsmiles = Required(str, lazy=True)
    signature = Required(bytes, unique=True, lazy=True)
    fingerprint = Required(IntArray, lazy=True)
    fingerprint_len = Required(int)
    _structure = Required(bytes, lazy=True)
    CGRReaction = set('CGRReaction')

    @cached_property
    def structure(self) -> CGRContainer:
        return loads(self._structure)

    @staticmethod
    def get(cls, mol: Optional[CGRContainer] = None, **kwargs):  # fix this
        if mol is None:
            return type(cls).get(cls, **kwargs)
        elif not isinstance(mol, CGRContainer):
            raise ValueError("CGRtools.CGRContainer should be provided")
        return CGR.get(signature=bytes(mol))

    def similar(self):

        raise NotImplemented

    def substructures(self):
        raise NotImplemented


class CGRReaction(Entity):
    id = PrimaryKey(int, auto=True)
    reactions = set('Reaction')
    cgr_data = Required('CGRData')

    @cached_property
    def structure(self) -> CGRContainer:
        return self.cgr_data.structure

    def __init__(self, structure: ReactionContainer, /, reaction=None):
        if reaction:
            cgr_structure = ~structure
            try:
                cgr_data = CGR.get(cgrsmiles=str(cgr_structure))
            except Exception:
                fingerprint = LinearFingerprint(**config.fingerprint).transform_bitset([cgr_structure])[0]
                cgr_data = CGR(cgrsmiles=str(cgr_structure), fingerprint=fingerprint,
                                   fingerpint_len=len(fingerprint), _structure=dumps(cgr_structure))
            super().__init__(reactions=reaction, cgr_data=cgr_data)

