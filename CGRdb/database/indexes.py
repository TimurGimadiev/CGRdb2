from .config import db
from pony.orm import Required, PrimaryKey, IntArray

class MoleculeIndex(db.Entity):
    band = Required(int)
    key = Required(int, size=64)
    records = Required(IntArray)
    PrimaryKey(band, key)


__all__ = ['MoleculeIndex']