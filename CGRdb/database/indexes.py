from .config import db
from pony.orm import Required, PrimaryKey, IntArray

class MoleculeIndex(db.Entity):
    band = Required(int)
    key = Required(int, size=64)
    records = Required(IntArray)
    PrimaryKey(band, key)


#class ReactionIndex(db.Entity):
#    id = PrimaryKey(int, auto=True)
#    reaction = Required('Reaction')
#    signature = Required(bytes, unique=True,  column='signature')
#    fingerprint = Required(IntArray, optimistic=False, index=False, volatile=True, column='fingerprint')
#    structures = Required(IntArray, optimistic=False, index=False,  volatile=True, column='structures')


__all__ = ['MoleculeIndex']