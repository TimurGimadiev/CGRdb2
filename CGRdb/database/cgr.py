# -*- coding: utf-8 -*-
#
#  Copyright 2021 Timur Gimadiev <timur.gimadiev@gmail.com>
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of CGRdb.
#
#  CGRdb is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from functools import cached_property
from typing import Optional, Union
from CGRtools.containers import ReactionContainer, CGRContainer
from pony.orm import PrimaryKey, Required, IntArray
from pickle import loads, dumps
from .config import Entity
from StructureFingerprint import LinearFingerprint
from . import config
from .indexes import CursorHolder, RequestPack
from functools import partial
import json
from datasketch import MinHash


class CGR(Entity):
    id = PrimaryKey(int, auto=True)
    cgrsmiles = Required(str, lazy=True)
    fingerprint = Required(IntArray, lazy=True)
    fingerprint_len = Required(int)
    _structure = Required(bytes, lazy=True)
    reaction = Required('Reaction')

    def __init__(self, structure: ReactionContainer, /, **kwargs):
        cgr_structure = ~structure
        fingerprint = LinearFingerprint(**config.fingerprint).transform_bitset([cgr_structure])[0]
        super().__init__(cgrsmiles=str(cgr_structure), fingerprint=fingerprint,
                         fingerprint_len=len(fingerprint), _structure=dumps(cgr_structure), **kwargs)

    @cached_property
    def structure(self) -> CGRContainer:
        return loads(self._structure)

    def unload(self):
        self._vals_.pop(CGR._structure, None)
        self._vals_.pop(CGR.fingerprint, None)
        self._vals_.pop(CGR.cgrsmiles, None)
        self.__dict__.pop('structure', None)

    @staticmethod
    def __postprocess_exact_match(result):
        return CGR[result]

    @classmethod
    def get(cls, cgr: Optional[CGRContainer] = None, **kwargs):  # fix this
        if cgr is None:
            return type(cls).get(cls, **kwargs)
        elif not isinstance(cgr, CGRContainer):
            raise ValueError("CGRtools.CGRContainer should be provided")
        request = f"""
                    SELECT x.id
                    FROM cgr x
                    WHERE x.cgrsmiles = '{str(cgr)}'"""
        request_pack = RequestPack(request, cls.__postprocess_exact_match,
                           prefetch_map=(None, None, None))
        return CursorHolder(request_pack)

    @staticmethod
    def __postprocess_unordered_cgrs(result, *, fingerprint=None, substr=None):
        cgr = CGR[result[0]]
        fps = set(cgr.fingerprint)
        tanimoto = len(fps.intersection(fingerprint)) / len(fps.union(fingerprint))
        del cgr._vals_[CGR.fingerprint]

        if substr is None:
            return cgr, tanimoto
        else:
            if substr <= cgr.structure:
                return cgr, tanimoto
            # clean cache for memory keeping
            del cgr._vals_[CGR._structure]
            del cgr.__dict__['structure']

    @staticmethod
    def __postprocess_ordered_cgrs(result, substr=None):
        cgr = CGR[result[0]]
        tanimoto = result[1]
        if substr is None:
            return cgr, tanimoto
        else:
            if substr <= cgr.structure:
                return cgr, tanimoto
            # clean cache for memory keeping
            del cgr._vals_[CGR._structure]
            del cgr.__dict__['structure']

    @classmethod
    def __similarity_filtering(cls, fingerprint: list) -> str:
        h = MinHash(num_perm=config.cgr_lsh_num_permute, hashfunc=hash)
        hashranges = json.loads(config.cgr_hashranges)
        h.update_batch(fingerprint)
        request_end = "\nOR\n".join(f"""band={n} AND key='{hash(int.from_bytes(h.hashvalues[i:j].byteswap().data,
                                                                               'big'))}'""" for n, (i, j) in
                                    enumerate(hashranges, 1))
        return request_end

    @classmethod
    def _similarity_unorderd(cls, fingerprint: list) -> RequestPack:
        request_end = cls.__similarity_filtering(fingerprint)
        request = f"""
            SELECT x.id
            FROM cgr x
            WHERE x.id IN (
            SELECT distinct unnest(records) m FROM cgrsimilarityindex WHERE {request_end})"""
        return RequestPack(request, partial(cls.__postprocess_unordered_cgrs, fingerprint=fingerprint),
                           prefetch_map=(CGR, 0, [CGR.fingerprint]))

    @classmethod
    def _similarity_ordered(cls, fingerprint: list) -> RequestPack:
        request_end = cls.__similarity_filtering(fingerprint)
        request = f'''
           SELECT x.id, icount(x.fingerprint & '{set(fingerprint)}') / icount(x.fingerprint |
                '{set(fingerprint)}')::float t
           FROM cgr x
           WHERE x.id IN (
           SELECT distinct unnest(records) m FROM cgrsimilarityindex WHERE {request_end})
           ORDER BY t DESC
           '''
        return RequestPack(request, cls.__postprocess_ordered_cgrs,
                           prefetch_map=(CGR, 0, None))

    def similars(self: Union[CGRContainer, "CGR"], ordered=True, request_only=False):
        if isinstance(self, CGR):
            fingerprint = self.fingerprint
        elif isinstance(self, CGRContainer):
            fingerprint = LinearFingerprint(**config.fingerprint).transform_bitset([self])[0]
            self = CGR
        else:
            raise ValueError(" Only CGRtools.CGRContainer or CGRdb.CGR")
        if ordered:
            request_pack = self._similarity_ordered(fingerprint)
        else:
            request_pack = self._similarity_unorderd(fingerprint)
        if request_only:
            return request_pack.request
        return CursorHolder(request_pack)
# substructure search block

    @classmethod
    def _substructure_ordered(cls, fingerprint: list, substr=None) -> RequestPack:
        len_fp = len(fingerprint)
        #x = set(fingerprint)
        #request = select((c.id, len_fp / c.fingerprint_len) for c in CGR if
        #   between(c.fingerprint_len, len_fp, len_fp//10+5) and raw_sql(f"c.fingerprint @> {x}"))
        #request = request.get_sql()
        request = f'''
        SELECT x.id, {len_fp}/x.fingerprint_len::float t
        FROM cgr x
        WHERE x.fingerprint_len BETWEEN {len_fp} AND 
        {int(len_fp // 0.1)+5} AND x.fingerprint @> '{set(fingerprint)}'
        '''
        return RequestPack(request, partial(cls.__postprocess_ordered_cgrs, substr=substr),
                           prefetch_map=(CGR, 0, [CGR._structure]))

    @classmethod
    def _substructure_unordered(cls, fingerprint: list, substr=None):
        request = f'''
            SELECT x.id
            FROM cgr x
            WHERE x.fingerprint @> '{set(fingerprint)}'
            '''
        return RequestPack(request, partial(cls.__postprocess_unordered_cgrs,
                                            fingerprint=fingerprint, substr=substr),
                           prefetch_map=(
                           CGR, 0, [CGR.fingerprint, CGR._structure]))

    def substructures(self: Union[CGRContainer, "CGR"], ordered=True, request_only=False):
        if isinstance(self, CGR):
            fingerprint = self.fingerprint
            mol = self.structure
        elif isinstance(self, CGRContainer):
            fingerprint = LinearFingerprint(**config.fingerprint).transform_bitset([self])[0]
            mol = self
            self = CGR
        else:
            raise ValueError(" Only CGRtools.MoleculeContainer or CGRdb.Molecule")
        if ordered:
            request_pack = self._substructure_ordered(fingerprint, substr=mol)
        else:
            request_pack = self._substructure_unordered(fingerprint, substr=mol)
        if request_only:
            return request_pack.request
        return CursorHolder(request_pack)


__all__ = ['CGR']
