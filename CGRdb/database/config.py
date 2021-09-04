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
from pony.orm import PrimaryKey, Required, Json, Database, db_session
from types import MethodType

db = Database()
Entity = db.Entity
config = {}


def load_schema(self=None, provider='postgres', user='postgres', host='localhost', password="example", database='test',
                port=5444):
    if not self:
        self = db
    self.bind(provider=provider, user=user, host=host, password=password, database=database,
              port=port)
    self.generate_mapping(create_tables=True)
    return self


class Config(Entity):
    key = PrimaryKey(str)
    value = Required(Json, index=False, optimistic=False)


def __getattr__(key):
    if key in ('__path__',):
        raise AttributeError
    try:
        return config[key]
    except KeyError:
        with db_session:
            value = Config[key].value
        config[key] = value
        return value


db.load_schema = MethodType(load_schema, db)

__all__ = ['Config', 'Entity', 'db', 'load_schema']
