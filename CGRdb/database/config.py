from pony.orm import PrimaryKey, Required, Json, Database, db_session


db = Database()
Entity = db.Entity

config = {}


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


__all__ = ['Config', 'Entity', 'db']
