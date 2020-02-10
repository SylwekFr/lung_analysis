from aenum import Enum

class SeriesType(Enum):

    _init_ = 'value string'

    BEFORE = 1, 'Before'
    AFTER = 2, 'After'

    def __str__(self):
        return self.string

    @classmethod
    def _missing_value_(cls, value):
        for member in cls:
            if member.string == value:
                return member