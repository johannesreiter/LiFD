""" Abstract class database from which all other databases inherit """
from abc import ABCMeta, abstractmethod

__author__ = 'Johannes Reiter'
__date__ = 'Jan 17, 2019'


class Database(metaclass=ABCMeta):

    # which columns are required to evaluate the database
    REQUIRED_COLS = list()

    def __init__(self, name, db_source, threshold=True, weight=1.0):
        """
        Constructor
        :param name: name of database that will also be the column in the extended pandas dataframe that LiFD generates
        :param db_source: path to database file
        :param threshold: threshold such that database would predict functionality
        :param weight: support weight when database meets functionality prediction threshold
        """

        self.name = name
        # path to source of database
        self.db_source = db_source

        # pandas dataframe with database content
        self.db_df = None

        # threshold such that database would predict functionality
        self.threshold = threshold
        # support weight when database meets functionality prediction threshold
        self.weight = weight

        # column in dataframe where database presence is stored
        self.pred_col = 'In_'+name

    def __repr__(self):
        return f'{self.__class__.__name__}({self.db_source!r})'

    @abstractmethod
    def in_database(self, nt_var_key=None, pt_var_key=None):
        return False
