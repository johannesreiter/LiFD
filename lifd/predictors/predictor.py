
from abc import ABCMeta, abstractmethod


class Predictor(metaclass=ABCMeta):

    # which columns are required to run the predictor
    REQUIRED_COLS = []

    # which columns are required to evaluate the prediction
    PREDICTION_COLS = []

    # which reference genomes are supported
    SUPP_REF_GENOMES = []

    @staticmethod
    @abstractmethod
    def get_annotation(var_df, ds_name, input_fp=None, output_prefix='', output_dir=None, output_fp=None,
                       filter_condition=None, reference_genome='hg19'):
        """
        Abstract method that has to be implemented by all subclasses
        :param var_df:
        :param ds_name: name of dataset, possibly related to given filter;
                        used for default naming of input and output files
        :param input_fp:
        :param output_prefix:
        :param output_dir:
        :param output_fp:
        :param filter_condition: select only a subset of given dataframe
        :param reference_genome: human reference genome, e.g. hg19 or hg20
        :return:
        """
        pass

    @staticmethod
    @abstractmethod
    def predict_functionality(row):
        """
        Predict functionality
        :param row:
        :return:
        """
        raise NotImplementedError('Predictor did not implement abstract function to predict functionality!')
