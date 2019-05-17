
import traceback

from qcelemental.models import ComputeError, FailedOperation

class QCEngineException(Exception):
    """
    Base QCEngine exception, should never be called explicitly.
    """
    error_type = "base_error"
    header = "QCEngine Base Error:"

    def __init__(self, message, input_data, return_data):

        # Call the base class constructor with the parameters it needs
        super().__init__(message)

        # Now for your custom code...
        self.input_data = input_data
        self.return_data = return_data
        self.traceback = traceback.format_exc()

    def format_failed_operation(self):
        rdata = self.return_data.copy()

        message = f"{self.header}\n{str(self)}"
        rdata["error"] = ComputeError(error_type=self.error_type, error_message=self.message)

        return FailedOperation(input_data=self.input_data, **rdata)

class UnknownError(QCEngineException):
    """
    Unknown QCEngine error, the type was not able to be specified.
    """
    error_type = "unknown_error"
    header = "QCEngine Unknown Error:"



class InputError(QCEngineException):
    """
    Incorrect user parameters, not recoverable. Also may indicate a version issue.
    """
    error_type = "input_error"
    header = "QCEngine Input Error:"


class ResourceError(QCEngineException):
    """
    Not enough resources for computation such as not enough memory, cores, or disk was not available.
    Recoverable on different compute hardware.
    """
    error_type = "resource_error"
    header = "QCEngine Resource Error:"


class ConvergenceError(QCEngineException):
    """
    Failed iteration convergence error, likely recoverable with tweak parameters.
    """
    error_type = "convergence_error"
    header = "QCEngine Convergence Error:"


class RandomError(QCEngineException):
    """
    Likely recoverable errors such as segmentation faults or disk io problems.
    """
    error_type = "random_error"
    header = "QCEngine Random Error:"
