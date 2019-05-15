from qcelemental.models import ComputeError, FailedOperation


class UnknownError(Exception):
    """
    Unknown QCEngine error, the type was not able to be specified.
    """
    error_type = "unknown_error"

    def __init__(self, message, input_data, return_data):

        # Call the base class constructor with the parameters it needs
        super().__init__(message)

        # Now for your custom code...
        self.input_data = input_data
        self.return_data = return_data

    def format_failed_operation():
        rdata = self.return_data.copy()
        rdata["error"] = ComputeError(error_type=self.error_type, error_message=self.message)

        return FailedOperation(input_data=self.input_data, **rdata)


class InputError(QCEngineError):
    """
    Incorrect user parameters, not recoverable. Also may indicate a version issue.
    """
    error_type = "input_error"


class ResourceError(QCEngineError):
    """
    Not enough resources for computation such as not enough memory, cores, or disk was not available.
    Recoverable on different compute hardware.
    """
    error_type = "resource_error"


class ConvergenceError(QCEngineError):
    """
    Failed iteration convergence error, likely recoverable with tweak parameters.
    """
    error_type = "convergence_error"


class RandomError(QCEngineError):
    """
    Likely recoverable errors such as segmentation faults or disk io problems.
    """
    error_type = "random_error"
