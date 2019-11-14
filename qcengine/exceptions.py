import traceback


class QCEngineException(Exception):
    """
    Base QCEngine exception, should never be called explicitly.
    """

    error_type = "base_error"
    header = "QCEngine Base Error"

    def __init__(self, message: str):

        # Call the base class constructor with the parameters it needs
        super().__init__(message)

        # Now for your custom code...
        self.raw_message = message
        self.traceback = traceback.format_exc()

    @property
    def error_message(self) -> str:
        return f"{self.header}: {self.raw_message}"


class UnknownError(QCEngineException):
    """
    Unknown QCEngine error, the type was not able to be specified.
    """

    error_type = "unknown_error"
    header = "QCEngine Unknown Error"


class InputError(QCEngineException):
    """
    Incorrect user parameters, not recoverable. Also may indicate a version issue.
    """

    error_type = "input_error"
    header = "QCEngine Input Error"


class ResourceError(QCEngineException):
    """
    Not enough resources for computation such as not enough memory, cores, or disk was not available.
    Recoverable on different compute hardware.
    """

    error_type = "resource_error"
    header = "QCEngine Resource Error"


class ConvergenceError(QCEngineException):
    """
    Failed iteration convergence error, likely recoverable with tweaked parameters.
    """

    error_type = "convergence_error"
    header = "QCEngine Convergence Error"


class RandomError(QCEngineException):
    """
    Likely recoverable errors such as segmentation faults or disk io problems.
    """

    error_type = "random_error"
    header = "QCEngine Random Error"
