import abc

from typing import Optional, List

from pydantic import BaseModel

class ProgramExecutor(BaseModel, abc.ABC):

    name: str
    requires_folder: bool
    requires_scratch: bool
    single_node: bool
    thread_safe: bool
    max_cores: Optional[int]
    max_memory: Optional[float]

    requires_gpu: bool=False
    has_gpu: bool=False

    class Config:
        allow_mutation: False
        extra: "forbid"

    @abc.abstractmethod
    def compute(self, input_data: 'ResultInput', config: 'Config') -> 'Result':
        pass

    def build_input(self, input: 'ResultInput', config: 'Config', template: Optional[str]=None):
        raise ValueError("build_input is not implemented for {}.", self.__class__)

    def execute(self, inputs, extra_outfiles, extra_commands, scratch_name, timeout):
        raise ValueError("execute is not implemented for {}.", self.__class__)

    def parse_output(self):
        raise ValueError("parse_output is not implemented for {}.", self.__class__)


