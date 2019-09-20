"""
Integrates the computes together
"""
from typing import Any, Dict, Optional, Union

import qcelemental as qcel

from .compute import compute

try:
    from mdi import MDI_Init, MDI_Accept_Communicator, MDI_Recv_Command
    from mdi import MDI_Recv, MDI_Send
    from mdi import MDI_DOUBLE, MDI_INT
    from mdi import MDI_Get_Intra_Code_MPI_Comm
    use_mdi = True
except ImportError:
    use_mdi = False

try:
    from mpi4py import MPI
    use_mpi4py = True
except ImportError:
    use_mpi4py = False


class MDIEngine():
    def __init__(self, input_data: Union[Dict[str, Any], 'ResultInput'],
            program: str,
            raise_error: bool = False,
            local_options: Optional[Dict[str, Any]] = None):
        """ Initialize an MDIEngine object for communication with MDI

        Parameters
        ----------
        input_data : Union[Dict[str, Any], 'ResultInput']
            A QCSchema input specification in dictionary or model from QCElemental.models
        program : str
            The program to execute the input with.
        raise_error : bool, optional
            Determines if compute should raise an error or not.
        local_options : Optional[Dict[str, Any]], optional
            A dictionary of local configuration options
        """

        # Input variables
        self.input_data = input_data
        self.program = program
        self.raise_error = raise_error
        self.local_options = local_options

        # Output of most recent compute call
        self.compute_return = None

        # MPI variables
        self.mpi_world = None
        self.world_rank = 0

        # Get correct intra-code MPI communicator
        if use_mpi4py:
            self.mpi_world = MDI_Get_Intra_Code_MPI_Comm()
            self.world_rank = self.mpi_world.Get_rank()

            # QCEngine does not currently support multiple MPI ranks
            if self.mpi_world.Get_size() != 1:
                MPI.COMM_WORLD.Abort()

        # Flag to stop listening for MDI commands
        self.stop_listening = False

        print( "========================================================" )
        print( self.input_data )
        print( "========================================================" )

        # Dictionary of all supported MDI commands
        self.commands = {
            "<ENERGY": self.send_energy,
            "<FORCES": self.send_forces,
            ">COORDS": self.recv_coords,
            "SCF": self.run_scf,
            "EXIT": self.exit
        }

        # Accept a communicator to the driver code
        self.comm = MDI_Accept_Communicator()

    # Respond to the <ENERGY command
    def send_energy(self):
        """ Send the total energy through MDI

        :returns: *energy* Energy of the system
        """
        self.run_scf()
        energy = self.compute_return.return_result
        MDI_Send(energy, 1, MDI_DOUBLE, self.comm)
        return energy

    # Respond to the <FORCES command
    def send_forces(self):
        """ Send the nuclear forces through MDI

        :returns: *forces* Atomic forces
        """
        input = qcel.models.ResultInput(
            molecule = self.input_data.molecule, 
            driver = "gradient", 
            model = self.input_data.model, 
            keywords = self.input_data.keywords
            )
        self.compute_return = compute(input, self.program, self.raise_error, self.local_options)

        forces = self.compute_return.return_result
        MDI_Send(forces, len(forces), MDI_DOUBLE, self.comm)
        return forces

    # Respond to the SCF command
    def run_scf(self):
        """ Run an energy calculation
        """
        self.compute_return = compute(self.input_data, self.program, self.raise_error, self.local_options)

    # Respond to the >COORDS command
    def recv_coords(self, coords=None):
        """ Receive a set of nuclear coordinates through MDI and assign them to the atoms in the current molecule

        Parameters
        ----------
        coords: New nuclear coordinates. If None, receive through MDI.
        """
        natom = len(self.input_data.molecule.geometry)
        if coords is None:
            coords = MDI_Recv(3 * natom, MDI_DOUBLE, self.comm)
        for iatom in range(natom):
            for icoord in range(3):
                self.input_data.molecule.geometry[iatom][icoord] = coords[3*iatom + icoord]

    # Respond to the EXIT command 
    def exit(self):
        """ Stop listening for MDI commands
        """
        self.stop_listening = True

    # Enter server mode, listening for commands from the driver
    def listen_for_commands(self):
        """ Receive commands through MDI and respond to them as defined by the MDI Standard
        """

        while not self.stop_listening:
            if self.world_rank == 0:
                command = MDI_Recv_Command(self.comm)
            else:
                command = None
            if use_mpi4py:
                command = self.mpi_world.bcast(command, root=0)
            if self.world_rank == 0:
                print("MDI command received: " + str(command))

            # Search for this command in self.commands
            found_command = False
            for supported_command in self.commands:
                if not found_command and command == supported_command:
                    # Run the function corresponding to this command
                    self.commands[supported_command]()
                    found_command = True
            if not found_command:
                raise Exception('Unrecognized command: ' + str(command))


def mdi(input_data: Union[Dict[str, Any], 'ResultInput'],
            program: str,
            raise_error: bool = False,
            local_options: Optional[Dict[str, Any]] = None):
    """Executes a single quantum chemistry program given a QC Schema input.

    The full specification can be found at:
        http://molssi-qc-schema.readthedocs.io/en/latest/index.html#

    Parameters
    ----------
    input_data : Union[Dict[str, Any], 'ResultInput']
        A QCSchema input specification in dictionary or model from QCElemental.models
    program : str
        The program to execute the input with.
    raise_error : bool, optional
        Determines if compute should raise an error or not.
    local_options : Optional[Dict[str, Any]], optional
        A dictionary of local configuration options
    """

    if not use_mdi:
        raise Exception("Unable to find MDI package")

    engine = MDIEngine(input_data, program, raise_error, local_options)
    engine.listen_for_commands()

    return



def mdi_init(mdi_arguments):
    """ Initialize the MDI Library

    Arguments:
        mdi_arguments: MDI configuration options
    """

    mpi_world = None
    if use_mpi4py:
        mpi_world = MPI.COMM_WORLD
    MDI_Init(mdi_arguments, mpi_world)
