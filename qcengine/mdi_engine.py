"""Support for using QCEngine as an MDI engine.
For details regarding MDI, see https://molssi.github.io/MDI_Library/html/index.html.

"""
from typing import Any, Dict, Optional, Union

import qcelemental as qcel

from .compute import compute

try:
    from mdi import MDI_Init, MDI_Accept_Communicator, MDI_Recv_Command
    from mdi import MDI_Recv, MDI_Send
    from mdi import MDI_DOUBLE, MDI_CHAR, MDI_INT, MDI_COMMAND_LENGTH
    from mdi import MDI_Get_Intra_Code_MPI_Comm
    use_mdi = True
except ImportError:
    use_mdi = False

try:
    from mpi4py import MPI
    use_mpi4py = True
except ImportError:
    use_mpi4py = False


class MDIServer():
    def __init__(self, mdi_options: str, 
            program: str,
            molecule,
            model,
            keywords,
            raise_error: bool = False,
            local_options: Optional[Dict[str, Any]] = None):
        """ Initialize an MDIServer object for communication with MDI

        Parameters
        ----------
        mdi_options: str
            Options used during MDI initialization.
        program : str
            The program to execute the input with.
        molecule
            The initial state of the molecule.
        model
            The simulation model to use.
        keywords
            Program-specific keywords.
        raise_error : bool, optional
            Determines if compute should raise an error or not.
        local_options : Optional[Dict[str, Any]], optional
            A dictionary of local configuration options
        """

        # Confirm that the MDI library has been located
        if not use_mdi:
            raise Exception("Unable to find MDI package")

        # Initialize MDI
        mpi_world = None
        if use_mpi4py:
            mpi_world = MPI.COMM_WORLD
        MDI_Init(mdi_options, mpi_world)

        # Input variables
        self.molecule = molecule
        self.model = model
        self.keywords = keywords
        self.program = program
        self.raise_error = raise_error
        self.local_options = local_options

        # The MDI interface does not currently support multiple fragments
        if len(self.molecule.fragments) != 1:
            raise Exception('The MDI interface does not support multiple fragments')

        # Molecule charge and multiplicity
        self.total_charge = self.molecule.molecular_charge
        self.multiplicity = self.molecule.molecular_multiplicity

        # Flag to track whether the latest molecule specification has been validated
        self.molecule_validated = True

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

        # Dictionary of all supported MDI commands
        self.commands = {
            "<NATOMS": self.send_natoms,
            "<COORDS": self.send_coords,
            "<ENERGY": self.send_energy,
            "<FORCES": self.send_forces,
            ">COORDS": self.recv_coords,
            "SCF": self.run_scf,
            "<NCOMMANDS": self.send_ncommands,
            "<COMMANDS": self.send_commands,
            "<ELEMENTS": self.send_elements,
            ">ELEMENTS": self.recv_elements,
            "<MASSES": self.send_masses,
            "<TOTCHARGE": self.send_total_charge,
            ">TOTCHARGE": self.recv_total_charge,
            "<ELEC_MULT": self.send_multiplicity,
            ">ELEC_MULT": self.recv_multiplicity,
            "EXIT": self.exit
        }

        # Accept a communicator to the driver code
        self.comm = MDI_Accept_Communicator()

    def update_molecule(self, key, value):
        """ Update the molecule

        Parameters
        ----------
        key: Key of the molecular element to update
        value: Update value
        """
        if key == "molecular_charge" or key == "molecular_multiplicity":
            # In order to validate correctly, the charges and multiplicities must be set simultaneously
            try:
                self.molecule = qcel.models.Molecule(**{**self.molecule.dict(),
                                      **{"molecular_charge": self.total_charge},
                                      **{"fragment_charges": [self.total_charge]},
                                      **{"molecular_multiplicity": self.multiplicity},
                                      **{"fragment_multiplicities": [self.multiplicity]}})
                self.molecule_validated = True
            except qcel.exceptions.ValidationError:
                # The molecule didn't validate correctly, but a future >TOTCHARGE or >ELEC_MULT command might fix it
                self.molecule_validated = False
        else:
            try:
                self.molecule = qcel.models.Molecule(**{**self.molecule.dict(), **{key: value}})
                self.molecule_validated = True
            except qcel.exceptions.ValidationError:
                if self.molecule_validated:
                    # This update caused the validation error
                    raise Exception('MDI command caused a validation error')

                    

    # Respond to the <NATOMS command
    def send_natoms(self):
        """ Send the number of atoms through MDI

        :returns: *natom* Number of atoms
        """
        natom = len(self.molecule.geometry)
        MDI_Send(natom, 1, MDI_INT, self.comm)
        return natom

    # Respond to the <COORDS command
    def send_coords(self, coords=None):
        """ Send the nuclear coordinates through MDI

        :returns: *coords* Nuclear coordinates
        """
        natom = len(self.molecule.geometry)

        coords = [ 0.0 for i in range(3 * natom) ]
        for iatom in range(natom):
            for icoord in range(3):
                coords[3 * iatom + icoord] = self.molecule.geometry[iatom][icoord]

        return coords

    # Respond to the >COORDS command
    def recv_coords(self, coords=None):
        """ Receive a set of nuclear coordinates through MDI and assign them to the atoms in the current molecule

        Parameters
        ----------
        coords: New nuclear coordinates. If None, receive through MDI.
        """
        natom = len(self.molecule.geometry)
        if coords is None:
            coords = MDI_Recv(3 * natom, MDI_DOUBLE, self.comm)
        for iatom in range(natom):
            for icoord in range(3):
                self.molecule.geometry[iatom][icoord] = coords[3*iatom + icoord]

    # Respond to the <ENERGY command
    def send_energy(self):
        """ Send the total energy through MDI

        :returns: *energy* Energy of the system
        """
        # Ensure that the molecule currently passes validation
        if not self.molecule_validated:
            raise Exception('MDI attempting to compute energy on an unvalidated molecule')
        self.run_scf()
        energy = self.compute_return.return_result
        MDI_Send(energy, 1, MDI_DOUBLE, self.comm)
        return energy

    # Respond to the <FORCES command
    def send_forces(self):
        """ Send the nuclear forces through MDI

        :returns: *forces* Atomic forces
        """
        # Ensure that the molecule currently passes validation
        if not self.molecule_validated:
            raise Exception('MDI attempting to compute gradients on an unvalidated molecule')

        input = qcel.models.ResultInput(
            molecule = self.molecule, 
            driver = "gradient", 
            model = self.model, 
            keywords = self.keywords
            )
        self.compute_return = compute(input, self.program, self.raise_error, self.local_options)

        forces = self.compute_return.return_result
        MDI_Send(forces, len(forces), MDI_DOUBLE, self.comm)
        return forces

    # Respond to the SCF command
    def run_scf(self):
        """ Run an energy calculation
        """
        input = qcel.models.ResultInput(
            molecule = self.molecule, 
            driver = "energy", 
            model = self.model, 
            keywords = self.keywords
            )
        self.compute_return = compute(input, self.program, self.raise_error, self.local_options)

    # Respond to the <NCOMMANDS command
    def send_ncommands(self):
        """ Send the number of supported MDI commands through MDI

        :returns: *ncommands* Number of supported commands
        """
        ncommands = len(self.commands)
        MDI_Send(ncommands, 1, MDI_INT, self.comm)
        return ncommands

    # Respond to the <COMMANDS command
    def send_commands(self):
        """ Send the supported MDI commands through MDI

        :returns: *command_string* String containing the name of each supported command
        """
        command_string = "".join([f"{c:{MDI_COMMAND_LENGTH}}" for c in self.commands.keys()])

        # confirm that command_string is the correct length
        if len(command_string) != len(self.commands) * MDI_COMMAND_LENGTH:
            raise Exception('Programming error: MDI command_string is incorrect length')

        MDI_Send(command_string, len(command_string), MDI_CHAR, self.comm)
        return command_string

    # Respond to the <ELEMENTS command
    def send_elements(self):
        """ Send the atomic number of each nucleus through MDI

        :returns: *elements* Element of each atom
        """
        natom = len(self.molecule.geometry)
        elements = [ qcel.periodictable.to_atomic_number(self.molecule.symbols[iatom]) for iatom in range(natom) ]
        MDI_Send(elements, natom, MDI_INT, self.comm)
        return elements

    # Respond to the >ELEMENTS command
    def recv_elements(self):
        """ Receive a set of atomic numbers through MDI and assign them to the atoms in the current molecule

        Arguments:
            elements: New element numbers. If None, receive through MDI.
        """
        natom = len(self.molecule.geometry)
        if elements is None:
            elements = MDI_Recv(natom, MDI_DOUBLE, self.comm)

        for iatom in range(natom):
            self.molecule.symbols[iatom] = qcel.to_symbol(elements[iatom])

        return elements

    # Respond to the <MASSES command
    def send_masses(self):
        """ Send the nuclear masses through MDI

        :returns: *masses* Atomic masses
        """
        natom = len(self.molecule.geometry)
        masses = self.molecule.masses
        MDI_Send(masses, natom, MDI_DOUBLE, self.comm)
        return masses

    # Respond to the >MASSES command
    def recv_masses(self, masses=None):
        """ Receive a set of nuclear masses through MDI and assign them to the atoms in the current molecule

        Arguments:
            masses: New nuclear masses. If None, receive through MDI.
        """
        natom = len(self.molecule.geometry)
        if masses is None:
            masses = MDI_Recv(natom, MDI_DOUBLE, self.comm)
        self.update_molecule("masses", masses)

    # Respond to the <TOTCHARGE command
    def send_total_charge(self):
        """ Send the total system charge through MDI

        :returns: *charge* Total charge of the system
        """
        charge = self.molecule.molecular_charge
        MDI_Send(charge, 1, MDI_DOUBLE, self.comm)
        return charge

    # Respond to the >TOTCHARGE command
    def recv_total_charge(self, charge=None):
        """ Receive the total system charge through MDI

        Arguments:
            charge: New charge of the system. If None, receive through MDI.
        """
        if charge is None:
            charge = MDI_Recv(1, MDI_DOUBLE, self.comm)
        self.total_charge = charge

        # Allow a validation error here, because a future >ELEC_MULT command might resolve it
        try:
            self.update_molecule("molecular_charge", self.total_charge)
        except qcel.exceptions.ValidationError:
            pass

    # Respond to the <ELEC_MULT command
    def send_multiplicity(self):
        """ Send the electronic multiplicity through MDI

        :returns: *multiplicity* Multiplicity of the system
        """
        multiplicity = self.molecule.molecular_multiplicity
        MDI_Send(multiplicity, 1, MDI_INT, self.comm)
        return multiplicity

    # Respond to the >ELEC_MULT command
    def recv_multiplicity(self, multiplicity=None):
        """ Receive the electronic multiplicity through MDI

        Arguments:
            multiplicity: New multiplicity of the system. If None, receive through MDI.
        """
        if multiplicity is None:
            multiplicity = MDI_Recv(1, MDI_INT, self.comm)
        self.multiplicity = multiplicity

        # Allow a validation error here, because a future >TOTCHARGE command might resolve it
        try:
            self.update_molecule("molecular_multiplicity", self.multiplicity)
        except qcel.exceptions.ValidationError:
            pass

    # Respond to the EXIT command 
    def exit(self):
        """ Stop listening for MDI commands
        """
        self.stop_listening = True

    # Enter server mode, listening for commands from the driver
    def start(self):
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
