##############################################################################################
# MDTraj extension file for XYZ trajectory file of IPI program
#
##############################################################################################

try:
    import mdtraj as md
    import numpy as np
    MDTRAJ_IMPORTED = 1
except:
    pass


@FormatRegistry.register_loader('.xyzipi')
@FormatRegistry.register_loader('.xyzipi.gz')
def load_xyzipi(filename, top=None, stride=None, atom_indices=None, frame=None):
    """Load a xyz ipi trajectory file.

    load xyz file loaded by us

    Parameters
    ----------
    filename : str
        String filename of xyz trajectory file.
    top : {str, Trajectory, Topology}
        The xyz format does not contain topology information. Pass in
        either the path to a pdb file, a trajectory, or a topology to supply
        this information.
    stride : int, default=None
        Only read every stride-th frame
    atom_indices : array_like, optional
        If not none, then read only a subset of the atoms coordinates from the
        file.
    frame : int, optional
        Use this option to load only a single frame from a trajectory on disk.
        If frame is None, the default, the entire trajectory will be loaded.
        If supplied, ``stride`` will be ignored.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.XYZTrajectoryFile :  Low level interface to xyz files
    """
    from mdtraj.core.trajectory import _parse_topology, Trajectory

    # We make `top` required. Although this is a little weird, its good because
    # this function is usually called by a dispatch from load(), where top comes
    # from **kwargs. So if its not supplied, we want to give the user an
    # informative error message.
    if top is None:
        raise ValueError('"top" argument is required for load_xyz')

    if not isinstance(filename, string_types):
        raise TypeError('filename must be of type string for load_xyz. '
                        'you supplied %s'.format(type(filename)))

    topology = _parse_topology(top)
    atom_indices = cast_indices(atom_indices)

    with XYZTrajectoryFile(filename) as f:
        if frame is not None:
            f.seek(frame)
            n_frames = 1
        else:
            n_frames = None
        return f.read_as_traj(topology, n_frames=n_frames, stride=stride,
                              atom_indices=atom_indices)


class XYZIPITrajectoryFile(md.formats.XYZTrajectoryFile):

    def read_as_traj(self, topology, n_frames=None, stride=None, atom_indices=None):
        """Read a trajectory from a XYZ file

        Parameters
        ----------
        topology : Topology
            The system topology
        n_frames : int, optional
            If positive, then read only the next `n_frames` frames. Otherwise read all
            of the frames in the file.
        stride : np.ndarray, optional
            Read only every stride-th frame.
        atom_indices : array_like, optional
            If not none, then read only a subset of the atoms coordinates from the
            file. This may be slightly slower than the standard read because it required
            an extra copy, but will save memory.

        Returns
        -------
        trajectory : Trajectory
            A trajectory object containing the loaded portion of the file.
        """
        from mdtraj.core.trajectory import Trajectory
        if atom_indices is not None:
            topology = topology.subset(atom_indices)

        initial = int(self._frame_index)
        xyz,cell_lengths = self.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        if len(xyz) == 0:
            return Trajectory(xyz=np.zeros((0, topology.n_atoms, 3)), topology=topology)

        in_units_of(xyz, self.distance_unit, Trajectory._distance_unit, inplace=True)
        in_units_of(cell_lengths, self.distance_unit, Trajectory._distance_unit, inplace=True)

        if cell_lengths is None:
            cell_angles = None
        else:
            # Assume that its a rectilinear box
            cell_angles = 90.0 * np.ones_like(cell_lengths)

        if stride is None:
            stride = 1
        time = (stride*np.arange(len(xyz))) + initial
        t = Trajectory(xyz=xyz, topology=topology, time=time)
        t.unitcell_lengths = cell_lengths
        t.unitcell_angles = cell_angles
        return t
    
    def read(self, n_frames=None, stride=None, atom_indices=None):
        """Read data from a xyz file.

        Parameters
        ----------
        n_frames : int, None
            The number of frames you would like to read from the file.
            If None, all of the remaining frames will be loaded.
        stride : np.ndarray, optional
            Read only every stride-th frame.
        atom_indices : array_like, optional
            If not none, then read only a subset of the atoms coordinates
            from the file.

        Returns
        -------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=np.float32
        cell_lengths : {np.ndarray, None}
        """

        if not self._mode == 'r':
            raise ValueError('read() is only available when file is opened '
                             'in mode="r"')

        if n_frames is None:
            frame_counter = itertools.count()
        else:
            frame_counter = xrange(n_frames)

        if stride is None:
            stride = 1

        all_coords, boxes = [],[] 
        for i in frame_counter:
            try:
                frame_coords, box = self._read()
                if atom_indices is not None:
                    frame_coords = frame_coords[atom_indices, :]
            except _EOF:
                break

            all_coords.append(frame_coords)
            boxes.append(box)

            for j in range(stride - 1):
                # throw away these frames
                try:
                    self._read()
                except _EOF:
                    break

        all_coords = np.array(all_coords)
        if all(b is None for b in boxes):
            # if there was no box information in any frame, that's cool
            return all_coords, None

        if not all(b is not None for b in boxes):
            # but if some of them had box information and others didn't
            # that probably means there was a bug in the parsing.
            raise IOError('Inconsistent box information. Try manually '
                          'setting has_box? Your xyz ipi file might be '
                          'corrupt.')

        return all_coords, np.array(boxes, dtype=np.float32)


    def _read(self):
        """Read a single frame. """

        first = self._fh.readline()  # Number of atoms.
        if first == '':
            raise _EOF()
        else:
            self._n_atoms = int(first)
        
        cmt = self._fh.readline()  # Comment line.
        split_cmt = cmt.split()
        box = [float(x) for x in split_line[3:9]]

        self._line_counter += 2

        xyz = np.empty(shape=(self._n_atoms, 3))
        types = np.empty(shape=self._n_atoms, dtype=str)

        for i in xrange(self._n_atoms):
            line = self._fh.readline()
            if line == '':
                raise _EOF()
            split_line = line.split()
            try:
                types[i] = split_line[0]
                xyz[i] = [float(x) for x in split_line[1:4]]
            except Exception:
                raise IOError('xyz parse error on line {0:d} of "{1:s}". '
                              'This file does not appear to be a valid '
                              'xyz file.'.format(
                        self._line_counter,  self._filename))
            self._line_counter += 1
        # --- end body ---

        self._frame_index += 1
        return xyz,box