def diffusionConvectionKernel(mesh, a=None, b=0.0,
                              uB=None, duB=None,
                              vel=0,
                              # u0=0,
                              fn=None,
                              scheme='CDS', sparse=False, time=0.0,
                              userData=None):
    """
    Generate system matrix for diffusion and convection in a velocity field.

    Particle concentration u inside a velocity field.

    Peclet Number - ratio between convection/diffusion = F/D
        F = velocity flow trough volume boundary,
        D = diffusion coefficient

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        Mesh represents spatial discretization of the calculation domain
    a: value | array | callable(cell, userData)
        Diffusion coefficient per cell
    b: value | array | callable(cell, userData)
        TODO What is b
    fn: iterable(cell)
        TODO What is fn
    vel: ndarray (N,dim) | RMatrix(N,dim)
        velocity field [[v_i,]_j,] with i=[1..3] for the mesh dimension
        and j = [0 .. N-1] per Cell or per Node so N is either
        mesh.cellCount() or mesh.nodeCount()
    scheme: str [CDS]
        Finite volume scheme
        * CDS -- Central Difference Scheme.
            maybe irregular for Peclet no. |F/D| > 2
            Diffusion dominant. Error of order 2
        * UDS -- Upwind Scheme.
            Convection dominant. Error of order 1
        * HS -- Hybrid Scheme.
            Diffusion dominant for Peclet-number |(F/D)| < 2
            Convection dominant else.
        * PS -- Power Law Scheme.
            Identical to HS for Peclet-number |(F/D)| > 10 and near to ES else
            Convection dominant.
        * ES -- Exponential scheme
            Only stationary one-dimensional but exact solution
    Returns
    -------
    S: :gimliapi:`GIMLI::SparseMatrix` | numpy.ndarray(nCells, nCells)
        Kernel matrix, depends on vel, a, b, scheme, uB, duB .. if some of this
        has been changed you cannot cache these matrix
    rhsBoundaryScales: ndarray(nCells)
        RHS offset vector
    """
    if a is None:
        a = pg.Vector(mesh.boundaryCount(), 1.0)

    AScheme = None
    if scheme == 'CDS':
        AScheme = lambda peclet_: 1.0 - 0.5 * abs(peclet_)
    elif scheme == 'UDS':
        AScheme = lambda peclet_: 1.0
    elif scheme == 'HS':
        AScheme = lambda peclet_: max(0.0, 1.0 - 0.5 * abs(peclet_))
    elif scheme == 'PS':
        AScheme = lambda peclet_: max(0.0, (1.0 - 0.1 * abs(peclet_))**5.0)
    elif scheme == 'ES':
        AScheme = lambda peclet_: (peclet_) / (np.exp(abs(peclet_)) - 1.0) \
            if peclet_ != 0.0 else 1
    else:
        raise BaseException("Scheme unknwon:" + scheme)

    useHalfBoundaries = False

    dof = mesh.cellCount()

    if not uB:
        uB = []
    if not duB:
        duB = []

    if useHalfBoundaries:
        dof = mesh.cellCount() + len(uB)

    S = None
    if sparse:
        S = pg.matrix.SparseMapMatrix(dof, dof, stype=0) + identity(dof, scale=b)
    else:
        S = np.zeros((dof, dof))

    rhsBoundaryScales = np.zeros(dof)

    # we need this to fast identify uBoundary and value by boundary
    uBoundaryID = []
    uBoundaryVals = [None] * mesh.boundaryCount()

    for [b, val] in uB:

        if isinstance(b, pg.core.Boundary):
            uBoundaryID.append(b.id())
            uBoundaryVals[b.id()] = val
        elif isinstance(b, pg.core.Node):
            for _b in b.boundSet():
                if _b.rightCell() is None:
                    pg.warn('Dirichlet for one node considered for the nearest boundary.', _b.id())
                    uBoundaryID.append(_b.id())
                    uBoundaryVals[_b.id()] = val
                    break
        else:
            raise BaseException("Please give boundary, value list")


    duBoundaryID = []
    duBoundaryVals = [None] * mesh.boundaryCount()

    for [boundary, val] in duB:
        if not isinstance(boundary, pg.core.Boundary):
            raise BaseException("Please give boundary, value list")

        duBoundaryID.append(boundary.id())
        duBoundaryVals[boundary.id()] = val

    # iterate over all cells
    for cell in mesh.cells():
        cID = cell.id()
        for bi in range(cell.boundaryCount()):
            boundary = pg.core.findBoundary(cell.boundaryNodes(bi))

            ncell = boundary.leftCell()
            if ncell == cell:
                ncell = boundary.rightCell()

            v = findVelocity(mesh, vel, boundary, cell, ncell)

            # Convection part
            F = boundary.norm(cell).dot(v) * boundary.size()
            # print(F, boundary.size(), v, vel)
            # Diffusion part
            D = findDiffusion(mesh, a, boundary, cell, ncell)

            # print(F, D, F/D)
            # print((1.0 - 0.1 * abs(F/D))**5.0)
            if D > 0:
                aB = D * AScheme(F / D) + max(-F, 0.0)
            else:
                aB = max(-F, 0.0)

            aB /= cell.size()

            # print(cell.center(), boundary.center(), boundary.norm(cell), aB)
            if ncell:
                # no boundary
                if sparse:
                    S.addVal(cID, ncell.id(), -aB)
                    S.addVal(cID, cID, +aB)
                else:
                    S[cID, ncell.id()] -= aB
                    S[cID, cID] += aB

            elif not useHalfBoundaries:

                if boundary.id() in uBoundaryID:
                    val = pg.solver.generateBoundaryValue(boundary,
                                                uBoundaryVals[boundary.id()],
                                                time=time,
                                                userData=userData)

                    if sparse:
                        S.addVal(cID, cID, aB)
                    else:
                        S[cID, cID] += aB

                    #print(aB, uBoundaryVals[boundary.id()])
                    rhsBoundaryScales[cID] += aB * np.mean(val)

                if boundary.id() in duBoundaryID:
                    # Neumann boundary condition
                    val = pg.solver.generateBoundaryValue(
                        boundary,
                        duBoundaryVals[boundary.id()],
                        time=time,
                        userData=userData)

                    val = np.mean(val)

                    # amount of flow through the boundary .. maybe buggy
                    # fill be replaced by suitable FE solver
                    outflow = -val * boundary.size() / cell.size()

                    if sparse:
                        S.addVal(cID, cID, outflow)
                    else:
                        S[cID, cID] += outflow

        if fn is not None:
            if sparse:
                # * cell.shape().domainSize())
                S.addVal(cell.id(), cell.id(), -fn[cell.id()])
            else:
                # * cell.shape().domainSize()
                S[cell.id(), cell.id()] -= fn[cell.id()]

    return S, rhsBoundaryScales



[docs]
def solveFiniteVolume(mesh, a=1.0, b=0.0, f=0.0, fn=0.0, vel=None, u0=0.0,
                      times=None,
                      uL=None, relax=1.0,
                      ws=None, scheme='CDS', **kwargs):
    r"""Solve partial differential equation with Finite Volumes.

    This function is a syntactic sugar proxy for using the Finite Volume
    functionality of the library core to solve elliptic and parabolic partial
    differential of the following type:

    .. math::

        \frac{\partial u}{\partial t} + \mathbf{v}\cdot\nabla u & = \nabla\cdot(a \nabla u) + b u + f(\mathbf{r},t)\\
        u(\mathbf{r}, t) & = u_B  \quad\mathbf{r}\in\Gamma_{\text{Dirichlet}}\\
        \frac{\partial u(\mathbf{r}, t)}{\partial \mathbf{n}} & = u_{\partial \text{B}}  \quad\mathbf{r}\in\Gamma_{\text{Neumann}}\\
        u(\mathbf{r}, t=0) & = u_0 \quad\text{with} \quad\mathbf{r}\in\Omega

    The Domain :math:`\Omega` and the Boundary :math:`\Gamma` are defined
    through the given mesh with appropriate boundary marker.

    The solution :math:`u(\mathbf{r}, t)` is given for each cell in the mesh.

    TODO:

     * Refactor with solver class and Runga-Kutte solver
     * non steady boundary conditions

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        Mesh represents spatial discretization of the calculation domain
    a: value | array | callable(cell, userData)
        Stiffness weighting per cell values.
    b: value | array | callable(cell, userData)
        Scale for mass values b
    f: iterable(cell)
        Load vector
    fn: iterable(cell)
        TODO What is fn
    vel: ndarray (N,dim) | RMatrix(N,dim)
        Velocity field :math:`\mathbf{v}(\mathbf{r}, t=\text{const}) = \{[v_i]_j,\}`
        with :math:`i=[1\ldots 3]` for the mesh dimension
        and :math:`j = [0\ldots N-1]` with N either the amount of cells,
        nodes, or boundaries.
        Velocities per boundary are preferred and will be interpolated
        on demand.
    u0: value | array | callable(cell, userData)
        Starting field
    times: iterable
        Time steps to calculate for.
    ws Workspace
        This can be an empty class that will used as an Workspace to store and
        cache data.

        If ws is given: The system matrix is taken from ws or
        calculated once and stored in ws for further usage.

        The system matrix is cached in this Workspace as ws.S
        The LinearSolver with the factorized matrix is cached in
        this Workspace as ws.solver
        The rhs vector is only stored in this Workspace as ws.rhs
    scheme: str [CDS]
        Finite volume scheme:
        :py:mod:`pygimli.solver.diffusionConvectionKernel`
    **kwargs:
        * bc : Boundary Conditions dictionary, see pg.solver
        * uB : Dirichlet boundary conditions DEPRECATED
        * duB : Neumann boundary conditions DEPRECATED

    Returns
    -------
    u: ndarray(nTimes, nCells)
        Solution field for all time steps.
    """
    verbose = kwargs.pop('verbose', False)
    # The Workspace is to hold temporary data or preserve matrix rebuild
    # swatch = pg.Stopwatch(True)
    sparse = True

    workspace = pg.solver.WorkSpace()
    if ws:
        workspace = ws

    a = pg.solver.parseArgToArray(a, [mesh.cellCount(), mesh.boundaryCount()])
    b = pg.solver.parseArgToArray(b, mesh.cellCount())
    f = pg.solver.parseArgToArray(f, mesh.cellCount())
    fn = pg.solver.parseArgToArray(fn, mesh.cellCount())

    boundsDirichlet = None
    boundsNeumann = None

    # BEGIN check for Courant-Friedrichs-Lewy
    if vel is not None:

        if isinstance(vel, float):
            print("Warning! .. velocity is float and no vector field")

        # we need velocities for boundaries
        if len(vel) is not mesh.boundaryCount():
            if len(vel) == mesh.cellCount():
                vel = pg.meshtools.cellDataToNodeData(mesh, vel)
            elif len(vel) == mesh.nodeCount():
                vel = pg.meshtools.nodeDataToBoundaryData(mesh, vel)
            else:
                print("mesh:", mesh)
                print("vel:", vel.shape)

                raise Exception("Cannot determine data format for velocities")

        if times is not None:
            pg.solver.checkCFL(times, mesh, np.max(pg.abs(vel)))

    if not hasattr(workspace, 'S'):

        boundsDirichlet = []
        if 'bc' in kwargs:
            bct = dict(kwargs['bc'])
            if 'Dirichlet' in bct:
                boundsDirichlet += pg.solver.parseArgToBoundaries(bct.pop('Dirichlet'), mesh)

            if 'Node' in bct:
                n = bct.pop('Node')
                boundsDirichlet.append([mesh.node(n[0]), n[1]])

            if 'Neumann' in bct:
                boundsNeumann = pg.solver.parseArgToBoundaries(bct.pop('Neumann'), mesh)


        if 'uB' in kwargs:
            pg.deprecated('use new bc dictionary')
            boundsDirichlet = pg.solver.parseArgToBoundaries(kwargs['uB'],
                                                             mesh)

        if 'duB' in kwargs:
            pg.deprecated('use new bc dictionary')
            boundsNeumann = pg.solver.parseArgToBoundaries(kwargs['duB'], mesh)

        workspace.S, workspace.rhsBCScales = diffusionConvectionKernel(
            mesh=mesh,
            a=a,
            b=b,
            uB=boundsDirichlet,
            duB=boundsNeumann,
            # u0=u0,
            fn=fn,
            vel=vel,
            scheme=scheme,
            sparse=sparse,
            userData=kwargs.pop('userData', None))

        dof = len(workspace.rhsBCScales)
        workspace.ap = np.zeros(dof)

        # for nonlinears
        if uL is not None:
            for i in range(dof):
                val = 0.0
                if sparse:
                    val = workspace.S.getVal(i, i) / relax
                    workspace.S.setVal(i, i, val)
                else:
                    val = workspace.S[i, i] / relax
                    workspace.S[i, i] = val

                workspace.ap[i] = val

        # print('FVM kernel 2:', swatch.duration(True))
    # endif: not hasattr(workspace, 'S'):

    workspace.rhs = np.zeros(len(workspace.rhsBCScales))
    workspace.rhs[0:mesh.cellCount()] = f  # * mesh.cellSizes()

    workspace.rhs += workspace.rhsBCScales

    # for nonlinear: relax progress with scaled last result
    if uL is not None:
        workspace.rhs += (1. - relax) * workspace.ap * uL
    # print('FVM: Prep:', swatch.duration(True))

    if not hasattr(times, '__len__'):

        if sparse and not hasattr(workspace, 'solver'):
            Sm = pg.matrix.SparseMatrix(workspace.S)
            # hold Sm until we have reference counting,
            # loosing Sm here will kill LinSolver later
            workspace.Sm = Sm
            workspace.solver = pg.core.LinSolver(Sm, verbose=verbose)

        u = None
        if sparse:
            u = workspace.solver.solve(workspace.rhs)
        else:
            u = np.linalg.solve(workspace.S, workspace.rhs)
        # print('FVM solve:', swatch.duration(True))
        return u[0:mesh.cellCount():1]
    else:
        theta = kwargs.pop('theta', 0.5 + 1e-6)

        if sparse:
            I = pg.solver.identity(len(workspace.rhs))
        else:
            I = np.diag(np.ones(len(workspace.rhs)))

        progress = None
        if verbose:
            from pygimli.utils import ProgressBar
            progress = ProgressBar(its=len(times))

            print("Solve timesteps with Crank-Nicolson.")

        return pg.solver.crankNicolson(times, workspace.S, I,
                    f=workspace.rhs,
                    u0=pg.solver.cellValues(mesh, u0),
                    theta=theta,
                    progress=progress)

