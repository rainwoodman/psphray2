import numpy
import os.path
from convert import readblock
from trilinearinterp import trilinearinterp

args = None

GHEADER = numpy.dtype([
      ('N', ('u4', 6)),
      ('mass', ('f8', 6)),
      ('time', 'f8'),
      ('redshift', 'f8'),
      ('flag_sfr', 'i4'),     # unused 
      ('flag_feedback', 'i4'),  # unused
      ('Ntot_low', ('u4', 6)),
      ('flag_cool', 'i4'),  # unused
      ('Nfiles', 'i4'),
      ('boxsize', 'f8'),
      ('OmegaM', 'f8'),
      ('OmegaL', 'f8'),
      ('h', 'f8'),
      ('flag_sft', 'i4'),
      ('flag_met', 'i4'),
      ('Ntot_high', ('u4', 6)),
      ('flag_entropy', 'i4'),
      ('flag_double', 'i4'),
      ('flag_ic_info', 'i4'),
      ('flag_lpt_scalingfactor', 'i4'),
      ('unused', ('u4', 12)),
  ])

def gadget(_args):
  global args 
  args = _args
  this = setup_level(0)
  FID = 0
  IDstart = 0
  H = []
  while this is not None:
    if this.ilevel < len(args.Levels) - 1:
      next = setup_level(this.ilevel + 1)
      dig(this, next)
    else: 
      next = None
    bigdispcount = (this.data['disp'] > args.BoxSize / this.Nmesh).sum()
    if bigdispcount > 0:
      print bigdispcount, ' > mesh division', args.BoxSize / this.Nmesh, 'expecting ugly effects'
    HH, IDstart = write_gadget(this, FID, IDstart)
    FID += len(HH)
    H += HH
    this = next

  Ntot = numpy.zeros(6, dtype='i8')
  masstab = numpy.zeros(6, dtype='f8')

  for fname, header in H:
    Ntot += header['N']
    mask = header['mass'] != 0
    masstab[mask] = header['mass'][mask]
  print 'total particles written', Ntot
  for fname, header in H:
    header['mass'] = masstab
    header['Ntot_low'][:] = (Ntot & 0xffffffff)
    header['Ntot_high'][:] = (Ntot >> 32)
    header['Nfiles'] = len(H)
  
    print fname, header['N']
    rewrite_header(fname, header)

def setup_level(ilevel):
  Nmesh = args.Levels[ilevel]
  rt = lambda : None
  if args.Scale[Nmesh] == 0: 
    rt.data = setup_baselevel(Nmesh)
  else:
    rt.data = setup_innerlevel(Nmesh)
  rt.scale = args.Scale[Nmesh]
  rt.Nmesh = Nmesh
  rt.DownSample = args.DownSample[Nmesh]
  rt.makegas = args.MakeGas[Nmesh]
  rt.ilevel = ilevel
  rt.ptype = args.DMptype[Nmesh]
  rt.usemasstab = numpy.sum([args.DMptype[n] == rt.ptype for n in args.Levels]) == 1
  return rt

def setup_baselevel(Nmesh):
  data = numpy.empty((Nmesh, Nmesh, Nmesh), dtype=[
           ('gps', ('i4', 3)), 
           ('disp', ('f4', 3)),
           ('delta', 'f4')])

  grid1d = numpy.arange(Nmesh);

  data['gps'][..., 2] = grid1d[None, None, :]
  data['gps'][..., 1] = grid1d[None, :, None]
  data['gps'][..., 0] = grid1d[:, None, None]

  data = data.reshape(-1)
  for i, block in enumerate(['dispx', 'dispy', 'dispz']):
    content = readblock(Nmesh, block, 'f4', args)
    data['disp'][..., i] = content

  data['delta'] = readblock(Nmesh, 'delta', 'f4', args)

  return data

def setup_innerlevel(Nmesh):
  index = readblock(Nmesh, 'index', 'i8', args)
  data = numpy.empty(len(index), dtype=[
         ('gps', ('i4', 3)), 
         ('disp', ('f4', 3)),
         ('delta', 'f4'),
         ('region', 'u4')
         ])
  data['region'][:] = readblock(Nmesh, 'region', 'u4', args)
  data['disp'][:, 0] = readblock(Nmesh, 'dispx', 'f4', args)
  data['disp'][:, 1] = readblock(Nmesh, 'dispy', 'f4', args)
  data['disp'][:, 2] = readblock(Nmesh, 'dispz', 'f4', args)
  data['delta'][:] = readblock(Nmesh, 'delta', 'f4', args)
  Size = numpy.array([args.Meta[Nmesh]['Size'][i] for i in range(len(args.Regions))])
  Stride1 = Size[:, 1] * Size[:, 2]
  Stride2 = Size[:, 2]
  Offset = numpy.array([args.Meta[Nmesh]['Offset'][i] for i in range(len(args.Regions))])
  data['gps'][:, 0] = index // Stride1[data['region']]
  data['gps'][:, 1] = (index \
           - data['gps'][:, 0] * Stride1[data['region']]) \
            // Stride2[data['region']]
  data['gps'][:, 2] = index % Stride2[data['region']]
  data['gps'] += Offset[data['region']]
  # now gonna remove the overlaps and make sure points are within the mesh.
  numpy.remainder(data['gps'], Nmesh, data['gps'])
  index = numpy.ravel_multi_index(data['gps'].T, (Nmesh, Nmesh, Nmesh))
  u, ui = numpy.unique(index, return_index=True)
  return data[ui]

def dig(this, next):
  # first select active mesh from next level, aligned to this level
  next.data = next.data[activemask(next, this.Nmesh)]
  # calculate the index of the covered mesh on this level
  index_from_next = numpy.ravel_multi_index(
              next.data['gps'].T // (next.Nmesh // this.Nmesh),
              (this.Nmesh, this.Nmesh, this.Nmesh))
  # only use unique ones
  index_from_next, label = numpy.unique(index_from_next, return_inverse=True)
  print index_from_next, label
  # check mass conservation 
  assert (numpy.bincount(label) == (next.Nmesh // this.Nmesh) ** 3).all()
  if next.DownSample > 1:
    # downsample levels need to carry on the interpolation
    # of previous level
    next.data['delta'] += trilinearinterp(next.data['gps'] * (1.0 * this.Nmesh / next.Nmesh),
                    this.data['gps'], 
                    this.data['delta'], 
                   (next.Nmesh, next.Nmesh, next.Nmesh))
    next.data['disp'] += trilinearinterp(next.data['gps'] * (1.0 * this.Nmesh / next.Nmesh),
                    this.data['gps'], 
                    this.data['disp'], 
                   (next.Nmesh, next.Nmesh, next.Nmesh))
  # index this level
  index = numpy.ravel_multi_index(this.data['gps'].T, (this.Nmesh, this.Nmesh, this.Nmesh))
  # locate those to be removed
  ind = index_from_next.searchsorted(index)
  covered_mask = (index_from_next.take(ind, mode='clip') == index)

  assert covered_mask.sum() == len(index_from_next)
  preserve_mask = ~covered_mask
  this.data = this.data[preserve_mask]

def activemask(level, align):
  sep = args.BoxSize / align
  regions = args.Regions
  gridcenter = regions['center'] / sep
  halfsize = 1. / (regions['size'] * level.scale * 0.5)
  assert level.Nmesh % align == 0
  dis = level.data['gps'] // (level.Nmesh // align) - gridcenter[level.data['region']]
  dis[dis < -0.5] += 1
  dis *= sep
  dis += args.BoxSize * 0.5
  numpy.remainder(dis, args.BoxSize, dis)
  dis -= args.BoxSize * 0.5
  dis *= halfsize[level.data['region']]
  if not args.cubic:
    dis **= 2
    dis = dis.sum(axis=-1)
    mask = dis < 1.0
  else:
    numpy.abs(dis, dis)
    mask = (dis < 1.0).all(axis=-1)
  return mask

def writerecord(file, data):
    foo = numpy.empty(1, 'i4')
    foo[...] = data.nbytes
    foo.tofile(file)
    data.tofile(file)
    foo.tofile(file)

def rewrite_header(fname, header):
  with file(fname, 'r+') as icfile:
    writerecord(icfile, header)
  
def write_gadget(level, FID, IDstart):
  nchunks = int(numpy.ceil(len(level.data) / (1.0 * args.Npar)))
  if nchunks < 1: nchunks = 1
  datasp = numpy.array_split(level.data, nchunks)
  H = []
  for k in range(nchunks):
    fname = '%s.%d' % (args.prefix, FID + k)
    header = write_gadget_one(fname, level, datasp[k], IDstart)
    H.append((fname, header))
    IDstart += header['N'].sum()
  return H, IDstart

def write_gadget_one(fname, level, data, IDStart=None):
  Nmesh = level.Nmesh
  with file(fname, 'w') as icfile:
    critical_mass = 3 * args.H ** 2 / (8 * numpy.pi * args.G) * (args.BoxSize / Nmesh) ** 3

    header = numpy.zeros(None, GHEADER)
    header['time'] = args.a
    header['redshift'] = 1. / args.a - 1
    header['boxsize'] = args.BoxSize
    header['OmegaM'] = args.OmegaM
    header['OmegaL'] = args.OmegaL
    header['h'] = args.h
    Npar = len(data)

    header['N'][level.ptype] = Npar
    if level.makegas:
      header['N'][0] = Npar
      if level.usemasstab:
        header['mass'][level.ptype] = (args.OmegaM - args.OmegaB) * critical_mass
    else:
      if level.usemasstab:
        header['mass'][level.ptype] = args.OmegaM * critical_mass

    writerecord(icfile, header)

    # position
    pos = (numpy.float32(data['gps']) + 0.5) * (1.0 * args.BoxSize / Nmesh)
    if not args.nodisp:
      pos += data['disp']
    if level.makegas: 
      pos = numpy.tile(pos, (2, 1))
      # gas
      pos[:Npar] -= 0.5 * (1 - args.OmegaB / args.OmegaM) * args.BoxSize / Nmesh
      # dm
      pos[Npar:] += 0.5 * args.OmegaB / args.OmegaM * args.BoxSize / Nmesh

    numpy.remainder(pos, args.BoxSize, pos)

    writerecord(icfile, pos)
    pos = None

    # velocity
    vel = data['disp'] * args.Meta[Nmesh]['vfact']
    if level.makegas: 
      vel = numpy.tile(vel, (2, 1))
    writerecord(icfile, vel)
    vel = None

    # id
    # first four bits are left for gadget
    # the coming bit is gas or dm.
    # the coming 7 bits is ilevel.
    # then is the grid position.
    # this system supports a mesh size of 512K.
    if level.makegas:
      id = numpy.empty(Npar + Npar, dtype=args.iddtype)
    else:
      id = numpy.empty(Npar, dtype=args.iddtype)

    bits = args.iddtype.itemsize * 8
    if args.idscheme == 'mesh':
      id[:Npar] = numpy.ravel_multi_index(data['gps'].T, (Nmesh, Nmesh, Nmesh))
      if level.makegas:
        id[Npar:] = id[:Npar]
        id[:Npar] |= args.iddtype.type(1L << (bits - 5))

      id[:] |= args.iddtype.type(level.ilevel << (bits - 12))
    else:
      id[:] = numpy.arange(IDstart, len(id) + IDstart)


    writerecord(icfile, id)
    id = None

    # mass
    if level.usemasstab:
      if level.makegas:
        mass = numpy.empty(Npar, dtype='f4')
        mass[:] = args.OmegaB * critical_mass
        writerecord(icfile, mass)
        mass = None 
    else:
      if level.makegas:
        mass = numpy.empty(Npar + Npar, dtype='f4')
        mass[:Npar] = args.OmegaB * critical_mass
        mass[Npar:] = (args.OmegaM - args.OmegaB) * critical_mass
      else:
        mass = numpy.empty(Npar, dtype='f4')
        mass[:Npar] = args.OmegaM * critical_mass
      writerecord(icfile, mass)
      mass = None 
    # ie ( all zeros)
    if level.makegas:
      ie = numpy.zeros(Npar, dtype='f4')
      writerecord(icfile, ie)
      ie = None
  return header

