import numpy
import os.path

def yieldfilename(H):
  width2 = len(str(H['NTask']))
  if H['DownSample'] > 1:
    width1 = len(str(H['DownSample']))
    for d in range(H['DownSample']):
      for fid in range(H['NTask']):
        fname = '%s/%%s-%d.%0*d-%0*d' % (
            H['datadir'], H['Nmesh'], width1, d, width2, fid)
        if os.path.exists(fname % 'delta'):
          yield fname
  else:
    for fid in range(H['NTask']):
      fname = '%s/%%s-%d.%0*d' % (
          H['datadir'], H['Nmesh'], width2, fid)
      if os.path.exists(fname % 'delta'):
        yield fname

def readchunk(H, fname, bottom=None, top=None, extra=''):
  fid = int(fname.replace('.', '-').split('-')[-1])
  dispx = numpy.memmap(fname % ('dispx' + extra), mode='r', dtype='f4')
  dispy = numpy.memmap(fname % ('dispy' + extra), mode='r', dtype='f4')
  dispz = numpy.memmap(fname % ('dispz' + extra), mode='r', dtype='f4')
  delta = numpy.memmap(fname % ('delta' + extra), mode='r', dtype='f4')
  if H['Scale'] == 0.0:
    assert H['DownSample'] == 1
    xstart = fid * H['Nmesh'] // H['NTask']
    xend = (fid + 1) * H['Nmesh'] // H['NTask']
    index = numpy.arange(
               xstart * H['Nmesh'] * H['Nmesh'],
               xend * H['Nmesh'] * H['Nmesh'])
  else:
    index = numpy.memmap(fname % 'index', mode='r', dtype='i8')
  ipos = numpy.array(numpy.unravel_index(index, H['Size']), dtype='i4').T
  ipos += H['Offset']
  if bottom is not None and top is not None:
    includemask = ipos_in_box(ipos, bottom, top)
  else:
    includemask = numpy.ones(len(ipos), dtype='?')
  result = numpy.empty(includemask.sum(), dtype=[('ipos', ('i4', 3)), ('disp', ('f4', 3)), ('delta', 'f4')])
  result['ipos'] = ipos[includemask]
  result['disp'][:, 0] = dispx[includemask]
  result['disp'][:, 1] = dispy[includemask]
  result['disp'][:, 2] = dispz[includemask]
  result['delta'] = delta[includemask]
  return result

def writechunk(H, fname, chunk):
  assert H['DownSample'] > 1
  # we shall never write anything to a non downsampling level
  chunk['disp'][:, 0].tofile(fname % 'dispx')
  chunk['disp'][:, 1].tofile(fname % 'dispy')
  chunk['disp'][:, 2].tofile(fname % 'dispz')
  chunk['delta'].tofile(fname % 'delta')

def ipos_in_box(ipos, x1, x2):
  return ((ipos >= x1[None, :]) & \
          (ipos < x2[None, :])).all(axis=-1)

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

def write_gadget_one(fname, h, chunk):
  A = h['A']
  mass_c = 3 * A.H ** 2 / (8 * numpy.pi * A.G) * \
                 (1.0 * A.BoxSize / h['Nmesh']) ** 3
  N = len(chunk)

  dmptype = h['ptype']
  if h['makegas']:
    gasptype = 0
    dmmass = mass_c * (A.OmegaM - A.OmegaB)
    gasmass = mass_c * A.OmegaB
  else:
    gasptype = None
    dmmass = mass_c * A.OmegaM

  Spacing = h['A'].BoxSize / h['Nmesh']

  with file(fname, 'w') as icfile:
    header = numpy.zeros(None, GHEADER)
    header['time'] = A.a
    header['redshift'] = 1. / A.a - 1
    header['boxsize'] = A.BoxSize
    header['OmegaM'] = A.OmegaM
    header['OmegaL'] = A.OmegaL
    header['h'] = A.h
    header['N'][dmptype] = N
    header['flag_double'] = 0
    if gasptype is not None:
      header['N'][gasptype] = N

    writerecord(icfile, header)

    # position
    pos = (chunk['ipos'] + 0.5) * Spacing
    if not A.nodisp:
      pos += chunk['disp']

    pos = numpy.float32(pos)
    if h['makegas']: 
      pos = numpy.tile(pos, (2, 1))
      # gas
      pos[:N] -= 0.5 * (1 - A.OmegaB / A.OmegaM) * Spacing
      # dm
      pos[N:] += 0.5 * A.OmegaB / A.OmegaM * Spacing

    numpy.remainder(pos, A.BoxSize, pos)

    writerecord(icfile, pos)
    pos = None

    # velocity
    vel = numpy.float32(chunk['disp'] * h['vfact'])
    if h['makegas']: 
      vel = numpy.tile(vel, (2, 1))
    writerecord(icfile, vel)
    vel = None

    # id
    # first four bits are left for gadget
    # the coming bit is gas or dm.
    # the coming 7 bits is ilevel.
    # then is the grid position.
    # this system supports a mesh size of 512K.
    id = numpy.uint64(numpy.ravel_multi_index(chunk['ipos'].T, (h['Nmesh'],) * 3))
    id += h['ilevel'] << 55
    if h['makegas']:
      id = numpy.tile(id, (2, 1))
      id[N:] += 1 << 54
    writerecord(icfile, id)
    id = None

    # mass
    if h['makegas']:
      mass = numpy.empty(N * 2, dtype='f4')
      mass[:N] = gasmass
      mass[N:] = dmmass
    else:
      mass = numpy.empty(N, dtype='f4')
      mass[:] = dmmass
    writerecord(icfile, mass)
    mass = None 

    # ie ( all zeros)
    if h['makegas']:
      ie = numpy.zeros(N, dtype='f4')
      writerecord(icfile, ie)
      ie = None
  return fname, header

def rewrite_header(fname, header):
  with file(fname, 'r+') as icfile:
    writerecord(icfile, header)

def writerecord(file, data):
    # do not write empty
    if data.nbytes == 0: return
    foo = numpy.empty(1, 'i4')
    foo[...] = data.nbytes
    foo.tofile(file)
    data.tofile(file)
    foo.tofile(file)
