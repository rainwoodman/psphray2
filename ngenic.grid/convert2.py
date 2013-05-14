import numpy
from scipy.ndimage import map_coordinates
# this requires numpy > 1.6.0
from args import parseargs
from myio import readchunk, writechunk, yieldfilename, write_gadget_one, rewrite_header, writerecord

def joinchunks(h, exbot=None, extop=None):
  overlay_list = []
  for fn in list(yieldfilename(h)):
    chunk = readchunk(h, fn, exbot, extop)
    ipos = chunk['ipos']
    if len(ipos) > 0:
      print 'ipos stat', 'min', ipos.min(axis=0), 'max', ipos.max(axis=0)
    else:
      print 'ipos empty'
    overlay_list.append(chunk)
  overlay = numpy.concatenate(overlay_list)
  overlay_list = []
  index = numpy.ravel_multi_index(overlay['ipos'].T - h['Offset'][:, None], h['Size'])
  overlay = overlay[index.argsort()]
  return overlay

def interpolate(h, hnext):
  """ interpolate h into hnext """
  # round down
  exbot = (hnext['Offset'] * h['Nmesh']) // hnext['Nmesh']
  # round up
  extop = ((hnext['Offset'] + hnext['Size']) * h['Nmesh'] + h['Nmesh'] - 1)// hnext['Nmesh']

  print h['Nmesh'], 'box', exbot, extop

  overlay = joinchunks(h, exbot, extop)
  src = overlay.reshape(*(extop-exbot))

  for fn in list(yieldfilename(hnext, lookfor='delta1')):
    # read in the current level displacement
    chunk = readchunk(hnext, fn, None, None, extra='1')
    ipos = chunk['ipos']
    ipos = 1.0 * ipos * h['Nmesh'] / hnext['Nmesh'] - exbot
    chunk['delta'] += map_coordinates(src['delta'], ipos.T)
    for d in range(3):
      chunk['disp'][:, d] += map_coordinates(src['disp'][..., d], ipos.T)
    # write to the full displacement
    writechunk(hnext, fn, chunk)

def main():
  A = parseargs()
  if A.action == 'interpolate':
    main_interpolate(A)
  if A.action == 'gadget':
    main_gadget(A)

def main_interpolate(A):
  Nmesh = A.L['Nmesh']
  for i in range(1, len(A.L)):
    h = A.META[Nmesh[i - 1]]
    hnext = A.META[Nmesh[i]]
    if hnext['DownSample'] > 1:
      interpolate(h, hnext)

def main_gadget(A):
  Nmesh = A.L['Nmesh']
  fid = 0
  F = []
  for i in range(len(A.L)):
    h = A.META[Nmesh[i]]

    overlay = joinchunks(h).reshape(h['Size'])
    for fn in list(yieldfilename(h)):
      chunk = readchunk(h, fn, None, None)
      if i != 0:
        IncludeCenter = h['ZoomCenter']
        IncludeRadius = h['ZoomRadius']
        includemask = selectmask(chunk, 
               A.BoxSize, 
               Nmesh[i],
               Nmesh[i - 1],
               IncludeCenter, IncludeRadius, A.cubic)
        chunk = chunk[includemask]
        print 'using', includemask.sum(), 'particles'
        ipos = chunk['ipos']
        print 'ipos stat', 'min', ipos.min(axis=0), 'max', ipos.max(axis=0)
      if i != len(A.L) - 1:
        hnext = A.META[Nmesh[i + 1]]
        ExcludeCenter = hnext['ZoomCenter']
        ExcludeRadius = hnext['ZoomRadius']

        excludemask = selectmask(chunk,
               A.BoxSize,
               Nmesh[i], 
               Nmesh[i],
               ExcludeCenter, ExcludeRadius, A.cubic)
        exclude = chunk[excludemask]
        chunk = chunk[~excludemask]
        print 'excluding', excludemask.sum(), 'particles'
        ipos = exclude['ipos']
        if len(ipos):
          print 'exclude stat', 'min', ipos.min(axis=0), 'max', ipos.max(axis=0)

      # here we use the center cell value for these particles.
      # the locations are shifted later in write_gadget_one.
      chunk['delta'] = map_coordinates(overlay['delta'], 
              chunk['ipos'].T + 0.5 - h['Offset'][:, None], mode='wrap')
      for d in range(3):
          chunk['disp'][:, d] = map_coordinates(overlay['disp'][..., d], 
              chunk['ipos'].T + 0.5 - h['Offset'][:, None], mode='wrap')

      for start in range(0, len(chunk), A.Npar):
        F.append(write_gadget_one('%s.%d' % (A.prefix, fid), h, chunk[start:start+A.Npar]))
        print F[-1][1]['mass']
        fid = fid + 1

  Ntot = numpy.zeros(6, dtype='i8')
  mass = numpy.zeros(6, dtype='f8')

  Ntot[:] = numpy.sum([ header['N'] for fname, header in F], axis=0)
  mass[:] = [numpy.unique([header['mass'][ptype] for fname, header in F])[-1] for ptype in range(6)]
  print [numpy.unique([header['mass'][ptype] for fname, header in F]) for ptype in range(6)]
  print 'Ntot', Ntot
  print 'mass', mass

  for fname, header in F:
    header['mass'] = mass
    header['Ntot_low'][:] = (Ntot & 0xffffffff)
    header['Ntot_high'][:] = (Ntot >> 32)
    header['Nfiles'] = len(F)
  
    print 'touching up',  fname
    rewrite_header(fname, header)

def selectmask(chunk, BoxSize, Nmesh, Nalign, center, radius, cubic):
  assert Nmesh % Nalign == 0
  AlignSpacing = BoxSize / Nalign
  gridcenter = center / AlignSpacing
  print Nmesh, Nalign, AlignSpacing, gridcenter
  halfsize = 1. / radius
  dis = (chunk['ipos'] // (Nmesh // Nalign) + 0.5) - gridcenter
  dis *= AlignSpacing
  dis += BoxSize * 0.5
  numpy.remainder(dis, BoxSize, dis)
  dis -= BoxSize * 0.5
  dis *= halfsize
  if not cubic:
    dis **= 2
    dis = dis.sum(axis=-1)
    mask = dis < 1.0
  else:
    numpy.abs(dis, dis)
    mask = (dis < 1.0).all(axis=-1)
  return mask
main()

def main_ramses(A):
  def write_ramses(header, Nmesh, name, buffer):
    with file('%s%d-%s' % (A.prefix, Nmesh, name), 'w') as icfile:
      buffer = buffer.reshape(header['np'])
      writerecord(icfile, header)
      for row in buffer.transpose((2, 0, 1)):
        writerecord(icfile, row.ravel('F'))

  RHEADER = numpy.dtype([
   ('np', ('i4', 3)), ('dx', 'f4'), 
   ('xo', ('f4', 3)), ('astart', 'f4'), 
   ('omegam', 'f4'), ('omegav', 'f4'), ('h0', 'f4')])
  # to MPC
  fac = 1.0 / A.h * args.UnitLength_in_cm / 3.08567758e24
  boxsize = A.BoxSize * fac

  Nmesh = A.L['Nmesh']
  fid = 0
  F = []
  for i in range(len(A.L)):
    h = A.META[Nmesh[i]]

    header['dx'] = boxsize / Nmesh[i]
    header['astart'] = A.a
    header['omegam'] = A.OmegaM
    header['omegav'] = A.OmegaL
    header['h0'] = A.h * 100.0
    if i == 0:
      header['np'][:] = Nmesh
      header['xo'][:] = 0
    else:
      header['np'][:] = h['Size']
      header['xo'][:] = header['dx'] * h['Offset']

    overlay = joinchunks(h)
    write_ramses(header, Nmesh, 'deltab', overlay['delta'])
    for i, block in enumerate(['dispx', 'dispy', 'dispz']):
      disp = overlay['disp'][:, i] * fac
      write_ramses(header, Nmesh, block, disp)
      disp = overlay['disp'][:, i]
      # is ramses vel unit in km/s?
      disp * (h['vfact'] * A.a ** 0.5)
      write_ramses(header, Nmesh, ['velcx', 'velcy', 'velcz'][i], disp)
