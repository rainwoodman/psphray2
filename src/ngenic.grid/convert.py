import ConfigParser
import argparse
import numpy

parser = argparse.ArgumentParser()
parser.add_argument("paramfile", help="use the -used file generated from genic")
parser.add_argument("format", choices=['gadget', 'ramses'])
parser.add_argument("-o", "--prefix", help="outputprefix of the IC files", default="IC")
parser.add_argument("-q", '--cubic', action='store_true', default=False)

def parseargs(parser):
  args = parser.parse_args()
  config = ConfigParser.ConfigParser()
  config.read(args.paramfile)

  args.datadir = config.get("IO", "datadir")
  args.BoxSize = config.getfloat("IC", "BoxSize")
  args.a = config.getfloat("IC", "a")
  args.OmegaM = config.getfloat("Cosmology", "OmegaM")
  args.OmegaB = config.getfloat("Cosmology", "OmegaB")
  args.OmegaL = config.getfloat("Cosmology", "OmegaL")
  args.h = config.getfloat("Cosmology", "h")
  args.UnitLength_in_cm = config.getfloat("Unit", "UnitLength_in_cm")
  args.UnitVelocity_in_cm_per_s = config.getfloat("Unit", "UnitVelocity_in_cm_per_s")
  args.UnitMass_in_g = config.getfloat("Unit", "UnitMass_in_g")
  args.UnitTime_in_s = args.UnitLength_in_cm / args.UnitVelocity_in_cm_per_s;
  args.HUBBLE = 3.2407789e-18
  args.GRAVITY = 6.672e-8
  args.G = args.GRAVITY / args.UnitLength_in_cm ** 3 * args.UnitMass_in_g * args.UnitTime_in_s ** 2
  args.Hubble = args.HUBBLE * args.UnitTime_in_s

  Regions = []
  for name, value in config.items("Regions"):
    all = [float(t) for t in value.replace(';', ',').split(',')]
    if len(all) == 4:
      x, y, z, s = all
      sx = s
      sy = s
      sz = s
    elif len(all) == 6:
      x, y, z, sx, sy, sz = all
  
    Regions.append(((x, y, z), (sx, sy, sz)))
  args.Regions = numpy.array(Regions, dtype=[('center', ('f8', 3)), ('size', ('f8', 3))])
  
  Levels = []
  
  args.Scale = {}
  args.Meta = {}
  for name, value in config.items("Levels"):
     n, s = [float(t) for t in value.replace(';', ',').split(',')]
     Levels.append(n)
     args.Scale[n] = s
     args.Meta[n] = read_meta(n, args.datadir)

  args.Levels = numpy.array(Levels, dtype='i8')
  args.Levels.sort()

  return args

def read_meta(Nmesh, datadir):
    meta = {'Offset':{}, 'Size':{}}
    metalines = file('%s/filters-%d' % (datadir, Nmesh)).read()
    exec(metalines, meta)
    return meta


args = parseargs(parser)

def setup_level(Nmesh):
  if args.Scale[Nmesh] == 0: 
    return setup_baselevel(Nmesh)
  else:
    return setup_innerlevel(Nmesh)
  
def readblock(Nmesh, block, dtype):
  meta = args.Meta[Nmesh]
  NTask = meta['NTask']
  ZeroPadding = int(numpy.ceil(numpy.log10(NTask+0.1)))
  content = []
  for fid in range(NTask):
    content.append(numpy.fromfile('%s/%s-%d.%0*d' %(
                        args.datadir, 
                        block, Nmesh, ZeroPadding, fid), dtype))
  return numpy.concatenate(content)

def setup_innerlevel(Nmesh):
  index = readblock(Nmesh, 'index', 'i8')
  r = readblock(Nmesh, 'region', 'u1')
  dispx = readblock(Nmesh, 'dispx', 'f4')
  dispy = readblock(Nmesh, 'dispy', 'f4')
  dispz = readblock(Nmesh, 'dispz', 'f4')
  delta = readblock(Nmesh, 'delta', 'f4')
  meta = args.Meta[Nmesh]

  alldata = []
  for i, region in enumerate(args.Regions):
    mask = r == i
    data = numpy.empty(mask.sum(), dtype=[
           ('gps', ('i4', 3)), 
           ('disp', ('f4', 3)),
           ('delta', 'f4')])
    data['gps'][:, 0], data['gps'][:, 1], data['gps'][:, 2] = numpy.unravel_index(index[mask], meta['Size'][i])
    data['gps'] += meta['Offset'][i]
    data['disp'][:, 0] = dispx[mask]
    data['disp'][:, 1] = dispy[mask]
    data['disp'][:, 2] = dispz[mask]
    data['delta'] = delta[mask]
    alldata.append(data)
  data = numpy.concatenate(alldata)
  # now gonna remove the overlaps and make sure points are within the mesh.
  numpy.remainder(data['gps'], Nmesh, data['gps'])
  index = numpy.ravel_multi_index(data['gps'].T, (Nmesh, Nmesh, Nmesh))
  u, ui = numpy.unique(index, return_index=True)
  return data[ui]

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
    content = readblock(Nmesh, block, 'f4')
    data['disp'][..., i] = content

  data['delta'] = readblock(Nmesh, 'delta', 'f4')

  return data

def ramses():
  pass

def gadget():
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

  def dig(data, Nmesh, scale, Align=None):
    """ dig a hole on Nmesh points data, for NmeshDig
        returns the mask of the hole,
        align to Align, if omitted use Nmesh
    """
    if Align is None: Align = Nmesh
    sep = args.BoxSize / Align
    assert Nmesh % Align == 0
    mask2 = numpy.zeros(len(data), dtype='?')
    for region in args.Regions:
      gridcenter = region['center'] / sep
      dis = data['gps'] // (Nmesh // Align) - gridcenter
      # use the if any part of the mesh overlaps the sphere,
      dis[dis < - 0.5] += 1
      dis *= sep
      dis += args.BoxSize * 0.5
      numpy.remainder(dis, args.BoxSize, dis)
      dis -= args.BoxSize * 0.5
      dis /= (0.5 * region['size'] * scale)
      if not args.cubic:
        dis **= 2
        dis = dis.sum(axis=-1)
        mask2 |= (dis < 1.0)
      else:
        numpy.abs(dis, dis)
        mask2 |= (dis < 1.0).all(axis=-1)
    return mask2

  def writerecord(file, data):
    foo = numpy.empty(1, 'i4')
    foo[...] = data.nbytes
    foo.tofile(file)
    data.tofile(file)
    foo.tofile(file)

  for i, Nmesh in enumerate(args.Levels):
    data = setup_level(Nmesh)
    icfile = file('%s.%d' % (args.prefix, i), 'w')
    if i < len(args.Levels) - 1:
      # if there is a refine level, dig a hole from the current level
      # leaving space for the next level.
      NmeshNext = args.Levels[i + 1]
      assert NmeshNext % Nmesh == 0
      holemask = dig(data, Nmesh, scale=args.Scale[NmeshNext])
      data = data[~holemask]
    if i > 0:
      # if there is a parent level, clip the current level
      # to fill in the space of the parent level
      NmeshPrev = args.Levels[i - 1]
      assert Nmesh % NmeshPrev == 0
      fillmask = dig(data, Nmesh, scale=args.Scale[Nmesh], Align=NmeshPrev)
      data = data[fillmask]

    header = numpy.zeros(None, GHEADER)
    header['time'] = args.a
    header['redshift'] = 1. / args.a + 1
    header['boxsize'] = args.BoxSize
    header['OmegaM'] = args.OmegaM
    header['OmegaL'] = args.OmegaL
    header['h'] = args.h
    header['N'][0] = len(data)
    header['Nfiles'] = len(args.Levels)
    header['Ntot_low'][0] = len(data)
    writerecord(icfile, header)

    pos = data['gps'] * args.BoxSize
    pos /= Nmesh
    pos += data['disp']
    numpy.remainder(pos, args.BoxSize, pos)
    writerecord(icfile, pos)
    pos = None
    vel = data['disp'] * args.Meta[Nmesh]['vfact']
    writerecord(icfile, vel)
    vel = None
    id = numpy.empty(len(data), dtype=('i4', 2))
    id[:, 1] = Nmesh
    id[:, 0] = numpy.arange(len(data))
    writerecord(icfile, id)
    id = None
    mass = numpy.empty(len(data), dtype='f4')
    critical_mass = 3 * args.Hubble ** 2 / (8 * numpy.pi * args.G) * (args.BoxSize / Nmesh) ** 3
    mass[:] = args.OmegaM * critical_mass
    writerecord(icfile, mass)
    mass[:] = 0
    writerecord(icfile, mass)

def main():
  if args.format == 'gadget':
    gadget()
  elif args.format == 'ramses':
    ramses()

main()
